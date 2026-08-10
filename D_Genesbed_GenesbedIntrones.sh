#!/usr/bin/env bash
#
# generar_bed_intrones.sh
#
# Genera un BED de INTRONES (con margen/padding hacia las regiones
# flanqueantes de splicing) para los genes de interés, a partir de:
#
#   1. genes_CMT.bed          -> BED de genes completos (script genes_to_bed.py)
#   2. GTF de MANE Select     -> anotación de exones por transcrito canónico
#
# Flujo:
#   GTF MANE -> exones por gen (restringidos a genes_CMT.bed)
#            -> exones_mane.bed (ordenado, fusionado)
#   genes_CMT.bed - exones_mane.bed = intrones "crudos" (bedtools subtract)
#   intrones "crudos" + padding hacia el exón (bedtools slop, recortado
#   para no invadir el exón) = intrones_mane_padded.bed
#
# Requisitos: bedtools, awk, sort, un archivo .genome (tamaños de cromosoma)
#             para bedtools slop.
#
# Descarga previa necesaria (una sola vez, no por muestra):
#   GTF de MANE: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#                MANE.GRCh38.vX.X.ensembl_genomic.gtf.gz
#   Tamaños de cromosoma GRCh38 (.genome / .fai):
#                samtools faidx GRCh38.fa && cut -f1,2 GRCh38.fa.fai > GRCh38.genome

set -euo pipefail

# ---------------------------------------------------------------------------
# PARSEO DE FLAGS con getopts (opciones cortas; getopts no soporta --flags
# largos de forma nativa en bash)
# ---------------------------------------------------------------------------
mostrar_ayuda() {
    cat <<EOF
Uso: $(basename "$0") -g FILE -m FILE -f FILE [-p N] [-o DIR]

Obligatorios:
  -g FILE   BED de genes completos (ej. genes_CMT.bed)
  -m FILE   GTF de MANE Select (GRCh38)
  -f FILE   Archivo .genome (chrom<TAB>size), ver samtools faidx + cut -f1,2

Opcionales:
  -p N      Pares de bases de margen hacia el exón (default: 15)
  -o DIR    Carpeta de salida (default: ./bed_intrones)
  -h        Muestra esta ayuda

También se pueden definir como variables de entorno (GENES_BED, MANE_GTF,
GENOME_FILE, PADDING, OUTDIR); los flags, si se indican, tienen prioridad.

Ejemplo:
  $(basename "$0") \\
      -g genes_CMT.bed \\
      -m MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz \\
      -f GRCh38.genome \\
      -p 15 \\
      -o resultados/bed_intrones
EOF
}

while getopts ":g:m:f:p:o:h" opt; do
    case "$opt" in
        g) GENES_BED="$OPTARG" ;;
        m) MANE_GTF="$OPTARG" ;;
        f) GENOME_FILE="$OPTARG" ;;
        p) PADDING="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) mostrar_ayuda; exit 0 ;;
        \?) echo "ERROR: opción no reconocida: -$OPTARG" >&2; mostrar_ayuda; exit 1 ;;
        :)  echo "ERROR: la opción -$OPTARG requiere un valor" >&2; mostrar_ayuda; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# ---------------------------------------------------------------------------
# CONFIGURACIÓN (valores por flag, si no, por variable de entorno; falla si
# ninguno de los dos está definido en los campos obligatorios)
# ---------------------------------------------------------------------------
GENES_BED="${GENES_BED:?Debe definir GENES_BED (-g o variable de entorno)}"
MANE_GTF="${MANE_GTF:?Debe definir MANE_GTF (-m o variable de entorno)}"
GENOME_FILE="${GENOME_FILE:?Debe definir GENOME_FILE (-f o variable de entorno)}"
PADDING="${PADDING:-15}"
OUTDIR="${OUTDIR:-$(pwd)/bed_intrones}"

mkdir -p "$OUTDIR" "$OUTDIR/logs"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

EXONES_RAW="$OUTDIR/exones_mane_raw.bed"
EXONES_BED="$OUTDIR/exones_mane.bed"
INTRONES_CRUDOS="$OUTDIR/intrones_mane_raw.bed"
INTRONES_PADDED="$OUTDIR/intrones_mane_padded.bed"

# ---------------------------------------------------------------------------
# 1. Extraer nombres de gen desde el BED de genes (columna 4)
# ---------------------------------------------------------------------------
GENES_LISTA="$OUTDIR/lista_genes.txt"
cut -f4 "$GENES_BED" | sort -u > "$GENES_LISTA"
N_GENES=$(wc -l < "$GENES_LISTA")
log "Paso 1/5: ${N_GENES} genes leídos desde $GENES_BED"

# ---------------------------------------------------------------------------
# 2. Extraer exones del GTF de MANE, restringidos a esos genes
# ---------------------------------------------------------------------------
log "Paso 2/5: extrayendo exones del GTF de MANE para los genes de interés"

# Detecta si el GTF está comprimido (gzip) y lo descomprime al vuelo.
# awk no interpreta gzip directamente: si se le pasa un .gtf.gz sin
# descomprimir, no encuentra ninguna línea "exon" y el script falla más
# adelante con un error de "no se encontraron exones".
if [[ "$MANE_GTF" == *.gz ]]; then
    GTF_READER="zcat"
else
    GTF_READER="cat"
fi

"$GTF_READER" "$MANE_GTF" \
    | awk -F'\t' '$3=="exon"' \
    | grep -Ff <(sed 's/^/gene_name "/; s/$/"/' "$GENES_LISTA") \
    > "$OUTDIR/exones_mane_gtf.tmp" \
    2> "$OUTDIR/logs/02_grep_exones.log"

if [[ ! -s "$OUTDIR/exones_mane_gtf.tmp" ]]; then
    log "  [ERROR] No se encontraron exones para los genes indicados."
    log "  Comprobar: (1) que MANE_GTF no esté vacío/corrupto, (2) que use el"
    log "  atributo 'gene_name \"SIMBOLO\"' en la columna 9, (3) que los símbolos"
    log "  de $GENES_LISTA coincidan exactamente (mayúsculas/minúsculas) con el GTF."
    log "  Prueba manual: zcat \"$MANE_GTF\" | awk -F'\\t' '\$3==\"exon\"' | head -1"
    exit 1
fi

# Convertir a BED (GTF es 1-based inclusive; BED es 0-based semi-abierto)
# Nota: la extracción del gene_name usa split()/array de awk (funciones
# POSIX estándar, compatibles con mawk/gawk/awk BSD), evitando tanto el
# match() de 3 argumentos (específico de gawk) como un sed frágil.
awk -F'\t' 'BEGIN{OFS="\t"}
    {
        chrom=$1
        if (chrom !~ /^chr/) chrom="chr"chrom
        start=$4-1
        end=$5
        strand=$7

        gene=""
        n = split($9, attrs, ";")
        for (i=1; i<=n; i++) {
            if (attrs[i] ~ /gene_name/) {
                split(attrs[i], kv, "\"")
                gene = kv[2]
            }
        }
        print chrom, start, end, gene, ".", strand
    }' "$OUTDIR/exones_mane_gtf.tmp" > "$EXONES_RAW"

rm -f "$OUTDIR/exones_mane_gtf.tmp"

N_EXONES_RAW=$(wc -l < "$EXONES_RAW")
log "  Exones extraídos (sin fusionar): $N_EXONES_RAW"

# ---------------------------------------------------------------------------
# 3. Ordenar y fusionar exones solapantes/adyacentes por gen
# ---------------------------------------------------------------------------
log "Paso 3/5: ordenando y fusionando exones"
sort -k1,1 -k2,2n "$EXONES_RAW" > "$OUTDIR/exones_mane_sorted.bed"

# Fusionar exones solapantes, conservando el nombre de gen (columna 4)
bedtools merge -i "$OUTDIR/exones_mane_sorted.bed" -c 4,6 -o distinct,distinct \
    2> "$OUTDIR/logs/03_bedtools_merge.log" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,".",$5}' \
    > "$EXONES_BED"

N_EXONES=$(wc -l < "$EXONES_BED")
log "  Exones tras fusionar: $N_EXONES"

# ---------------------------------------------------------------------------
# 4. Intrones "crudos" = gen completo - exones (bedtools subtract)
# ---------------------------------------------------------------------------
log "Paso 4/5: calculando intrones (gen completo - exones)"
sort -k1,1 -k2,2n "$GENES_BED" > "$OUTDIR/genes_sorted.bed"

bedtools subtract \
    -a "$OUTDIR/genes_sorted.bed" \
    -b "$EXONES_BED" \
    2> "$OUTDIR/logs/04_bedtools_subtract.log" \
    | sort -k1,1 -k2,2n \
    > "$INTRONES_CRUDOS"

N_INTRONES_CRUDOS=$(wc -l < "$INTRONES_CRUDOS")
log "  Intrones crudos generados: $N_INTRONES_CRUDOS"

# ---------------------------------------------------------------------------
# 5. Añadir padding hacia el exón (bedtools slop), recortando para no
#    invadir el propio intrón vecino (-i tratado como límite del gen)
# ---------------------------------------------------------------------------
log "Paso 5/5: ampliando intrones ±${PADDING}pb (bedtools slop) hacia las regiones flanqueantes"

# bedtools slop expande simétricamente ambos lados del intervalo.
# -g: archivo de tamaños de cromosoma, evita salirse de los límites del genoma.
bedtools slop \
    -i "$INTRONES_CRUDOS" \
    -g "$GENOME_FILE" \
    -b "$PADDING" \
    2> "$OUTDIR/logs/05_bedtools_slop.log" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - -c 4 -o distinct \
    > "$INTRONES_PADDED"

N_FINAL=$(wc -l < "$INTRONES_PADDED")
log "Intrones finales (con padding ±${PADDING}pb, fusionados): $N_FINAL"
log "BED final: $INTRONES_PADDED"

# ---------------------------------------------------------------------------
# Aviso importante sobre el padding y el exón
# ---------------------------------------------------------------------------
# bedtools slop expande el intrón hacia AMBOS lados, incluyendo unos pb hacia
# dentro del exón vecino. Esto es intencionado: cubre splice_donor_variant,
# splice_acceptor_variant y splice_region_variant, que incluyen posiciones
# tanto intrónicas como las últimas/primeras bases exónicas adyacentes.
# Si se prefiere un BED estrictamente intrónico (sin invadir el exón), usar
# bedtools slop con -l/-r asimétrico o bedtools subtract final contra
# exones_mane.bed tras el slop. Documentar la decisión tomada en la memoria.
