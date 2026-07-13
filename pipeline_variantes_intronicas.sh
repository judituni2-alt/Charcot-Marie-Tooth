#!/usr/bin/env bash
#
# pipeline_variantes_intronicas.sh
#
# Flujo (rama intrónica únicamente):
#
#   VCF crudo (post-QC)
#         |
#   1. Filtro de genes (BED genes CMT, gen completo)
#         |
#   2. Separar intrónico (BED intrones_mane_padded.bed)
#         |
#   3. Filtro MAF (gnomAD)
#         |
#   4. Anotación VEP + SpliceAI (solo anotación, ya no filtra consecuencia)
#      (sin filtro de sinónimas: no aplica a intrones)
#         |
#   Variantes intrónicas candidatas (VCF final)
#
# Requisitos: bcftools, bedtools, tabix, vep (con cache offline + plugin SpliceAI),
#             recurso gnomAD (AF por variante), recursos SpliceAI (snv/indel).
#
# El BED de intrones (intrones_mane_padded.bed) debe generarse PREVIAMENTE
# con generar_bed_intrones.sh, a partir de genes_CMT.bed + GTF de MANE.

set -euo pipefail

# ---------------------------------------------------------------------------
# PARSEO DE FLAGS con getopts 
# ---------------------------------------------------------------------------
mostrar_ayuda() {
    cat <<EOF
Uso: $(basename "$0") -v FILE -g FILE -i FILE -r FILE -c DIR -n FILE -x FILE -y FILE [opciones]

Obligatorios:
  -v FILE   VCF crudo de entrada (post-QC general)
  -g FILE   BED de genes completos (ej. genes_CMT.bed)
  -i FILE   BED de intrones MANE + padding (generado con generar_bed_intrones.sh)
  -r FILE   Genoma de referencia FASTA (GRCh38, indexado)
  -c DIR    Carpeta de cache de VEP (--dir_cache)
  -n FILE   VCF de recurso gnomAD (con AF por variante, GRCh38)
  -x FILE   Recurso SpliceAI para SNV (spliceai_scores.raw.snv.hg38.vcf.gz)
  -y FILE   Recurso SpliceAI para indels (spliceai_scores.raw.indel.hg38.vcf.gz)

Opcionales:
  -s ID     ID de la muestra (default: muestra01)
  -m FLOAT  Umbral de MAF gnomAD (default: 0.0001)
  -t N      Hilos/threads para VEP (default: 2)
  -o DIR    Carpeta de salida (default: ./resultados_<SAMPLE_ID>)
  -h        Muestra esta ayuda

También se pueden definir como variables de entorno (VCF_RAW, GENES_BED,
INTRONES_BED, REF_GENOME, VEP_CACHE_DIR, GNOMAD_VCF, SPLICEAI_SNV,
SPLICEAI_INDEL, SAMPLE_ID, MAF_THRESHOLD, THREADS, OUTDIR); los flags, si
se indican, tienen prioridad.

Ejemplo:
  $(basename "$0") \\
      -s paciente01 \\
      -v paciente01.raw.vcf.gz \\
      -g genes_CMT.bed \\
      -i intrones_mane_padded.bed \\
      -r GRCh38.fa \\
      -c /ruta/.vep \\
      -n gnomad.genomes.af_only.vcf.gz \\
      -x spliceai_scores.raw.snv.hg38.vcf.gz \\
      -y spliceai_scores.raw.indel.hg38.vcf.gz \\
      -m 0.0001 -t 4 -o resultados_paciente01
EOF
}

while getopts ":s:v:g:i:r:c:n:x:y:m:t:o:h" opt; do
    case "$opt" in
        s) SAMPLE_ID="$OPTARG" ;;
        v) VCF_RAW="$OPTARG" ;;
        g) GENES_BED="$OPTARG" ;;
        i) INTRONES_BED="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        c) VEP_CACHE_DIR="$OPTARG" ;;
        n) GNOMAD_VCF="$OPTARG" ;;
        x) SPLICEAI_SNV="$OPTARG" ;;
        y) SPLICEAI_INDEL="$OPTARG" ;;
        m) MAF_THRESHOLD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) mostrar_ayuda; exit 0 ;;
        \?) echo "ERROR: opción no reconocida: -$OPTARG" >&2; mostrar_ayuda; exit 1 ;;
        :)  echo "ERROR: la opción -$OPTARG requiere un valor" >&2; mostrar_ayuda; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# ---------------------------------------------------------------------------
# CONFIGURACIÓN (valores por flag, si no, por variable de entorno)
# Nota: SAMPLE_ID y OUTDIR se resuelven ANTES que VCF_RAW/GENES_BED, porque
# estos últimos usan $OUTDIR en su valor por defecto (orden corregido
# respecto a la versión anterior, donde OUTDIR se definía después de usarse).
# ---------------------------------------------------------------------------
SAMPLE_ID="${SAMPLE_ID:-muestra01}"
OUTDIR="${OUTDIR:-$(pwd)/resultados_${SAMPLE_ID}}"

VCF_RAW="${VCF_RAW:?Debe definir VCF_RAW (-v o variable de entorno)}"
GENES_BED="${GENES_BED:?Debe definir GENES_BED (-g o variable de entorno)}"
INTRONES_BED="${INTRONES_BED:?Debe definir INTRONES_BED (-i o variable de entorno)}"
REF_GENOME="${REF_GENOME:?Debe definir REF_GENOME (-r o variable de entorno)}"
VEP_CACHE_DIR="${VEP_CACHE_DIR:?Debe definir VEP_CACHE_DIR (-c o variable de entorno)}"
GNOMAD_VCF="${GNOMAD_VCF:?Debe definir GNOMAD_VCF (-n o variable de entorno)}"
SPLICEAI_SNV="${SPLICEAI_SNV:?Debe definir SPLICEAI_SNV (-x o variable de entorno)}"
SPLICEAI_INDEL="${SPLICEAI_INDEL:?Debe definir SPLICEAI_INDEL (-y o variable de entorno)}"
MAF_THRESHOLD="${MAF_THRESHOLD:-0.0001}"
THREADS="${THREADS:-2}"

mkdir -p "$OUTDIR"/{01_genes,02_intronicas,03_maf,04_spliceai,logs}

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# 1. FILTRO DE GENES (BED de genes completos)
# ---------------------------------------------------------------------------
VCF_GENES="$OUTDIR/01_genes/${SAMPLE_ID}.genes_CMT.vcf.gz"

log "Paso 1/4: filtrando variantes por lista de genes CMT ($GENES_BED)"
bedtools intersect -header -a "$VCF_RAW" -b "$GENES_BED" -u \
    2> "$OUTDIR/logs/01_bedtools_genes.log" \
    | bgzip > "$VCF_GENES"
tabix -p vcf "$VCF_GENES"

N_GENES=$(bcftools view -H "$VCF_GENES" | wc -l)
log "  Variantes tras filtro de genes: $N_GENES"

# ---------------------------------------------------------------------------
# 2. SEPARAR RAMA INTRÓNICA (BED de intrones MANE + padding)
# ---------------------------------------------------------------------------
VCF_INTRONICAS="$OUTDIR/02_intronicas/${SAMPLE_ID}.intronicas.vcf.gz"

log "Paso 2/4: separando variantes intrónicas ($INTRONES_BED)"
bedtools intersect -header -a "$VCF_GENES" -b "$INTRONES_BED" -u \
    2> "$OUTDIR/logs/02_bedtools_intronicas.log" \
    | bgzip > "$VCF_INTRONICAS"
tabix -p vcf "$VCF_INTRONICAS"

N_INTRON=$(bcftools view -H "$VCF_INTRONICAS" | wc -l)
log "  Variantes intrónicas (rama seleccionada): $N_INTRON"

# ---------------------------------------------------------------------------
# 3. FILTRO DE FRECUENCIA POBLACIONAL (MAF, gnomAD)
# ---------------------------------------------------------------------------
VCF_MAF_ANOT="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.gnomad_annot.vcf.gz"
VCF_MAF="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.maf_filtered.vcf.gz"

log "Paso 3/4: anotando y filtrando por MAF gnomAD (< $MAF_THRESHOLD)"
bcftools annotate \
    -a "$GNOMAD_VCF" \
    -c INFO/AF_gnomad \
    -O z -o "$VCF_MAF_ANOT" \
    "$VCF_INTRONICAS" \
    2> "$OUTDIR/logs/03_bcftools_annotate_gnomad.log"
tabix -p vcf "$VCF_MAF_ANOT"

# Se conservan variantes sin AF_gnomad anotado (candidatas ausentes de gnomAD)
bcftools filter \
    -e "INFO/AF_gnomad > $MAF_THRESHOLD" \
    -O z -o "$VCF_MAF" \
    "$VCF_MAF_ANOT" \
    2> "$OUTDIR/logs/03_bcftools_filter_maf.log"
tabix -p vcf "$VCF_MAF"

N_MAF=$(bcftools view -H "$VCF_MAF" | wc -l)
log "  Variantes intrónicas tras filtro MAF: $N_MAF"

# ---------------------------------------------------------------------------
# 4. ANOTACIÓN VEP + SpliceAI (solo anotación: la separación intrón/exón ya
#    se hizo por BED en el paso 2, aquí VEP no filtra nada)
#    Sin filtro de sinónimas: no aplica a variantes intrónicas.
# ---------------------------------------------------------------------------
VCF_FINAL="$OUTDIR/04_spliceai/${SAMPLE_ID}.intronicas.spliceai.vcf.gz"

log "Paso 4/4: anotando con VEP (MANE Select) + SpliceAI"
vep --input_file "$VCF_MAF" \
    --output_file STDOUT \
    --vcf --compress_output bgzip \
    --cache --offline --dir_cache "$VEP_CACHE_DIR" \
    --fasta "$REF_GENOME" \
    --fork "$THREADS" \
    --mane_select \
    --pick --pick_order mane_select,canonical,rank \
    --plugin SpliceAI,snv="$SPLICEAI_SNV",indel="$SPLICEAI_INDEL" \
    2> "$OUTDIR/logs/04_vep_spliceai.log" \
    > "$VCF_FINAL"
tabix -p vcf "$VCF_FINAL"

N_FINAL=$(bcftools view -H "$VCF_FINAL" | wc -l)
log "Pipeline completado. Variantes intrónicas candidatas finales: $N_FINAL"
log "Resultado: $VCF_FINAL"
