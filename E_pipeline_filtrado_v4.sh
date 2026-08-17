

#!/usr/bin/env bash
#
#
# Requisitos: bcftools, bedtools, tabix, recurso gnomAD (AF por variante)
#

set -euo pipefail

# ---------------------------- 0. CONFIGURACIÓN ------------------------------

# ---------------------------- Declaración de variables --------------------------

VCF_RAW=""          # VCF con todas las variantes generado en el pipeline wgs_pipeline_vX
SAMPLE_ID=""        # ID de la muestra
GENES_BED=""        # Archivo BED de genes de interés
INTRONES_BED=""     # Archivo BED de regiones intrónicas de genes de interes
REF_GENOME=""       # Genoma de referencia
GNOMAD_VCF=""       # vcf de gnomad con frecuencia alelica (AF) por variante (específico según genoma referencia empleado)
MAF_THRESHOLD=0.0001 # MAF: minimum allele frequency
THREADS=8            # Hilos de procesador utilizados: 8 por defecto
OUTDIR=""            # Directorio de salida


# ------------------------------------ Se crea función de uso ------------------
usage() {
    cat <<EOF
Uso: $(basename "$0") -v CMT1234.vcf.gz -g CMT_genes.bed -i CMT_intrones_padded.bed -r ref.fasta -c vep_cache -n gnomad.genomes.af_only.vcf.gz -x spliceai_scores.raw.snv.hg38.vcf.gz -y spliceai_scores.raw.indel.hg38.vcf.gz)

Obligatorios:
  -s ID de la muestra
  -v VCF crudo de entrada
  -g BED de genes completos
  -i BED de intrones MANE + padding
  -r Genoma de referencia FASTA
  -n VCF de recurso gnomAD (con AF por variante, GRCh38)
  -o Directorio de salida
Opcionales:

  -t Opcional. Hilos de procesador que se quieren utilizar en las herramientas que lo permitan. Por defecto 8

Ejemplo:
  $(basename "$0") \\
      -s paciente01 \\
      -v paciente01.raw.vcf.gz \\
      -g genes_CMT.bed \\
      -i intrones_mane_padded.bed \\
      -r GRCh38.fa \\
      -n gnomad.genomes.vcf.gz \\
      -m 0.0001 -t 8 -o resultados_paciente01
EOF
}

# -------------------------- Asignación de valores a variables mediante flags ------------------------------
while getopts ":s:v:g:i:r:c:n:x:y:m:t:o:h" opt; do
    case "$opt" in
        s) SAMPLE_ID="$OPTARG" ;;
        v) VCF_RAW="$OPTARG" ;;
        g) GENES_BED="$OPTARG" ;;
        i) INTRONES_BED="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        n) GNOMAD_VCF="$OPTARG" ;;
        m) MAF_THRESHOLD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) mostrar_ayuda; exit 0 ;;
        \?) echo "ERROR: opción no reconocida: -$OPTARG" >&2; mostrar_ayuda; exit 1 ;;
        :)  echo "ERROR: la opción -$OPTARG requiere un valor" >&2; mostrar_ayuda; exit 1 ;;
    esac
done

# ---------------------- Comprobación de que las variables obligatorias están definidas ------------------------

[[ -z "$SAMPLE_ID" || -z "VCF_RAW" || -z "GENES_BED" || -z "INTRONES_BED" || -z "REF_GENOME" || -z "GNOMAD_VCF" || -z "MAF_THERESHOLD"  ]] && usage
# comprueba si las variables especificadas estan vacias (-z, zero length es TRUE)
# Si lo están ejecuta usage y se para el pipeline

# -------------------- Crear directorio y subsirectorios específicos para la muestra ---------------------

mkdir -p "$OUTDIR"/{01_genes,02_intronicas,03_maf,logs}

# ---------------------------- Creación de la función log ------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# 1. FILTRO DE GENES (BED de genes completos)
# ---------------------------------------------------------------------------
VCF_GENES="$OUTDIR/01_genes/${SAMPLE_ID}.genes_CMT.vcf.gz"

log "Paso 1/3: filtrando variantes por lista de genes CMT ($GENES_BED)"
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

log "Paso 2/3: separando variantes intrónicas ($INTRONES_BED)"
bedtools intersect -header -a "$VCF_GENES" -b "$INTRONES_BED" -u \
    2> "$OUTDIR/logs/02_bedtools_intronicas.log" \
    | bgzip > "$VCF_INTRONICAS"
tabix -p vcf "$VCF_INTRONICAS"

N_INTRON=$(bcftools view -H "$VCF_INTRONICAS" | wc -l)
log "  Variantes intrónicas (rama seleccionada): $N_INTRON"

# ------------------------------------------------------------------------------
# 3. ANOTAR CON MAF DE gnomAD Y FILTRAR
# ------------------------------------------------------------------------------
VCF_MAF_ANOT="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.gnomad_annot.vcf.gz"
VCF_MAF="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.maf_filtered.vcf.gz"
log "Paso 3/3 Anotando variantes con frecuencia alélica de gnomAD..."

bcftools annotate \
    -a "$GNOMAD_VCF" \
    -c CHROM,POS,REF,ALT,INFO/AF \
    -Oz \
    -o "$VCF_MAF_ANOT" \
    "$VCF_INTRONICAS"

tabix -f -p vcf "$VCF_MAF_ANOT"

log "Paso 3/3 Filtrando variantes con MAF < ${MAF_THRESHOLD}..."

bcftools view \
    -i "INFO/AF<${MAF_THRESHOLD} || INFO/AF='.'" \
    -Oz \
    -o "$VCF_MAF" \
    "$VCF_MAF_ANOT"

tabix -f -p vcf "$VCF_MAF"

log "Variantes tras filtrado por MAF:"
bcftools view -H "$VCF_MAF" | wc -l



