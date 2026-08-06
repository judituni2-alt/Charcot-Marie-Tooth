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
#   4. Anotación VEP + SpliceAI (solo anotación)
#         |
#   Variantes intrónicas candidatas (VCF final)
#
# Requisitos: bcftools, bedtools, tabix, vep (con cache offline + plugin SpliceAI),
#             recurso gnomAD (AF por variante), recursos SpliceAI (snv/indel).
#
# El BED de intrones (intrones_mane_padded.bed) debe generarse PREVIAMENTE
# con generar_bed_intrones.sh, a partir de genes_CMT.bed + GTF de MANE.

set -euo pipefail

# ---------------------------- 0. CONFIGURACIÓN ------------------------------

# ---------------------------- Declaración de variables --------------------------

VCF_RAW=""          # VCF con todas las variantes generado en el pipeline wgs_pipeline_vX
SAMPLE_ID=""        # ID de la muestra
GENES_BED=""        # Archivo BED de genes de interés
INTRONES_BED=""     # Archivo BED de regiones intrónicas de genes de interes
REF_GENOME=""       # Genoma de referencia
VEP_CACHE_DIR=""    # Directorio con caché vep
GNOMAD_VCF=""       # vcf de gnomad con frecuencia alelica (AF) por variante (específico según genoma referencia empleado)
SPLICEAI_SNV=""     # Recurso SpliceAI para SNV
SPLICEAI_INDEL=""   # Recurso SpliceAI para indels
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
  -c Carpeta de cache de VEP
  -n VCF de recurso gnomAD (con AF por variante, GRCh38)
  -x Recurso SpliceAI para SNV (spliceai_scores.raw.snv.hg38.vcf.gz)
  -y Recurso SpliceAI para indels (spliceai_scores.raw.indel.hg38.vcf.gz)
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
      -c /ruta/.vep \\
      -n gnomad.genomes.vcf.gz \\
      -x spliceai_scores.raw.snv.hg38.vcf.gz \\
      -y spliceai_scores.raw.indel.hg38.vcf.gz \\
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

# ---------------------- Comprobación de que las variables obligatorias están definidas ------------------------

[[ -z "$SAMPLE_ID" || -z "VCF_RAW" || -z "GENES_BED" || -z "INTRONES_BED" || -z "REF_GENOME" || -z "VEP_CACHE_DIR" || -z "GNOMAD_VCF" || -z "SPLICEAI_SNV"  || -z "SPLICEAI_INDEL"  || -z "MAF_THERESHOLD"  ]] && usage
# comprueba si las variables especificadas estan vacias (-z, zero length es TRUE)
# Si lo están ejecuta usage y se para el pipeline

# -------------------- Crear directorio y subsirectorios específicos para la muestra ---------------------

mkdir -p "$OUTDIR"/{01_genes,02_intronicas,03_maf,04_spliceai,logs}

# ---------------------------- Creación de la función log ------------------------------
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

# bedtools intersect: compara dos archivos de coordenadas genómicas, 
# el VCF original (-a $VCF_RAW) y el BED de genes (-b $GENES_BED).
# -header conserva las líneas de cabecera del VCF original en la salida
# -u asegura que las variantes que caen en genes solapante se cuenten 1 sola vex
# se comprime el resulfato con bgzip y se crea un archivo vcf.tbi con tabix

# bcftools view -H "$VCF_GENES" | wc -l quita la cabecera y cuenta el numero de lineas
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

# mismo procedimiento que antes para filtrar las regiones intrónicas

# ---------------------------------------------------------------------------
# 3. FILTRO DE FRECUENCIA POBLACIONAL (MAF, gnomAD)
# ---------------------------------------------------------------------------
VCF_MAF_ANOT="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.gnomad_annot.vcf.gz"
VCF_MAF="$OUTDIR/03_maf/${SAMPLE_ID}.intronicas.maf_filtered.vcf.gz"
# Primero se anota cada variante con su MAF
log "Paso 3/4: anotando y filtrando por MAF gnomAD (< $MAF_THRESHOLD)"
bcftools annotate \
    -a "$GNOMAD_VCF" \
    -c INFO/AF_gnomad \
    -O z -o "$VCF_MAF_ANOT" \
    "$VCF_INTRONICAS" \
    2> "$OUTDIR/logs/03_bcftools_annotate_gnomad.log"
tabix -p vcf "$VCF_MAF_ANOT"

# Se filtra segun el MAF threshold definido anteriormente
# Se conservan variantes sin AF_gnomad anotado (candidatas ausentes de gnomAD)
bcftools filter \
    -e "INFO/AF_gnomad > $MAF_THRESHOLD" \
    -O z -o "$VCF_MAF" \
    "$VCF_MAF_ANOT" \
    2> "$OUTDIR/logs/03_bcftools_filter_maf.log"
tabix -p vcf "$VCF_MAF"

# Conteo de numero de variantes tras el filtro
N_MAF=$(bcftools view -H "$VCF_MAF" | wc -l)
log "  Variantes intrónicas tras filtro MAF: $N_MAF"

# ---------------------------------------------------------------------------
# 4. ANOTACIÓN VEP
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
