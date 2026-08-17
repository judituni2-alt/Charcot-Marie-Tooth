#!/usr/bin/env bash




# ---------------------------------------------------------------------------
# ANOTACIÓN VEP + SpliceAI (solo anotación: la separación intrón/exón ya
# se hizo por BED en el paso 2, aquí VEP no filtra nada)
# Sin filtro de sinónimas: no aplica a variantes intrónicas.
# ---------------------------------------------------------------------------

set -euo pipefail

# ---------------------------- 0. CONFIGURACIÓN ------------------------------


VCF_MAF=""          # VCF con variantes filtradas por MAF y Genes de interés
VEP_CACHE_DIR=""    # Directorio con caché vep
REF_GENOME=""       # Genoma de referencia
THREADS=8            # Hilos de procesador utilizados: 8 por defecto
SPLICEAI_SNV=""     # Recurso SpliceAI para SNV
SPLICEAI_INDEL=""   # Recurso SpliceAI para indels
OUTDIR=""            # Directorio de salida


# -------------------------- Asignación de valores a variables mediante flags ------------------------------
while getopts "s:v:r:c:x:y:t:o:h" opt; do
    case "$opt" in
        s) SAMPLE_ID="$OPTARG" ;;
        v) VCF_MAF="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        c) VEP_CACHE_DIR="$OPTARG" ;;
        x) SPLICEAI_SNV="$OPTARG" ;;
        y) SPLICEAI_INDEL="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) mostrar_ayuda; exit 0 ;;
        \?) echo "ERROR: opción no reconocida: -$OPTARG" >&2; mostrar_ayuda; exit 1 ;;
        :)  echo "ERROR: la opción -$OPTARG requiere un valor" >&2; mostrar_ayuda; exit 1 ;;
    esac
done


# -------------------- Crear directorio y subsirectorios específicos para la muestra ---------------------

mkdir -p "$OUTDIR"/{logs,vep_spliceai}

# ---------------------------- Creación de la función log ------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

VCF_FINAL="$OUTDIR/vep_spliceai/${SAMPLE_ID}.intronicas.spliceai.vcf.gz"

log "Anotando con VEP (MANE Select) + SpliceAI"
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

