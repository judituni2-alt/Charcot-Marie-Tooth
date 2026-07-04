#!/usr/bin/env bash
###############################################################################
# wgs_pipeline.sh
#
# Pipeline de análisis de Secuenciación de Genoma Completo (WGS) - línea germinal
# Secuenciación -> FASTQ -> QC -> Alineamiento -> QC BAM -> Duplicados -> VCF
#
# Uso:
#   ./wgs_pipeline.sh -1 CMT1234_R1.fastq.gz -2 CMT1234_R2.fastq.gz -s CMT1234
#
#   Ejemplo uso: se cargan los archivos en el directorio actual y se añade
#   su ID de modo que. Este ID servira para identificar los archivos que 
#   correspondan a un mismo paciente.
# Requiere (instalar con conda/mamba idealmente):
#   fastqc, fastp, bwa, samtools, bamtools, picard, gatk4
###############################################################################

set -euo pipefail 
# activa 3 procedimientos de seguridad en caso de que el script falle
#  -e: si cualquier comando del script devuelve un error, el script se detiene inmediatamente
#  -u: si se usa una variable que no ha sido definida el script falla
#  -o pipefail: cuando se encadadena con tuberias, si cualquiera de los comandos de la cadena falla, toda la tubería
#   se considera fallida



# ---------------------------- 0. CONFIGURACIÓN ------------------------------

REF_GENOME=""          # ruta al genoma de referencia (.fasta, ya indexado con bwa index)
SAMPLE_ID=""
FASTQ_R1=""
FASTQ_R2=""
THREADS=4
OUTDIR="resultados"
MAX_TRIM_ROUNDS=2       # nº máximo de intentos de recorte antes de abortar

usage() {
    echo "Uso: $0 -r ref.fasta -1 R1.fastq.gz -2 R2.fastq.gz -s sample_id [-t threads] [-o outdir]"
    exit 1
}

while getopts "r:1:2:s:t:o:h" opt; do
    case $opt in
        r) REF_GENOME="$OPTARG" ;;
        1) FASTQ_R1="$OPTARG" ;;
        2) FASTQ_R2="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$REF_GENOME" || -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$SAMPLE_ID" ]] && usage

mkdir -p "$OUTDIR"/{qc_raw,trimmed,qc_trimmed,align,qc_align,dedup,variants,logs}

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$OUTDIR/logs/pipeline.log"; }

# ------------------------ 1. CARGA / VARIABLES -------------------------------

log "== Iniciando pipeline para muestra: $SAMPLE_ID =="
R1="$FASTQ_R1"
R2="$FASTQ_R2"

# ------------------------ 2. ANÁLISIS DE CALIDAD (FastQC) -------------------

run_fastqc() {
    local r1=$1 r2=$2 outdir=$3
    log "Ejecutando FastQC sobre $r1 y $r2"
    fastqc -t "$THREADS" -o "$outdir" "$r1" "$r2" &> "$OUTDIR/logs/fastqc_${outdir##*/}.log"
}

# Comprueba si FastQC ha marcado FAIL en algún módulo clave (adaptadores, calidad por base)
qc_pasa_criterios() {
    local outdir=$1
    local fail=0
    for zip in "$outdir"/*_fastqc.zip; do
        unzip -p "$zip" "*/summary.txt" > /tmp/summary.txt
        if grep -qE "^FAIL" /tmp/summary.txt; then
            fail=1
        fi
    done
    return $fail   # 0 = pasa (sin FAIL), 1 = no pasa
}

# ------------------------ 3. BUCLE QC -> RECORTE (¿Cumple criterios?) --------

round=0
run_fastqc "$R1" "$R2" "$OUTDIR/qc_raw"

while ! qc_pasa_criterios "$OUTDIR/qc_raw"; do
    round=$((round+1))
    if [[ $round -gt $MAX_TRIM_ROUNDS ]]; then
        log "ERROR: la calidad sigue sin cumplir criterios tras $MAX_TRIM_ROUNDS recortes. Abortando."
        exit 1
    fi
    log "No cumple criterios de calidad -> recortando con fastp (ronda $round)"
    fastp -i "$R1" -I "$R2" \
          -o "$OUTDIR/trimmed/${SAMPLE_ID}_R1.trim.fastq.gz" \
          -O "$OUTDIR/trimmed/${SAMPLE_ID}_R2.trim.fastq.gz" \
          --thread "$THREADS" \
          --json "$OUTDIR/trimmed/${SAMPLE_ID}_fastp.json" \
          --html "$OUTDIR/trimmed/${SAMPLE_ID}_fastp.html" \
          &> "$OUTDIR/logs/fastp_round${round}.log"

    R1="$OUTDIR/trimmed/${SAMPLE_ID}_R1.trim.fastq.gz"
    R2="$OUTDIR/trimmed/${SAMPLE_ID}_R2.trim.fastq.gz"

    rm -rf "$OUTDIR/qc_raw"/*
    run_fastqc "$R1" "$R2" "$OUTDIR/qc_raw"
done
log "Sí cumple criterios de calidad -> continuando con alineamiento"

# ------------------------ 4. ALINEAMIENTO (BWA) ------------------------------
# Nota: el diagrama contempla BWA / SOAP / Bowtie / MOSAIK. Se usa BWA-MEM
# por ser el estándar actual para WGS de lectura corta.

BAM_RAW="$OUTDIR/align/${SAMPLE_ID}.sorted.bam"

run_alineamiento() {
    log "Alineando con BWA-MEM contra $REF_GENOME"
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA" \
        "$REF_GENOME" "$R1" "$R2" 2> "$OUTDIR/logs/bwa.log" \
        | samtools sort -@ "$THREADS" -o "$BAM_RAW" -
    samtools index "$BAM_RAW"
}
run_alineamiento

# ------------------------ 5. ANÁLISIS DE CALIDAD DEL BAM (BAMtools) ---------

run_qc_bam() {
    log "Ejecutando control de calidad del BAM (bamtools stats + samtools flagstat)"
    bamtools stats -in "$BAM_RAW" > "$OUTDIR/qc_align/${SAMPLE_ID}_bamtools_stats.txt"
    samtools flagstat "$BAM_RAW" > "$OUTDIR/qc_align/${SAMPLE_ID}_flagstat.txt"
}
run_qc_bam

# Criterio simple de calidad del BAM: % de lecturas mapeadas >= 90%
bam_pasa_criterios() {
    local pct
    pct=$(grep "mapped (" "$OUTDIR/qc_align/${SAMPLE_ID}_flagstat.txt" | head -1 \
          | grep -oP '\(\K[0-9.]+(?=%)')
    log "Porcentaje de lecturas mapeadas: ${pct}%"
    awk -v p="$pct" 'BEGIN{exit !(p>=90)}'
}

if ! bam_pasa_criterios; then
    log "No cumple criterios de calidad del BAM -> se recomienda revisar el alineamiento (parámetros, referencia)"
    log "El pipeline continúa pero revisa $OUTDIR/qc_align antes de confiar en el VCF final."
else
    log "Sí cumple criterios de calidad del BAM -> continuando"
fi

# ------------------------ 6. ELIMINAR DUPLICADOS PCR (Picard) ---------------

BAM_DEDUP="$OUTDIR/dedup/${SAMPLE_ID}.dedup.bam"

log "Eliminando duplicados PCR con Picard MarkDuplicates"
picard MarkDuplicates \
    I="$BAM_RAW" \
    O="$BAM_DEDUP" \
    M="$OUTDIR/dedup/${SAMPLE_ID}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    &> "$OUTDIR/logs/picard_markdup.log"

samtools index "$BAM_DEDUP"

# ------------------------ 7. LLAMADO DE VARIANTES (GATK HaplotypeCaller) ----
# Nota: el diagrama contempla GATK HaplotypeCaller / FreeBayes / Platypus /
# Samtools-BCFtools. Se usa GATK por ser el estándar para línea germinal.

VCF_OUT="$OUTDIR/variants/${SAMPLE_ID}.g.vcf.gz"

log "Llamando variantes germinales con GATK HaplotypeCaller"
gatk HaplotypeCaller \
    -R "$REF_GENOME" \
    -I "$BAM_DEDUP" \
    -O "$VCF_OUT" \
    &> "$OUTDIR/logs/gatk_haplotypecaller.log"

log "== Pipeline completado. VCF final: $VCF_OUT =="
