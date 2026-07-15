#!/usr/bin/env bash
################################################################################
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
# Se crean las variables vacías que se necesitan

REF_GENOME=""          # ruta al genoma de referencia
SAMPLE_ID=""
FASTQ_R1=""
FASTQ_R2=""
THREADS=8              # Son los hilos, es decir, unidades de ejecución paralela. Tiene el valor por defecto 4 (evaluar cuantos poner según el ordenador que vaya a utilizar)
OUTDIR="resultados_$SAMPLE_ID"    # tiene el valor por defecto "resultados"
MAX_TRIM_ROUNDS=1      # nº máximo de intentos de recorte antes de abortar

# Módulos de FastQC que se consideran críticos para decidir si hace falta recortar.
# Todos los  módulos posibles
#  - Basic Statistics
#  - Per base sequence quality
#  - Per sequence quality scores
#  - Per base sequence content
#  - Per sequence GC content
#  - Per base N content
#  - Sequence Length Distribution
#  - Sequence Duplication Levels
#  - Overrepresented sequences
#  - Adapter Content
MODULOS_CRITICOS=(
    "Per base sequence quality"
    "Per sequence quality scores"
    "Adapter Content"
)

# --- Parámetros de recorte/filtrado de fastp  ---
# Ajusta estos valores según el criterio de tu pipeline.
FASTP_CUT_WINDOW_SIZE=4            #     : tamaño de ventana para recorte por calidad deslizante (default fastp: 4)
FASTP_CUT_MEAN_QUALITY=20          #     : calidad media mínima exigida en cada ventana (default fastp: 20, Q20)





# se crea la funcion de uso
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
        h) usage ;;                  # si se detecta -h (help) como flag en el comando inicial, imprime usage
        *) usage ;;                  # si se detecta una flag no definica *, imprime usage
    esac
done

[[ -z "$REF_GENOME" || -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$SAMPLE_ID" ]] && usage 
# comprueba que las variables que necesitan valor no estan vacias (-z)


mkdir -p "$OUTDIR"/{qc_raw,trimmed,qc_trimmed,align,qc_align,dedup,variants,logs} 
# se crea directorio con la variable OUTDIR (valor por defecto resultados)  y los subdirectorios

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$OUTDIR/logs/pipeline.log"; }
# se crea la función log para trackear el proceso que se está ejecutando, el día y la hora



# ------------------------ 1. CARGA / VARIABLES -------------------------------

log "== Iniciando pipeline para muestra: $SAMPLE_ID =="
R1="$FASTQ_R1"
R2="$FASTQ_R2"

# ------------------------ 2. ANÁLISIS DE CALIDAD (FastQC) -------------------

# Función para correr el analisis de calidad (fastqc)
# se definen las variables locales r1 r2 y outdir
run_fastqc() {
    local r1=$1 r2=$2 outdir=$3
    log "Ejecutando FastQC sobre $r1 y $r2"
    fastqc -t "$THREADS" -o "$outdir" "$r1" "$r2" &> "$OUTDIR/logs/fastqc_${outdir##*/}.log" 
}

# &> manda los mensajes de progreso a un archivo en la carpeta logs

# Se crea funcion qc_pasa_criterios
# Comprueba si FastQC ha marcado FAIL en alguno de los módulos definidos como críticos
# en MODULOS_CRITICOS (ver configuración arriba). Módulos no listados ahí no bloquean
# el pipeline aunque FastQC los marque como FAIL o WARN.
qc_pasa_criterios() {
    local outdir=$1
    local fail=0
    for zip in "$outdir"/*_fastqc.zip; do
        unzip -p "$zip" "*/summary.txt" > /tmp/summary.txt
        for modulo in "${MODULOS_CRITICOS[@]}"; do
            if grep -qP "^FAIL\t${modulo}\t" /tmp/summary.txt; then
                log "FastQC FAIL en módulo crítico '$modulo' ($(basename "$zip"))"
                fail=1
            fi
        done
    done
    return $fail   # 0 = pasa (sin FAIL en módulos críticos), 1 = no pasa
}


# en la funcion qc_pasa_criterios
# se definen las variables locales outdir y fail. fail por defecto es 0 (no falla nada)
# para todos los archivos que acaben en _fastqc.zip de la ruta definida en outdir, lo descomprime y extrae el summary.txt 
# con grep comrpueba si exite el patron FAIL seguido del nombre de cada uno de los  modulos criticos en el summary.text (comprueba cada modulo uno a uno)
# Si existe ejecuta el log  con el mensaje de error y cambia la variable fail a 1


# ------------------------ 3. BUCLE QC -> RECORTE (¿Cumple criterios?) --------

round=0
qc_dir="$OUTDIR/qc_raw"
# definimos la variables qc_dir que es el nombre de la ruta donde tiene que guardar los archivos la funcion run_fastqc
run_fastqc "$R1" "$R2" "$qc_dir"
# le damos valores a las variables locales de la funcion run_fastqc


# FASTP: es una herramienta utilizada concretamente para short reads (Illumina, Sanger)


while ! qc_pasa_criterios "$qc_dir"; do
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
          --cut_front \
          --cut_tail \
          --cut_window_size "$FASTP_CUT_WINDOW_SIZE" \
          --cut_mean_quality "$FASTP_CUT_MEAN_QUALITY" \
          --detect_adapter_for_pe \
          --json "$OUTDIR/trimmed/${SAMPLE_ID}_fastp.json" \
          --html "$OUTDIR/trimmed/${SAMPLE_ID}_fastp.html" \
          &> "$OUTDIR/logs/fastp_round${round}.log"


    R1="$OUTDIR/trimmed/${SAMPLE_ID}_R1.trim.fastq.gz"
    R2="$OUTDIR/trimmed/${SAMPLE_ID}_R2.trim.fastq.gz"

    # A partir de la primera ronda de recorte, el QC se evalúa en qc_trimmed,
    # dejando qc_raw intacto con el QC de los datos originales sin tocar.
    qc_dir="$OUTDIR/qc_trimmed"
    rm -rf "$qc_dir"/*
    run_fastqc "$R1" "$R2" "$qc_dir"
done
log "Sí cumple criterios de calidad -> continuando con alineamiento"

# Bucle while: mientras los datos no cumplan los criterios de calidad, 
# sigue recortando con fastp y volviendo a comprobar, hasta un máximo de intentos.




# ------------------------ 4. ALINEAMIENTO (BWA) ------------------------------
# Se usa BWA-MEM por ser el estándar actual para WGS de lectura corta.

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

# ------------------------ 6. MARCAR DUPLICADOS PCR (Picard) ---------------

BAM_DEDUP="$OUTDIR/dedup/${SAMPLE_ID}.dedup.bam"

log "Marcando duplicados PCR con Picard MarkDuplicates"
picard MarkDuplicates \
    I="$BAM_RAW" \
    O="$BAM_DEDUP" \
    M="$OUTDIR/dedup/${SAMPLE_ID}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    &> "$OUTDIR/logs/picard_markdup.log"

samtools index "$BAM_DEDUP"

# ------------------------ 7. LLAMADO DE VARIANTES (GATK HaplotypeCaller) ----


VCF_OUT="$OUTDIR/variants/${SAMPLE_ID}.g.vcf.gz"
log "Llamando variantes germinales con GATK HaplotypeCaller"
gatk HaplotypeCaller \
    -R "$REF_GENOME" \
    -I "$BAM_DEDUP" \
    -O "$VCF_OUT" \
    &> "$OUTDIR/logs/gatk_haplotypecaller.log"

log "== Pipeline completado. VCF final: $VCF_OUT =="
