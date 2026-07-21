#!/usr/bin/env bash
################################################################################
# wgs_pipeline.sh
#
# Pipeline de análisis de Secuenciación de Genoma Completo (WGS) - línea germinal
# Secuenciación -> Fusión FASTQ -> QC -> Alineamiento -> QC BAM -> Duplicados -> VCF
#
# Los archivos FASTQ de entrada se leen desde un disco duro/unidad externa
# y todos los archivos y carpetas que genera el pipeline se
# guardan también en ese mismo disco.
#
# Uso:
#   ./B_wgs_pipeline.sh
#       -1 CMT1234_R1.fastq.gz -2 CMT1234_R2.fastq.gz -s CMT1234 -r ref.fasta
#
#
#
#   Si la muestra viene repartida en varios archivos por carril de
#   secuenciación (lane splitting, típico de Illumina: L001, L002, L003...),
#   se pueden indicar varios archivos separados por coma en -1 y -2, y el
#   pipeline los fusiona automáticamente en un único R1 y un único R2 antes
#   de continuar:
#
#   ./B_wgs_pipeline.sh
#       -1 CMT1234_L001_R1.fastq.gz,CMT1234_L002_R1.fastq.gz \
#       -2 CMT1234_L001_R2.fastq.gz,CMT1234_L002_R2.fastq.gz \
#       -s CMT1234 -r ref.fasta
#       -o resultados
#
#   Para que el pipeline siga ejecutándose aunque se cierre la terminal
#   (por ejemplo, en una conexión SSH que se pueda cortar), añadir -b:
#
#   ./B_wgs_pipeline.sh -1 ... -2 ... -s CMT1234 -r ref.fasta -b
#
# Requiere instalar: fastqc, fastp, bwa, samtools, bamtools, picard, bcftools, tabix
# Para ello crear un ambiente de conda/mamba con los recursos necesarios
#
# NOTA SOBRE VELOCIDAD DE FASTQC:
# FastQC paraleliza por ARCHIVO (cada archivo lo procesa un único hilo), no
# por contenido interno del archivo, así que con solo R1+R2 nunca aprovecha
# más de 2 hilos por mucho -t que se indique, y leer/descomprimir el FASTQ
# completo es lo que más tiempo tarda. Por eso, antes de cada llamada a
# FastQC, se submuestrean SUBSAMPLE_READS reads por archivo (R1 y R2 en
# paralelo): unos pocos millones de lecturas ya dan métricas de calidad
# estadísticamente representativas y se evita leer el FASTQ entero.
# Para desactivar el submuestreo y analizar el archivo completo, poner
# SUBSAMPLE_READS=0.
###############################################################################

set -euo pipefail
# activa 3 procedimientos de seguridad en caso de que el script falle
#  -e: si cualquier comando del script devuelve un error (exit status distinto de 0) el script se detiene inmediatamente
#  -u: si se usa una variable que no ha sido definida previamente el script falla
#  -o pipefail: cuando se encadadena con tuberias, si cualquiera de los comandos de la cadena falla, toda la tubería
#   se considera fallida
# https://linuxize.com/post/bash-best-practices/  --> Investigar mas buenas practica


# ---------------------------- 0. CONFIGURACIÓN ------------------------------
# Se crean las variables vacías que se necesitan

REF_GENOME=""          # ruta al genoma de referencia
SAMPLE_ID=""           # ID de la muestra
FASTQ_R1=""            # uno o varios archivos separados por coma (lanes)
FASTQ_R2=""            # uno o varios archivos separados por coma (lanes)
THREADS=8              # Son los hilos, es decir, unidades de ejecución paralela. Valor por defecto: 8
OUTDIR=""              # se resuelve más abajo, DESPUÉS de parsear -s
MAX_TRIM_ROUNDS=1      # nº máximo de intentos de recorte antes de abortar
BACKGROUND=false       # si se activa con -b, el pipeline se relanza en segundo plano
SUBSAMPLE_READS=2000000  # nº de reads a submuestrear por archivo antes de FastQC (0 = desactivado, usa el archivo completo)

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
    "Overrepresented sequences"
)

# --- Parámetros de recorte/filtrado de fastp  ---
# Ajusta estos valores según el criterio de tu pipeline.
FASTP_CUT_WINDOW_SIZE=4            #     : tamaño de ventana para recorte por calidad deslizante (default fastp: 4)
FASTP_CUT_MEAN_QUALITY=20          #     : calidad media mínima exigida en cada ventana (default fastp: 20, Q20)





# se crea la funcion de uso
usage() {
    cat <<EOF
Uso: $0  -r ref.fasta -1 R1.fastq.gz -2 R2.fastq.gz -s sample_id [-t threads] [-o outdir] [-b] [-q subsample_reads]

  -1 / -2   Admiten varios archivos separados por coma (sin espacios) si la
            muestra está repartida en varios carriles/lanes, ej:
            -1 L001_R1.fastq.gz,L002_R1.fastq.gz
  -b        Ejecuta el pipeline en segundo plano (nohup); puedes cerrar la
            terminal sin que el proceso se interrumpa.
  -q        Nº de reads a submuestrear por archivo antes de FastQC
            (por defecto 2000000). Usar -q 0 para analizar el FASTQ completo.
EOF
    exit 1
}

while getopts "r:1:2:s:t:o:q:bh" opt; do
    case $opt in
        r) REF_GENOME="$OPTARG" ;;
        1) FASTQ_R1="$OPTARG" ;;
        2) FASTQ_R2="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        q) SUBSAMPLE_READS="$OPTARG" ;;
        b) BACKGROUND=true ;;
        h) usage ;;                  # si se detecta -h (help) como flag en el comando inicial, imprime usage
        *) usage ;;                  # si se detecta una flag no definica *, imprime usage
    esac
done

[[ -z "$REF_GENOME" || -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$SAMPLE_ID" ]] && usage
# comprueba si las variables REF_GENOME o FASTQ_R1 o FASTQ_R2 o SAMPLE_ID estan vacias (-z es TRUE) y si lo están ejecuta usage

mkdir -p "$OUTDIR"/{merged,qc_raw,trimmed,qc_trimmed,align,qc_align,dedup,variants,logs,tmp_picard}
# se crea directorio con la variable OUTDIR (dentro del disco externo) y los subdirectorios
# tmp_picard: carpeta de archivos intermedios para Picard MarkDuplicates. Por defecto
# Picard/Java usan el /tmp del sistema (normalmente en el disco interno, con poco
# espacio libre) para los ficheros temporales de ordenación; con BAMs de WGS grandes
# esto puede llenar /tmp y provocar errores de "No space left on device" aunque el
# disco externo de resultados tenga espacio de sobra. Al indicarle TMP_DIR más abajo,
# Picard usa esta carpeta del disco externo en su lugar.

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$OUTDIR/logs/pipeline.log"; }
# se crea la función log para registrar las salidas del pipeline, el día y la hora

log "FASTQ R1 resueltos: $FASTQ_R1"
log "FASTQ R2 resueltos: $FASTQ_R2"
log "Carpeta de resultados: $OUTDIR"


# --------------- 0.c VALIDACIÓN TEMPRANA DE ENTRADAS -------------------------
# Comprobar que el archivo de genoma de referencia y los archivos asocidados necesarios para el alineamiento
# esten disponibles -f
if [[ ! -f "$REF_GENOME" ]]; then
    log "ERROR: no se encuentra el genoma de referencia: $REF_GENOME"
    exit 1
fi
if [[ ! -f "${REF_GENOME}.fai" ]]; then
    log "ERROR: falta el índice ${REF_GENOME}.fai (ejecutar: samtools faidx $REF_GENOME)"
    exit 1
fi
for ext in amb ann bwt pac sa; do
    if [[ ! -f "${REF_GENOME}.${ext}" ]]; then
        log "ERROR: falta el índice de BWA ${REF_GENOME}.${ext} (ejecutar: bwa index $REF_GENOME)"
        exit 1
    fi
done

IFS=',' read -ra _r1_check <<< "$FASTQ_R1"
IFS=',' read -ra _r2_check <<< "$FASTQ_R2"
for f in "${_r1_check[@]}" "${_r2_check[@]}"; do
    if [[ ! -f "$f" ]]; then
        log "ERROR: archivo FASTQ no encontrado: $f"
        exit 1
    fi
done


# --------------- 0.d EJECUCIÓN EN SEGUNDO PLANO (opción -b) -----------------
# Si se pidió -b, el script se relanza a sí mismo con nohup y en background,
# y el proceso original (el que ve el usuario en su terminal) termina de
# inmediato. Así, aunque se cierre la terminal o se corte la conexión SSH,
# el pipeline real sigue corriendo en el servidor.
#
# Se usa "$(realpath "$0")" en vez de "$0" a secas: si el script se invocó
# con una ruta relativa (./wgs_pipeline.sh) y el directorio de trabajo
# cambiase antes de que nohup termine de lanzar el proceso, "$0" podría
# dejar de apuntar al script correctamente.
#
# La variable de entorno WGS_PIPELINE_BG evita que el proceso relanzado
# vuelva a relanzarse a sí mismo en un bucle infinito.
if $BACKGROUND && [[ "${WGS_PIPELINE_BG:-0}" != "1" ]]; then
    SCRIPT_PATH="$(realpath "$0")"
    BG_LOG="$OUTDIR/logs/background_nohup.log"
    log "Flag -b activado: relanzando el pipeline en segundo plano con nohup"
    log "Puedes cerrar la terminal sin problema. Progreso en: $BG_LOG"
    WGS_PIPELINE_BG=1 nohup "$SCRIPT_PATH" \
         -r "$REF_GENOME" -1 "$FASTQ_R1" -2 "$FASTQ_R2" -s "$SAMPLE_ID" \
        -t "$THREADS" -o "$OUTDIR" -q "$SUBSAMPLE_READS" \
        > "$BG_LOG" 2>&1 < /dev/null &
    disown
    log "Proceso en segundo plano lanzado (PID $!)"
    echo "PID: $!"
    echo "Sigue el progreso con: tail -f \"$BG_LOG\""
    exit 0
fi


# ------------------------ 1. CARGA / VARIABLES -------------------------------

log "== Iniciando pipeline para muestra: $SAMPLE_ID =="

# ------------------------ 1.b FUSIÓN DE FASTQ (multi-lane) ------------------
# Si la muestra viene repartida en varios archivos (uno por carril de
# secuenciación), se fusionan aquí en un único R1 y un único R2 antes de
# cualquier análisis. Los FASTQ comprimidos con gzip se pueden concatenar
# directamente con "cat" sin descomprimir: el resultado es un gzip válido
# "multi-miembro" que cualquier herramienta (fastqc, fastp, bwa) lee sin
# problema, exactamente igual que si fuera un único archivo.
#

merge_fastqs() {
    local input_list="$1"
    local output_file="$2"

    IFS=',' read -ra files <<< "$input_list"

    if [[ ${#files[@]} -eq 1 ]]; then
        log "  $(basename "${files[0]}"): un único FASTQ, no se fusiona"
        echo "${files[0]}"
    else
        log "  $(basename "$output_file"): fusionando ${#files[@]} archivos (lanes)"
        cat "${files[@]}" > "$output_file"
        echo "$output_file"
    fi
}



R1_MERGED="$OUTDIR/merged/${SAMPLE_ID}_R1.merged.fastq.gz"
R2_MERGED="$OUTDIR/merged/${SAMPLE_ID}_R2.merged.fastq.gz"

log "Fusionando FASTQ de entrada en un único R1/R2 por muestra"
merge_fastqs "$FASTQ_R1" "$R1_MERGED"
merge_fastqs "$FASTQ_R2" "$R2_MERGED"

R1="$R1_MERGED"
R2="$R2_MERGED"

# ------------------------ 2. ANÁLISIS DE CALIDAD (FastQC) -------------------

# Función para correr el analisis de calidad (fastqc)
# se definen las variables locales r1 r2 y outdir
#
# Antes de invocar a FastQC, si SUBSAMPLE_READS > 0 se genera una copia
# submuestreada de R1 y R2 (las primeras SUBSAMPLE_READS lecturas de cada
# uno) dentro de $OUTDIR/tmp_fastqc_subsample, lanzando R1 y R2 en paralelo
# (en background con &, sincronizados con wait) para aprovechar 2 hilos en
# vez de leer los FASTQ completos en serie. FastQC se ejecuta después sobre
# esas copias pequeñas, mucho más rápido, y las copias se borran al acabar.
run_fastqc() {
    local r1=$1 r2=$2 outdir=$3
    local fastqc_r1="$r1"
    local fastqc_r2="$r2"
    local subdir=""

    if [[ "$SUBSAMPLE_READS" -gt 0 ]]; then
        subdir="$OUTDIR/tmp_fastqc_subsample"
        mkdir -p "$subdir"
        rm -f "$subdir"/*
        local nlines=$((SUBSAMPLE_READS * 4))
        fastqc_r1="$subdir/$(basename "$r1")"
        fastqc_r2="$subdir/$(basename "$r2")"

        log "Submuestreando $SUBSAMPLE_READS reads de $r1 y $r2 para FastQC (en paralelo)"
        zcat "$r1" | head -n "$nlines" | gzip > "$fastqc_r1" &
        zcat "$r2" | head -n "$nlines" | gzip > "$fastqc_r2" &
        wait
    fi

    log "Ejecutando FastQC sobre $fastqc_r1 y $fastqc_r2"
    fastqc -t "$THREADS" -o "$outdir" "$fastqc_r1" "$fastqc_r2" &> "$OUTDIR/logs/fastqc_${outdir##*/}.log"

    [[ -n "$subdir" ]] && rm -rf "$subdir"
}

# &> manda los mensajes de progreso a un archivo en la carpeta logs

# Se crea funcion qc_pasa_criterios
# Comprueba si FastQC ha marcado FAIL en alguno de los módulos definidos como críticos
# en MODULOS_CRITICOS (ver configuración arriba). Módulos no listados ahí no bloquean
# el pipeline aunque FastQC los marque como FAIL o WARN.
qc_pasa_criterios() {
    local outdir=$1
    local fail=0
    # Se usa un archivo temporal ÚNICO por llamada (mktemp), no un nombre
    # fijo /tmp/summary.txt: con -b activo puede haber varias ejecuciones
    # del pipeline corriendo a la vez (distintas muestras), y un nombre fijo
    # en /tmp haría que se pisaran el archivo entre sí de forma intermitente.
    local tmp_summary
    tmp_summary="$(mktemp)"
    trap 'rm -f "$tmp_summary"' RETURN

    for zip in "$outdir"/*_fastqc.zip; do
        unzip -p "$zip" "*/summary.txt" > "$tmp_summary"
        for modulo in "${MODULOS_CRITICOS[@]}"; do
            if grep -qP "^FAIL\t${modulo}\t" "$tmp_summary"; then
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
# (Los índices de $REF_GENOME ya se validaron en el paso 0.c, antes de
# gastar tiempo en QC/fastp, así que aquí no debería fallar por eso.)

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

log "Marcando duplicados PCR con Picard MarkDuplicates (TMP_DIR=$OUTDIR/tmp_picard)"
picard MarkDuplicates \
    I="$BAM_RAW" \
    O="$BAM_DEDUP" \
    M="$OUTDIR/dedup/${SAMPLE_ID}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    TMP_DIR="$OUTDIR/tmp_picard" \
    &> "$OUTDIR/logs/picard_markdup.log"

# Limpieza de los temporales de Picard una vez terminado, para no dejar
# ocupado el disco externo innecesariamente entre ejecuciones.
rm -rf "${OUTDIR:?}/tmp_picard"/*

samtools index "$BAM_DEDUP"

# ------------------------ 7. LLAMADO DE VARIANTES (bcftools mpileup+call) ---
# Se sustituye GATK HaplotypeCaller por bcftools mpileup + bcftools call.
# Es un llamador más ligero (menos dependencias, más rápido) y suficiente
# para variantes germinales individuales; GATK suele preferirse cuando se
# necesita el procesado conjunto de cohortes (GVCF + GenotypeGVCFs) o los
# pasos de recalibración de variantes (VQSR), que aquí no se usan.
#
#   bcftools mpileup: calcula las probabilidades genotípicas por posición
#                      a partir del BAM (equivalente al primer paso de GATK).
#   bcftools call -mv: aplica el modelo multialélico (-m) y se queda solo
#                      con posiciones variantes (-v), no con todo el genoma.
#   -a AD,DP:          añade profundidad total (DP) y por alelo (AD) al VCF,
#                      útil para filtrado posterior y para el informe.

VCF_OUT="$OUTDIR/variants/${SAMPLE_ID}.vcf.gz"
log "Llamando variantes germinales con bcftools mpileup + call"
bcftools mpileup \
    -f "$REF_GENOME" \
    -a AD,DP \
    --threads "$THREADS" \
    -Ou \
    "$BAM_DEDUP" \
    2> "$OUTDIR/logs/bcftools_mpileup.log" \
    | bcftools call \
        -mv \
        --threads "$THREADS" \
        -Oz \
        -o "$VCF_OUT" \
        2> "$OUTDIR/logs/bcftools_call.log"

tabix -p vcf "$VCF_OUT"

log "== Pipeline completado. VCF final: $VCF_OUT =="
log "Todos los resultados quedan en el disco externo: $OUTDIR"
