#!/usr/bin/env bash
################################################################################
# B_fastqc_merge.sh
#
# Versión reducida del pipeline WGS: solo incluye los pasos de
# FastQC y fusión (merge) de FASTQ, EN ESTE ORDEN:
#
#   1) FastQC sobre los FASTQ de entrada (uno por carril/lane, tal cual llegan)
#   2) Fusión (merge) de esos FASTQ en un único R1 y un único R2 por muestra
#
# Los archivos FASTQ de entrada se leen desde un disco duro/unidad externa
# y todos los archivos y carpetas que genera el script se guardan también
# en ese mismo disco.
#
# Uso:
#   ./B_fastqc_merge.sh
#       -1 CMT1234_R1.fastq.gz -2 CMT1234_R2.fastq.gz -s CMT1234 -o resultados
#
#   Si la muestra viene repartida en varios archivos por carril de
#   secuenciación (lane splitting, típico de Illumina: L001, L002, L003...),
#   se pueden indicar varios archivos separados por coma en -1 y -2:
#
#   ./B_fastqc_merge.sh
#       -1 CMT1234_L001_R1.fastq.gz,CMT1234_L002_R1.fastq.gz \
#       -2 CMT1234_L001_R2.fastq.gz,CMT1234_L002_R2.fastq.gz \
#       -s CMT1234 -o resultados
#
# Requiere instalar: fastqc
# Para ello crear un ambiente de conda/mamba con los recursos necesarios
###############################################################################

set -euo pipefail
# activa 3 procedimientos de seguridad en caso de que el script falle
#  -e: si cualquier comando del script devuelve un error (exit status distinto de 0) el script se detiene inmediatamente
#  -u: si se usa una variable que no ha sido definida previamente el script falla
#  -o pipefail: cuando se encadena con tuberías, si cualquiera de los comandos de la cadena falla, toda la tubería
#   se considera fallida


# ---------------------------- 0. CONFIGURACIÓN ------------------------------
# Se crean las variables vacías que se necesitan

SAMPLE_ID=""           # ID de la muestra
FASTQ_R1=""            # uno o varios archivos separados por coma (lanes)
FASTQ_R2=""            # uno o varios archivos separados por coma (lanes)
THREADS=8              # Son los hilos, es decir, unidades de ejecución paralela. Valor por defecto: 8
OUTDIR=""              # se resuelve más abajo, tras parsear -s

# se crea la funcion de uso
usage() {
    cat <<EOF
Uso: $0  -1 R1.fastq.gz -2 R2.fastq.gz -s sample_id [-t threads] [-o outdir]

  -1 / -2   Admiten varios archivos separados por coma (sin espacios) si la
            muestra está repartida en varios carriles/lanes, ej:
            -1 L001_R1.fastq.gz,L002_R1.fastq.gz
EOF
    exit 1
}

while getopts "1:2:s:t:o:h" opt; do
    case $opt in
        1) FASTQ_R1="$OPTARG" ;;
        2) FASTQ_R2="$OPTARG" ;;
        s) SAMPLE_ID="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;                  # si se detecta -h (help) como flag en el comando inicial, imprime usage
        *) usage ;;                  # si se detecta una flag no definida *, imprime usage
    esac
done

[[ -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$SAMPLE_ID" ]] && usage
# comprueba si las variables FASTQ_R1 o FASTQ_R2 o SAMPLE_ID estan vacias (-z es TRUE) y si lo están ejecuta usage

mkdir -p "$OUTDIR"/{qc_raw,merged,logs}
# se crea directorio con la variable OUTDIR (dentro del disco externo) y los subdirectorios necesarios

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$OUTDIR/logs/pipeline.log"; }
# se crea la función log para registrar las salidas del pipeline, el día y la hora

log "FASTQ R1 resueltos: $FASTQ_R1"
log "FASTQ R2 resueltos: $FASTQ_R2"
log "Carpeta de resultados: $OUTDIR"


# --------------- 0.b VALIDACIÓN TEMPRANA DE ENTRADAS -------------------------
IFS=',' read -ra _r1_check <<< "$FASTQ_R1"
IFS=',' read -ra _r2_check <<< "$FASTQ_R2"
for f in "${_r1_check[@]}" "${_r2_check[@]}"; do
    if [[ ! -f "$f" ]]; then
        log "ERROR: archivo FASTQ no encontrado: $f"
        exit 1
    fi
done


# ------------------------ 1. CARGA / VARIABLES -------------------------------

log "== Iniciando FastQC + merge para muestra: $SAMPLE_ID =="

# ------------------------ 2. ANÁLISIS DE CALIDAD (FastQC) -------------------
# Se ejecuta FastQC sobre TODOS los archivos de entrada tal cual llegan
# (un archivo por carril/lane si la muestra viene repartida), ANTES de
# fusionarlos. Así se puede inspeccionar la calidad de cada lane por separado.

run_fastqc() {
    local outdir=$1
    shift
    local files=("$@")
    log "Ejecutando FastQC sobre: ${files[*]}"
    fastqc -t "$THREADS" -o "$outdir" "${files[@]}" &> "$OUTDIR/logs/fastqc_${outdir##*/}.log"
}
# &> manda los mensajes de progreso a un archivo en la carpeta logs

qc_dir="$OUTDIR/qc_raw"
run_fastqc "$qc_dir" "${_r1_check[@]}" "${_r2_check[@]}"

log "FastQC completado. Resultados en: $qc_dir"


# ------------------------ 3. FUSIÓN DE FASTQ (multi-lane) -------------------
# Si la muestra viene repartida en varios archivos (uno por carril de
# secuenciación), se fusionan aquí en un único R1 y un único R2. Los FASTQ
# comprimidos con gzip se pueden concatenar directamente con "cat" sin
# descomprimir: el resultado es un gzip válido "multi-miembro" que cualquier
# herramienta (fastqc, fastp, bwa) lee sin problema, exactamente igual que
# si fuera un único archivo.

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

log "== Completado. R1 fusionado: $R1 =="
log "== Completado. R2 fusionado: $R2 =="
