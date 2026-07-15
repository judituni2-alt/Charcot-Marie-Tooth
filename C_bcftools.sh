#!/usr/bin/env bash


log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$OUTDIR/logs/pipeline.log"; }

REF_GENOME="hg38_v0_Homo_sapiens_assembly38.fasta"
SAMPLE_ID="TEST"
THREADS=8              # Son los hilos, es decir, unidades de ejecución paralela. Valor por defecto: 8 (evaluar cuantos poner según el ordenador que vaya a utilizar)
OUTDIR="resultados"
BAM_DEDUP="$OUTDIR/dedup/${SAMPLE_ID}.dedup.bam"


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
