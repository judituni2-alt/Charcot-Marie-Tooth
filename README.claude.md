# Pipeline WGS - Línea germinal

Implementa el flujo del diagrama: **FASTQ → QC → Alineamiento → QC BAM → Eliminar duplicados PCR → Llamado de variantes → VCF**

## Estructura de archivos

| Archivo | Lenguaje | Función |
|---|---|---|
| `wgs_pipeline.sh` | bash | Orquesta todo el pipeline, paso a paso, con los bucles de decisión del diagrama |
| `plot_qc_metrics.R` | R | Genera gráficos de las métricas de calidad (alineamiento, duplicados) |
| `summarize_vcf.py` | Python | Resume el VCF final (nº variantes, SNPs vs indels, calidad, por cromosoma) |

## 1. Instalación de dependencias

Recomendado con conda/mamba:

```bash
conda create -n wgs_pipeline -c bioconda -c conda-forge \
    fastqc fastp bwa samtools bamtools picard gatk4 \
    r-ggplot2 r-dplyr r-stringr
conda activate wgs_pipeline

pip install cyvcf2 pandas --break-system-packages
```

## 2. Preparar el genoma de referencia

```bash
bwa index referencia.fasta
samtools faidx referencia.fasta
picard CreateSequenceDictionary R=referencia.fasta O=referencia.dict
```

## 3. Ejecutar el pipeline

```bash
chmod +x wgs_pipeline.sh
./wgs_pipeline.sh \
    -r referencia.fasta \
    -1 muestra_R1.fastq.gz \
    -2 muestra_R2.fastq.gz \
    -s muestra01 \
    -t 8 \
    -o resultados
```

Esto genera dentro de `resultados/`:

```
resultados/
├── qc_raw/          # FastQC de las lecturas crudas (y recortadas si aplicó)
├── trimmed/          # Lecturas recortadas con fastp (solo si el QC inicial falló)
├── align/            # BAM ordenado
├── qc_align/         # bamtools stats + samtools flagstat
├── dedup/             # BAM sin duplicados + métricas de Picard
├── variants/          # VCF final
└── logs/              # log de cada paso
```

### Lógica de las decisiones del diagrama

- **¿Cumple criterios de calidad? (FASTQ)** → el script comprueba si FastQC marca algún módulo como `FAIL`. Si es así, recorta con `fastp` y repite el QC (hasta 2 rondas antes de abortar).
- **¿Cumple criterios de calidad? (BAM)** → comprueba el % de lecturas mapeadas con `samtools flagstat` (umbral por defecto: 90%). Si no lo cumple, el pipeline avisa pero continúa (puedes ajustar esto a "abortar" si prefieres ser más estricto).

## 4. Visualizar la calidad (R)

```bash
Rscript plot_qc_metrics.R resultados muestra01
```

Genera:
- `resultados/qc_align/muestra01_align_metrics.png`
- `resultados/dedup/muestra01_duplication.png`

## 5. Resumir el VCF final (Python)

```bash
python summarize_vcf.py resultados/variants/muestra01.g.vcf.gz resumen_muestra01.csv
```

## Notas sobre las alternativas del diagrama

El diagrama muestra varias herramientas posibles en cada paso. El script usa las más estándar actualmente (BWA-MEM, GATK HaplotypeCaller), pero puedes sustituirlas fácilmente:

- **Alineamiento**: cambia la función `run_alineamiento()` por `bowtie2`, `soap`, o `mosaik` si tu proyecto lo requiere específicamente.
- **Llamado de variantes**: sustituye el bloque de `gatk HaplotypeCaller` por `freebayes`, `platypus`, o `bcftools mpileup | bcftools call` si tu TFM compara herramientas de variant calling.

## Siguiente paso sugerido

Si tu TFM requiere comparar varias herramientas (p. ej. varios alineadores o varios variant callers), lo más limpio es parametrizar esas dos funciones del bash con un flag (`-a bwa|bowtie2` y `-v gatk|freebayes`) para poder correr el mismo pipeline con distintas combinaciones y comparar resultados.
