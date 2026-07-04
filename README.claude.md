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

Opción 1: Utilizar archivos de genoma de referencia ya indexado extraído de 
https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
Archivos:
 **Homo_sapiens_assembly38.fasta** --> El genoma de referencia completo en formato FASTA. Archivo base contra el que se alinean las lecturas y se llaman variantes.
 **Homo_sapiens_assembly38.fasta.fai** --> Índice de samtools. Archivo que guarda, para cada cromosoma, dónde empieza y termina dentro del fasta, y cuántos caracteres tiene cada línea. Permite acceder rápidamente a cualquier posición del genoma sin tener que leer el archivo entero.
 **Homo_sapiens_assembly38.dict** --> Diccionario de secuencias formato SAM. Contiene lista de los cromosomas del genoma, su longitud y un checksum (MD5) de cada uno. GATK lo usa para verificar que el FASTA, el BAM y los archivos VCF sean todos compatibles entre sí (mismo genoma, mismas coordenadas).
 **Archivos para indice de alineamiento**
   - .64.amb: Guarda información sobre las regiones del genoma con bases ambiguas (N), es decir, huecos o zonas no secuenciadas con certeza.
   - .64.ann: Contiene anotaciones generales del genoma: nombres de los cromosomas, sus longitudes, y dónde empieza cada uno dentro del índice.
   - .64.bwt: La transformada de Burrows-Wheeler del genoma. Reordenación matemática de toda la secuencia que permite búsquedas de texto extremadamente rápidas.
   - .64.pac: El genoma comprimido en formato binario (2 bits por base en vez de 1 byte), para ahorrar espacio y acelerar el acceso.
   - .64.sa:  El "suffix array". Estructura de datos que, combinada con el .bwt, permite encontrar la posición exacta en el genoma de cualquier fragmento de secuencia de forma casi instantánea.
 **Contigs alternativos**
   - .64.alt: Lista de las regiones donde existe más de una versión conocida de la secuencia (por ejemplo, la región HLA, muy variable entre personas). Le dice a BWA que esas lecturas pueden alinear en dos sitios distintos legítimamente, evitando que las marque como error o las descarte.

Nota: preferible utilizar los archivos .64


Opción 2: Usar código para indexar otro genoma de referencia de forma que se generen estos archivos a partir de los 
archivos .fasta .fai y .dict
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
