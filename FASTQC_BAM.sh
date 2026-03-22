

#!/bin/bash

# Tener instalado una distribución de conda, ya sea Anaconda o Miniconda 
# 1. Crear un entorno con los programas que vamos a utilizar.

conda create -n FASTQC_BAM -c defaults -c conda-forge -c bioconda fastqc fastp multiqc bwa samtools qualimap htslib openjdk=17

# -n nombre del entorno
# -c se usa para indicar el canal de conda desde donde instalar los paquetes requeridos. Un canal de conda
# es un repositorio de paquetes.
# defaults: canal oficial de Anaconda, contiene los paquetes mas estables que distribuye anaconda
# conda-forge: canal comunitario de conda. Canal del que depence Bioconda
# Bioconda: canal de conda especializado en bioinformatica que contiene distintos paquetes/herramientas para
# el análisis de secuencias
# fastp: paquete para realizar control de calidad (QC) de las secuencias
# multiqc: para combinar todos los reportes generados con fastp
# bwa: paquete para alineamiento con genoma de referencia
# samtools: trabajar con los ficheros del mapeo SAM/BAM
# Qualimap: evaluar la calidad del mapeo

# CONTROL DE CALIDAD

# 2.  Activamos el entorno conda que acabamos de crear
conda activate FASTQC_BAM


# 3. En esta primera parte de control de calidad y filtrado, vamos a comenzar creando las carpetas que necesitaremos en un primer momento:

mkdir -p Quality/Raw Quality/Filtered Trimmed


# 5. Realizamos el control de calidad con FastQC:

fastqc *fastq.gz -o Quality/Raw/ -t 32


# 6. Ejecutamos el filtrado de secuencias por calidad y longitud
# Se realiza mediante bucle, en este caso habría que incluir todos los nombres de las muestras en un fichero, por ejemplo muestras.txt, un identificador por línea y ejecutar el bucle:

ls *fastq.gz | cut -d _ -f 1 | sort -u > muestras.txt
for i in $(cat muestras.txt); do fastp --in1 $i*R1* --in2 $i*R2* --out1 Trimmed/$i"_R1_filtered.fastq.gz" --out2 Trimmed/$i"_R2_filtered.fastq.gz" --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 30 --length_required 35 --json Trimmed/$i.json --html Trimmed/$i.html --thread 32; done


# 7. Ejecutamos de nuevo un segundo control de calidad para corroborar que hemos realizado correctamente el filtrado

fastqc Trimmed/*fastq.gz -o Quality/Filtered/ --threads 32


# 8. Podemos generar un informe resumen de los pasos que hemos realizado con multiqc:

multiqc .

 



# ALINEAMIENTO

# 1. Otro posible paso es el mapeo de las muestras contra un genoma de referencia.
# Para ello, vamos a enfrentar todas las lecturas de nuestros fastq para ver en qué ubicación exacta se localizan y poder tener información sobre la cobertura y profundidad de lectura.

### 1.2. La mayoría de programas dedicados al mapeo de lecturas requiere que indexemos el genoma de referencia antes de ejecutarlo.
bwa index Reference/NC_045512.2.fasta # Primero indexamos el genoma


## 1.3. Una vez indexada la referencia, enfrentamos nuestros ficheros a la referencia.

mkdir -p Mapped/Filtered
bwa mem -Y -M -t 32 -o Mapped/Muestra1.sam Reference/NC_045512.2.fasta Trimmed/Muestra1_R1_filtered.fastq.gz Trimmed/Muestra1_R2_filtered.fastq.gz
for i in $(cat muestras.txt); do bwa mem -Y -M -t 32 -o Mapped/$i".sam" Reference/NC_045512.2.fasta Trimmed/$i*R1* Trimmed/$i*R2*; done


## 1.4. El mapeo nos va a generar unos ficheros SAM donde se almacena la información del análisis, no obstante suelen ser ficheros pesados, desestructurados y conviene convertirlos a formato BAM para mejor eficiencia (visual y analitica)
samtools view -Sb Mapped/Muestra1.sam --threads 32 -o Mapped/Muestra1".bam"
for i in $(cat muestras.txt); do samtools view -Sb Mapped/$i.sam --threads 32 -o Mapped/$i".bam"; done 

## Además, se recomienda ordenarlo con la función 'sort' de samtools para ordenar las lecturas que sean similares, justamente para esa optimización computacional.

samtools sort Mapped/Muestra1*bam -o Mapped/Filtered/Muestra1"_sorted.bam"
for i in $(cat muestras.txt); do samtools sort Mapped/$i*bam -o Mapped/Filtered/$i"_sorted.bam"; done


## 1.5. Acto seguido podemos evaluar la calidad del mapeo, para verificar que es correcto.
qualimap bamqc -bam Mapped/Filtered/Muestra1_sorted.bam -outdir Mapped/Filtered/Quality/Muestra1
for i in $(cat muestras.txt); do qualimap bamqc -bam Mapped/Filtered/$i"_sorted.bam" -outdir Mapped/Filtered/Quality/$i; done
