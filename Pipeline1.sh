

#!/bin/bash

# Tener instalado una distribución de conda, ya sea Anaconda o Miniconda 
# Crear un entorno con los programas que vamos a utilizar.

conda create -n Pipeline1 -c defaults -c conda-forge -c bioconda fastp multiqc bwa samtools qualimap htslib openjdk=17 quast bandage

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
# Opciones de paquetes para ensamblado dependiendo de la tecnica
# Ensambladores recomendados
#   HiFi reads → Hifiasm (muy usado para genomas humanos completos).
#   Long reads → Canu, Flye, Shasta.
#   Los ensambladores híbridos ya son menos necesarios con PacBio HiFi, pero MaSuRCA o HybridSPAdes pueden usarse si se combinan short + long reads.
# Pulido y corrección
#   Pilon → corrige errores residuales usando Illumina.
#   Racon / Medaka → corrige errores específicos de Nanopore.
# Evaluación de calidad del ensamblado
#   QUAST → estadísticas de ensamblaje (N50, contigs, cobertura).
#   BUSCO → comprueba genes esenciales conservados para evaluar completitud
# htslib: Es una biblioteca en C que proporciona funciones para trabajar con archivos SAM/BAM/CRAM/VCF/BCF.
# Es la base de herramientas como Samtools y Bcftools.
# openjdk=17 es el Java Development Kit de código abierto. Necesario para trabajar con paquetes como ???
# bandage: para visualizar ensamblados



# TEORIA
# Lenguaje WDL: Workflow description language es un lenguaje diseñado para definir pipelines de analisis biologicos.
# Necesita un motor de ejecucion como Cromwell o MiniWDL para ejecurarlo desde bash


# 3. Activamos el entorno conda que acabamos de crear
conda activate


# 4. En esta primera parte de control de calidad y filtrado, vamos a comenzar creando las carpetas que necesitaremos en un primer momento:

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ mkdir -p Quality/Raw Quality/Filtered Trimmed


# 5. Realizamos el control de calidad con FastQC:

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fastqc --help

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fastqc *fastq.gz -o Quality/Raw/ -t 32


# 6. Ejecutamos el filtrado de secuencias por calidad y longitud

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fastp --help

### En este caso, lo podemos hacer o bien de forma individual, cuando tenemos pocas muestras: 

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fastp --in1 Muestra1_R1* --in2 Muestra1_R2* --out1 Trimmed/Muestra1_R1_filtered.fastq.gz --out2 Trimmed/Muestra1_R2_filtered.fastq.gz --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 20 --length_required 35 --json Trimmed/Muestra1.json --html Trimmed/Muestra1.html --thread 32

### O mediante un bucle cuando son muchas, en este caso habría que incluir todos los nombres de las muestras en un fichero, por ejemplo muestras.txt, un identificador por línea y ejecutar el bucle:

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ ls *fastq.gz | cut -d _ -f 1 | sort -u > muestras.txt

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ for i in $(cat muestras.txt); do fastp --in1 $i*R1* --in2 $i*R2* --out1 Trimmed/$i"_R1_filtered.fastq.gz" --out2 Trimmed/$i"_R2_filtered.fastq.gz" --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 30 --length_required 35 --json Trimmed/$i.json --html Trimmed/$i.html --thread 32; done


# 7. Ejecutamos de nuevo un segundo control de calidad para corroborar que hemos realizado correctamente el filtrado

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fastqc Trimmed/*fastq.gz -o Quality/Filtered/ --threads 32


# 8. Podemos generar un informe resumen de los pasos que hemos realizado con multiqc:

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ multiqc .


# BONUS. Además, podemos reproducir los pasos de la práctica de este tema con otras secuencias que podemos descargarnos. En mi caso he ido al NCBI y he seguido los siguientes pasos:

## 9.1. He buscado Neisseria gonorrhoeae y clicado en la sección de 'SRA' y seleccionado el primer resultado (https://www.ncbi.nlm.nih.gov/sra/SRX31029457[accn]). En esta página nos indica un SRR** que es el identificador de los fastq de la muestra. He buscado una segunda muestra, con el mismo procedimiento de Neisseria meningitidis. 

## 9.2. Con los identificadores SRR** realizamos la descarga, la forma más rápida es mediante la herramienta sra-tools. 

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ conda create -n sra -c bioconda -c conda-forge -c defaults -c r sra-tools

(Tema4) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ conda activate sra

(sra) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fasterq-dump --help

### Lo podemos hacer o bien de forma individual, cuando tenemos pocas muestras: 

(sra) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fasterq-dump SRR35980167 --split-files -e 32 -p -O .
(sra) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ fasterq-dump SRR35704398 --split-files -e 32 -p -O .

### O mediante un bucle cuando son muchas, en este caso habría que incluir todos los SRR** en un fichero, por ejemplo SRRs.txt, un identificador por línea y ejecutar el bucle:

(sra) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ for i in $(cat SRRs.txt); do fasterq-dump $i --split-files -e 32 -p -O .; done

## 9.3. Opcionalmente, podemos comprimir los fastq para ahorrar espacio de almacenamiento:

(sra) laura@laura-intel-nox:/media/8TB/UNIR/Tema4$ gzip *fastq
