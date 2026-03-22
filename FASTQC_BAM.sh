

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

# TEORIA
# Lenguaje WDL: Workflow description language es un lenguaje diseñado para definir pipelines de analisis biologicos.
# Necesita un motor de ejecucion como Cromwell o MiniWDL para ejecurarlo desde bash


# 2.  Activamos el entorno conda que acabamos de crear
conda activate FASTQC_BAM 


# 3. En esta primera parte de control de calidad y filtrado, vamos a comenzar creando las carpetas que necesitaremos en un primer momento:

mkdir -p Quality/Raw Quality/Filtered Trimmed


# 5. Realizamos el control de calidad con FastQC:

fastqc *fastq.gz -o Quality/Raw/ -t 32


# 6. Ejecutamos el filtrado de secuencias por calidad y longitud

### En este caso, lo podemos hacer o bien de forma individual, cuando tenemos pocas muestras: 

 fastp --in1 Muestra1_R1* --in2 Muestra1_R2* --out1 Trimmed/Muestra1_R1_filtered.fastq.gz --out2 Trimmed/Muestra1_R2_filtered.fastq.gz --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 20 --length_required 35 --json Trimmed/Muestra1.json --html Trimmed/Muestra1.html --thread 32

### O mediante un bucle cuando son muchas, en este caso habría que incluir todos los nombres de las muestras en un fichero, por ejemplo muestras.txt, un identificador por línea y ejecutar el bucle:

ls *fastq.gz | cut -d _ -f 1 | sort -u > muestras.txt
for i in $(cat muestras.txt); do fastp --in1 $i*R1* --in2 $i*R2* --out1 Trimmed/$i"_R1_filtered.fastq.gz" --out2 Trimmed/$i"_R2_filtered.fastq.gz" --detect_adapter_for_pe --cut_front --cut_tail --cut_window_size 12 --cut_mean_quality 30 --length_required 35 --json Trimmed/$i.json --html Trimmed/$i.html --thread 32; done


# 7. Ejecutamos de nuevo un segundo control de calidad para corroborar que hemos realizado correctamente el filtrado

fastqc Trimmed/*fastq.gz -o Quality/Filtered/ --threads 32


# 8. Podemos generar un informe resumen de los pasos que hemos realizado con multiqc:

multiqc .

 



