# PIPELINE VARIANTES INTRÓNICAS CHARCOT-MARIE-TOOTH
El pipeline sigue el siguiente flujo de trabajo

Preprocesamiento

<img width="342" height="256" alt="TFM" src="https://github.com/user-attachments/assets/d75167e5-d971-4a8c-9ce6-04bfe6bbda21" />

Filtrado variantes candidatas

<img width="342" height="256" alt="TFM (2)" src="https://github.com/user-attachments/assets/f0b9739d-04f7-452c-9d7e-62245d3c1e6e" />

## Archivos necesarios para realizar el Preprocesamiento
- Archivos en formato FASTQ paired end  
- Genoma de referencia:  https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0  
- Indices del BWA  
.fasta.amb  
.fasta.ann  
.fasta.bwt  
.fasta.pac  
.fasta.sa  
Se pueden descargar o generar a partir del genoma de referencia
```sh
bwa index Homo_sapiens_assembly38.fasta
```  
- Indice
.fasta.fai  
Se pueden descargar o generar a partir del genoma de referencia
```sh
samtools faidx Homo_sapiens_assembly38.fasta
```  
- Diccionario
.fasta.dict  
Se pueden descargar o generar a partir del genoma de referencia
```sh
picard CreateSequenceDictionary \
    R=Homo_sapiens_assembly38.fasta \
    O=Homo_sapiens_assembly38.dict
``` 


## Pasos para realizar Preprocesamiento
```sh
conda create -n Preprocesamiento -c bioconda -c conda-forge fastqc fastp bwa samtools bamtools picard bcftools tabix pandas
conda activate Preprocesamiento
```
 Una vez en el entorno, ejecutar el script con los argumentos necesarios  
 
 Listado de argumentos obligatorios y opcionales:
  
  -1 / -2   Obligatorio.       Archivos formato FASTQ paired end (-1 R1 -2 R2)  
  -r        Obligatorio.       Genoma de referencia en formato FASTA  
  -s        Obligatorio.       Identificador de la muestra  
  -o        Obligatorio.       Nombre del directorio donde se almacenan resultados.  
  -t        Opcional.          Hilos de procesador que se quieren utilizar en las herramientas que lo permitan. Por defecto 8  
  -q        Opcional.          Nº de reads a submuestrear por archivo antes de FastQC  (por defecto 2000000). Usar -q 0 para analizar el FASTQ completo  
	-b        Opcional.          Para activar la función background y evitar que se detenga el proceso si se cierra la terminal
  
```sh
 ./A_Preprocesamiento_pipeline_v7.sh -1 CMT1234_R1.fastq.gz -2 CMT1234_R2.fastq.gz -r ref.fasta -r ref.fasta -o resultados_CMT1234
```
Si la muestra viene repartida en varios archivos por carril de secuenciación (lane splitting, L001, L002, L003...),se pueden indicar varios archivos separados por coma en -1 y -2, y el pipeline los fusiona automáticamente en un único R1 y un único R2 antes de continuar:  
```sh
 ./A_Preprocesamiento_pipeline_v7.sh
-1 CMT1234_L001_R1.fastq.gz,CMT1234_L002_R1.fastq.gz \
-2 CMT1234_L001_R2.fastq.gz,CMT1234_L002_R2.fastq.gz \
-s CMT1234 \
-r ref.fasta \
-o resultados_CMT1234
```
Una vez ejecutado el pipeline se obtendrá el Archvio VCF crudo que contiene todas las variantes resultantes del Variant Calling

## Pasos para generar resumen numérico de Archivo VCF
```sh
conda create -n Resumen_vcf -c conda-forge -c bioconda pandas cyvcf2 
```
```sh
python B_Resumen_vcf_preprocesamiento_v3.py   --vcf /home/juditdvm/Charcot-Marie-Tooth/TEST.vcf.gz    --log-level INFO
```
Genera un resumen del VCF crudo final del pipeline de preprocesamiento:
 - nº de variantes totales
 - nº de SNPs vs indels
 - variantes por cromosoma

## Pasos para realizar el filtrado de variantes candidatas

Archivos necesarios para el filtrado de variantes
- VCF con todas las variantes (ya generado)
- Archivo BED genes de interes (en el repositorio)
- Archivo BED genes de interés regiones intrónicas+15pb hacia el exon (en repositorio)
- Genoma de referencia (ya descargado)
- Directorio con cache vep (crear)
- VCF gnomAD con frecuencia alelica
- Recurso SpliceAI para SNV
- Recurso SpliceAI para indels

  
Crear directorio con cache vep
```sh
mkdir -p /ruta/vep_cache
cd /ruta/vep_cache
wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_113_GRCh38.tar.gz
```
Crear entorno conda para ejecutar el script
```sh
conda create -n filtrado
```
Para realizar el filtrado ejecutar:
```sh
./F_pipeline_filtrado_v3.sh
-s paciente01 \
      -v paciente01.raw.vcf.gz \
      -g genes.bed \
      -i intrones_padded.bed \
      -r GRCh38.fasta \
      -c /ruta/.vep \
      -n gnomad.genomes.vcf.gz \
      -x spliceai_scores.raw.snv.hg38.vcf.gz \
      -y spliceai_scores.raw.indel.hg38.vcf.gz \
      -m 0.01 -t 8 -o resultados_paciente01
```
Para crear un tsv a partir del vfc generado
```
{
    echo -e "CHROM\tPOS\tREF\tALT\tAllele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL_SOURCE\tHGNC_ID\tMANE\tMANE_SELECT\tSpliceAI_pred_DP_AG\tSpliceAI_pred_DP_AL\tSpliceAI_pred_DP_DG\tSpliceAI_pred_DP_DL\tSpliceAI_pred_DS_AG\tSpliceAI_pred_DS_AL\tSpliceAI_pred_DS_DG\tSpliceAI_pred_DS_DL\tSpliceAI_pred_SYMBOL"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' input.vcf.gz |
awk -F'\t' 'BEGIN {OFS="\t"} {
    n=split($5,a,"|")
    printf "%s\t%s\t%s\t%s", $1,$2,$3,$4
    for(i=1;i<=n;i++)
        printf "\t%s",a[i]
    printf "\n"
}'
} > output.tsv
```
Para filtrar el archivo tsv de manera que se genere otro con variantes que tengan algún valor SpliceAI_pred_DS>0.2
```
python -c "import pandas as pd; df=pd.read_csv('input.tsv', sep='\t'); cols=['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL']; df[cols]=df[cols].apply(pd.to_numeric, errors='coerce'); df[df[cols].max(axis=1)>0.2].to_csv('spliceAI_gt_0.2.tsv', sep='\t', index=False)"
```

