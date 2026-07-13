# Pipeline: genes → BED (API Ensembl)

Genera un archivo BED de coordenadas génicas a partir de una lista de
nombres de genes en texto plano, usando la API REST de Ensembl
(https://rest.ensembl.org).

## Archivos

- `genes_to_bed.py` — script principal.
- `test_conexion.py` — comprobación rápida de conectividad con la API antes de correr el pipeline completo.
- `genes_ejemplo.txt` — ejemplo de archivo de entrada con genes de CMT.

## Requisitos

- Python 3.9+ (usa solo librería estándar: `urllib`, `json`, `argparse` — no requiere instalar nada con pip).
- Acceso a internet hacia `rest.ensembl.org` (o `grch37.rest.ensembl.org` si se utiliza 
 GRCh37).

## Formato de entrada

Archivo de texto plano, un gen por línea (símbolos HGNC estándar):

```
PMP22
GJB1
MFN2
MPZ
```

## Uso básico

```bash
# 1. Comprobar conectividad (opcional pero recomendable)
python test_conexion.py

# 2. Generar el BED
python genes_to_bed.py -i genes_ejemplo.txt -o genes_CMT.bed
```

Salida (`genes_CMT.bed`), formato BED de 6 columnas, ordenado por cromosoma/posición:

```
chr1    11976030    12013535    MFN2    .    -
chr1    161276020   161284206   MPZ     .    +
chrX    70443100    70452553    GJB1    .    -
chr17   15133295    15185710    PMP22   .    +
```

Si algún gen no se encuentra en Ensembl (símbolo mal escrito, sinónimo no
reconocido, etc.), se reporta en pantalla y se guarda en
`genes_CMT.no_encontrados.txt` para que se pueda revisar manualmente

## Opciones útiles

```bash
# Usar GRCh37/hg19 en vez de GRCh38 (usa el servidor espejo de Ensembl)
python genes_to_bed.py -i genes.txt -o genes.bed --assembly GRCh37

# No añadir el prefijo "chr" (algunos pipelines/referencias usan solo "17", "X", etc.)
python genes_to_bed.py -i genes.txt -o genes.bed --no-chr-prefix

# Especificar dónde guardar el log de genes no encontrados
python genes_to_bed.py -i genes.txt -o genes.bed --fallidos-log fallos.txt
```

## Integración en el pipeline completo de CMT

Este BED de genes es el **primer filtro** del flujo de priorización de
variantes (filtro de genes → separación intrón/exón → filtro MAF → filtro
de sinónimas). Se usa así con `bedtools`:

```bash
bedtools intersect -a variantes_qc.vcf -b genes_CMT.bed -u > variantes_genes_CMT.vcf
```

Para la separación intrón/exón posterior (basada en el transcrito MANE),
se necesita un paso adicional distinto, no cubierto por este script, que
genera BEDs de intrones/exones a partir del GTF de MANE (ver conversación
previa / `GenomicFeatures::intronsByTranscript()` en R).

## Notas técnicas

- El endpoint usado es `/lookup/symbol/<especie>/<gen>`, que devuelve las
  coordenadas del **gen completo**, no de un transcrito o exón concreto.
- Ensembl es 1-based; el script convierte automáticamente a 0-based para
  cumplir el estándar BED.
- Se respeta el rate limit de la API de Ensembl (pausa entre peticiones +
  reintentos automáticos con backoff si se recibe un HTTP 429).
- Si tu red institucional bloquea `rest.ensembl.org` (cortafuegos de
  universidad/hospital), puedes ejecutar `test_conexion.py` para
  diagnosticarlo rápidamente, o alternativamente resolver los genes con un
  GTF local (Ensembl/GENCODE descargado) en vez de la API, si el firewall
  es persistente.
