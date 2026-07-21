#!/usr/bin/env python3
"""
genes_to_bed.py

Genera un archivo BED de genes a partir de una lista de nombres de genes
(uno por línea) en un archivo de texto plano, usando la API REST de
Ensembl (https://rest.ensembl.org). Solo GRCh38.

Uso:
    python genes_to_bed.py -i genes.txt -o genes_CMT.bed

Formato de entrada (genes.txt), un gen por línea:
    PMP22
    GJB1
    MFN2
    MPZ

Salida (BED, 0-based, estándar bedtools):
    chr17   15133295   15185710   PMP22   .   +
    chrX    70443100   70452553   GJB1    .   -
    ...
"""

import argparse
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

ENSEMBL_SERVER = "https://rest.ensembl.org"
REQUEST_DELAY_SECONDS = 0.1  # margen de seguridad bajo el rate limit de Ensembl


def leer_genes(path_txt: Path) -> list[str]:
    """Lee la lista de genes desde un archivo de texto plano, un gen por línea."""
    genes = []
    with open(path_txt, "r") as f:
        for linea in f:
            gen = linea.strip()
            if gen and not gen.startswith("#"):
                genes.append(gen)
    if not genes:
        sys.exit(f"ERROR: no se encontraron genes válidos en {path_txt}")
    return genes


def consultar_ensembl_gen(simbolo: str, species: str) -> dict | None:
    """Consulta el endpoint lookup/symbol de Ensembl (GRCh38) para un gen concreto."""
    url = f"{ENSEMBL_SERVER}/lookup/symbol/{species}/{simbolo}?content-type=application/json"
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code != 404:
            print(f"  [ERROR] fallo HTTP consultando {simbolo}: {e}", file=sys.stderr)
        return None
    except (urllib.error.URLError, TimeoutError) as e:
        print(f"  [ERROR] fallo de red consultando {simbolo}: {e}", file=sys.stderr)
        return None


def genes_a_bed(genes: list[str], species: str, chr_prefix: bool) -> tuple[list[str], list[str]]:
    """Convierte la lista de genes en líneas BED. Devuelve (lineas_bed, no_encontrados)."""
    lineas_bed = []
    no_encontrados = []

    for i, gen in enumerate(genes, start=1):
        print(f"[{i}/{len(genes)}] Consultando {gen} ...", file=sys.stderr)
        data = consultar_ensembl_gen(gen, species)

        if data is None:
            print(f"  [ERROR] {gen} no encontrado en Ensembl ({species}, GRCh38)", file=sys.stderr)
            no_encontrados.append(gen)
            time.sleep(REQUEST_DELAY_SECONDS)
            continue

        try:
            chrom = data["seq_region_name"]
            start = int(data["start"]) - 1  # Ensembl es 1-based; BED es 0-based
            end = int(data["end"])
            nombre = data.get("display_name", gen)
            strand = "+" if data.get("strand", 1) == 1 else "-"
        except KeyError as e:
            print(f"  [ERROR] respuesta incompleta de Ensembl para {gen} (falta {e})", file=sys.stderr)
            no_encontrados.append(gen)
            time.sleep(REQUEST_DELAY_SECONDS)
            continue

        if chr_prefix and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        lineas_bed.append(f"{chrom}\t{start}\t{end}\t{nombre}\t.\t{strand}")
        time.sleep(REQUEST_DELAY_SECONDS)

    return lineas_bed, no_encontrados


def escribir_bed(lineas_bed: list[str], path_out: Path):
    """Escribe las líneas BED ordenadas por cromosoma/posición."""
    def clave_orden(linea):
        campos = linea.split("\t")
        chrom = campos[0].replace("chr", "")
        if chrom.isdigit():
            return (0, int(chrom), int(campos[1]))
        orden_especial = {"X": 100, "Y": 101, "MT": 102, "M": 102}
        return (1, orden_especial.get(chrom, 999), int(campos[1]))

    lineas_ordenadas = sorted(lineas_bed, key=clave_orden)
    with open(path_out, "w") as f:
        for linea in lineas_ordenadas:
            f.write(linea + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Genera un archivo BED (GRCh38) a partir de una lista de símbolos génicos, usando la API REST de Ensembl."
    )
    parser.add_argument("-i", "--input", required=True, type=Path,
                         help="Archivo de texto plano con un nombre de gen por línea")
    parser.add_argument("-o", "--output", required=True, type=Path,
                         help="Archivo BED de salida")
    parser.add_argument("--species", default="homo_sapiens",
                         help="Especie para la consulta a Ensembl (default: homo_sapiens)")
    parser.add_argument("--no-chr-prefix", action="store_true",
                         help="No añadir el prefijo 'chr' a los cromosomas")

    args = parser.parse_args()

    genes = leer_genes(args.input)
    print(f"Se han leído {len(genes)} genes desde {args.input}", file=sys.stderr)

    lineas_bed, no_encontrados = genes_a_bed(genes, species=args.species, chr_prefix=not args.no_chr_prefix)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    escribir_bed(lineas_bed, args.output)
    print(f"\nBED generado: {args.output} ({len(lineas_bed)} genes)", file=sys.stderr)

    if no_encontrados:
        log_path = args.output.with_suffix(".no_encontrados.txt")
        with open(log_path, "w") as f:
            f.write("\n".join(no_encontrados) + "\n")
        print(f"AVISO: {len(no_encontrados)} genes no encontrados, ver {log_path}", file=sys.stderr)
        print("Genes no encontrados: " + ", ".join(no_encontrados), file=sys.stderr)


if __name__ == "__main__":
    main()
