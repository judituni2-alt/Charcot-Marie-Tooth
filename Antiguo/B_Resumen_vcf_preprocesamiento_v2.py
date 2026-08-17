#!/usr/bin/env python3
"""
summarize_vcf.py

Genera un resumen del VCF crudo final del pipeline de preprocesamiento:
 - nº de variantes totales
 - nº de SNPs vs indels
 - distribución de calidad (QUAL)
 - variantes por cromosoma

Uso:
    python summarize_vcf.py --vcf ruta/al/archivo.vcf.gz --out ruta/al/variantes.csv
    python summarize_vcf.py --vcf ruta/al/archivo.vcf.gz --out ruta/al/variantes.csv --log-level DEBUG
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd
from cyvcf2 import VCF

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resumen de un VCF (nº variantes, tipo, QUAL, por cromosoma)."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        type=Path,
        help="Ruta al VCF de entrada (puede ser .vcf o .vcf.gz).",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Ruta del CSV de salida con el detalle de variantes.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Nivel de logging (por defecto INFO).",
    )
    return parser.parse_args()


def summarize(vcf_path: Path, out_csv: Path) -> None:
    if not vcf_path.exists():
        logger.error("El VCF de entrada no existe: %s", vcf_path)
        sys.exit(1)

    t_start = time.perf_counter()
    logger.info("Abriendo VCF: %s", vcf_path)

    t0 = time.perf_counter()
    vcf = VCF(str(vcf_path))
    logger.info("VCF abierto en %.2f s", time.perf_counter() - t0)

    rows = []
    t0 = time.perf_counter()
    n = 0
    for variant in vcf:
        tipo = "SNP" if variant.is_snp else ("INDEL" if variant.is_indel else "OTRO")
        rows.append(
            {
                "chrom": variant.CHROM,
                "pos": variant.POS,
                "tipo": tipo,
                "qual": variant.QUAL,
                "filter": variant.FILTER or "PASS",
            }
        )
        n += 1
        if n % 100_000 == 0:
            logger.debug("Procesadas %d variantes...", n)

    logger.info("Lectura de %d variantes completada en %.2f s", n, time.perf_counter() - t0)

    df = pd.DataFrame(rows)

    if df.empty:
        logger.warning("No se encontraron variantes en el VCF: %s", vcf_path)
        return

    t0 = time.perf_counter()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    logger.info("CSV guardado en %s (%.2f s)", out_csv, time.perf_counter() - t0)

    resumen_tipo = df["tipo"].value_counts().rename("n_variantes")
    resumen_chrom = df["chrom"].value_counts().rename("n_variantes")

    print("=== Resumen general ===")
    print(f"Total de variantes: {len(df)}")
    print("\nPor tipo:")
    print(resumen_tipo.to_string())
    print("\nPor cromosoma (top 10):")
    print(resumen_chrom.head(10).to_string())

    if df["qual"].notna().any():
        print(
            f"\nQUAL medio: {df['qual'].mean():.2f}  |  QUAL mediana: {df['qual'].median():.2f}"
        )
    else:
        print("\nQUAL no disponible para ninguna variante.")

    logger.info("Resumen completado en %.2f s (total)", time.perf_counter() - t_start)


if __name__ == "__main__":
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    summarize(args.vcf, args.out)
