#!/usr/bin/env python3
"""
summarize_vcf.py

Genera un resumen del VCF crudo final del pipeline de preprocesamiento:
 - nº de variantes totales
 - nº de SNPs vs indels
 - distribución de calidad (QUAL)
 - variantes por cromosoma

"""

import sys
import pandas as pd

def summarize(vcf_path: str, out_csv: str) -> None:
    vcf = VCF(vcf_path)

    rows = []
    for variant in vcf:
        tipo = "SNP" if variant.is_snp else ("INDEL" if variant.is_indel else "OTRO")
        rows.append({
            "chrom": variant.CHROM,
            "pos": variant.POS,
            "tipo": tipo,
            "qual": variant.QUAL,
            "filter": variant.FILTER or "PASS",
        })

    df = pd.DataFrame(rows)

    if df.empty:
        print("No se encontraron variantes en el VCF.")
        return

    resumen_tipo = df["tipo"].value_counts().rename("n_variantes")
    resumen_chrom = df["chrom"].value_counts().rename("n_variantes")

    print("=== Resumen general ===")
    print(f"Total de variantes: {len(df)}")
    print("\nPor tipo:")
    print(resumen_tipo.to_string())
    print("\nPor cromosoma (top 10):")
    print(resumen_chrom.head(10).to_string())
    print(f"\nQUAL medio: {df['qual'].mean():.2f}  |  QUAL mediana: {df['qual'].median():.2f}")


if __name__ == "__main__":
    summarize(
        "/home/judit/TFM/Charcot-Marie-Tooth/resultados/variants/TEST2.g.vcf.gz",
        "/home/judit/TFM/Charcot-Marie-Tooth/resultados/variants/variantes.csv"
    )
