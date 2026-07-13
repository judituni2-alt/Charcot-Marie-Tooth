#!/usr/bin/env python3
"""
test_conexion.py

Comprobación rápida de que hay acceso a la API REST de Ensembl antes de
lanzar genes_to_bed.py sobre la lista completa de genes.

Uso:
    python test_conexion.py
"""
import json
import sys
import urllib.request

URL_PRUEBA = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/PMP22?content-type=application/json"

try:
    with urllib.request.urlopen(URL_PRUEBA, timeout=10) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    print("OK: conexión con Ensembl REST API funcionando correctamente.")
    print(f"Ejemplo (PMP22): {data['seq_region_name']}:{data['start']}-{data['end']} ({data.get('assembly_name', '?')})")
    sys.exit(0)
except Exception as e:
    print(f"ERROR: no se pudo contactar con la API de Ensembl: {e}", file=sys.stderr)
    print("Comprueba tu conexión a internet o si tu red/firewall institucional bloquea rest.ensembl.org", file=sys.stderr)
    sys.exit(1)
