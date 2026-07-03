#!/usr/bin/env Rscript
###############################################################################
# plot_qc_metrics.R
#
# Genera gráficos resumen a partir de las salidas del pipeline:
#   - flagstat de samtools (% lecturas mapeadas, duplicados, etc.)
#   - métricas de duplicados de Picard
#
# Uso:
#   Rscript plot_qc_metrics.R resultados/ muestra_id
#
# Paquetes necesarios: ggplot2, dplyr, stringr
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript plot_qc_metrics.R <outdir> <sample_id>")
}

outdir <- args[1]
sample_id <- args[2]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

## --- 1. Leer flagstat de samtools ------------------------------------------

flagstat_path <- file.path(outdir, "qc_align", paste0(sample_id, "_flagstat.txt"))
flagstat_lines <- readLines(flagstat_path)

extraer_pct <- function(patron, lineas) {
  linea <- lineas[str_detect(lineas, patron)][1]
  pct <- str_extract(linea, "(?<=\\()[0-9.]+(?=%)")
  as.numeric(pct)
}

metricas_align <- data.frame(
  metrica = c("Mapeadas", "Emparejadas correctamente", "Duplicadas"),
  porcentaje = c(
    extraer_pct("mapped \\(", flagstat_lines),
    extraer_pct("properly paired", flagstat_lines),
    extraer_pct("duplicates", flagstat_lines)
  )
)

p1 <- ggplot(metricas_align, aes(x = metrica, y = porcentaje, fill = metrica)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
  labs(
    title = paste("Métricas de alineamiento -", sample_id),
    subtitle = "Línea roja discontinua: umbral de referencia (90%)",
    x = NULL, y = "Porcentaje (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "qc_align", paste0(sample_id, "_align_metrics.png")),
       p1, width = 6, height = 4, dpi = 150)

## --- 2. Leer métricas de duplicados de Picard -------------------------------

dup_path <- file.path(outdir, "dedup", paste0(sample_id, "_dup_metrics.txt"))
dup_raw <- readLines(dup_path)

# El archivo de Picard tiene una sección de métricas con cabecera "LIBRARY"
header_idx <- which(str_detect(dup_raw, "^LIBRARY"))
dup_table <- read.table(text = dup_raw[header_idx:(header_idx + 1)],
                         header = TRUE, sep = "\t", fill = TRUE)

pct_dup <- as.numeric(dup_table$PERCENT_DUPLICATION[1]) * 100

p2 <- ggplot(data.frame(x = "Duplicados PCR", pct = pct_dup),
              aes(x = x, y = pct)) +
  geom_col(fill = "#c77c7c", width = 0.4) +
  ylim(0, 100) +
  labs(title = paste("Porcentaje de duplicados PCR -", sample_id),
       x = NULL, y = "Porcentaje (%)") +
  theme_minimal(base_size = 13)

ggsave(file.path(outdir, "dedup", paste0(sample_id, "_duplication.png")),
       p2, width = 4, height = 4, dpi = 150)

cat("Gráficos guardados en:\n")
cat(" -", file.path(outdir, "qc_align", paste0(sample_id, "_align_metrics.png")), "\n")
cat(" -", file.path(outdir, "dedup", paste0(sample_id, "_duplication.png")), "\n")
