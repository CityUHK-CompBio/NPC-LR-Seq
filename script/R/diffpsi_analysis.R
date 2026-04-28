# Analyze SUPPA2 differential splicing results, draw a volcano plot, and save processed output.

library(dplyr)
library(stringr)
library(svglite)
library(EnhancedVolcano)

suppa_diff_file <- "/path/to/output/suppa/project_diffSplice.dpsi"
figure_file <- "/path/to/output/figures/suppa_diffsplice.svg"
rdata_file <- "/path/to/output/rdata/suppa_diffsplice.Rdata"

suppa_diffsplice <- read.table(suppa_diff_file, header = TRUE)
suppa_diffsplice <- suppa_diffsplice %>% filter(!is.nan(deltaPSI))
names(suppa_diffsplice) <- c("delta_PSI", "p_value")
suppa_diffsplice$ASE <- stringr::word(rownames(suppa_diffsplice), 2, 2, sep = ";|:")

table((suppa_diffsplice %>% filter(abs(delta_PSI) > 0.1 & p_value < 0.05))$ASE)

svglite(figure_file, width = 10, height = 8)
p <- EnhancedVolcano(
  suppa_diffsplice,
  lab = rownames(suppa_diffsplice),
  x = "delta_PSI",
  y = "p_value",
  xlim = c(-1, 1),
  ylim = c(0, 4),
  xlab = "ΔPSI",
  boxedLabels = FALSE,
  axisLabSize = 16,
  ylab = bquote(~-Log[10] ~ italic(P)),
  subtitle = "",
  title = "",
  pCutoff = 0.05,
  col = c("grey30", "grey30", "grey30", "red2"),
  pointSize = 3.5,
  colAlpha = 0.4,
  FCcutoff = 0.1,
  labSize = 3.5
)
plot(p)
dev.off()

save(suppa_diffsplice, file = rdata_file)
