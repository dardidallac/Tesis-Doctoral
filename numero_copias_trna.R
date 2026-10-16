# ============================================================
# ANÁLISIS INTEGRAL: ΔRSCU vs tRNA GCN (Leguminosas)
# ============================================================

if (!require("coRdon")) BiocManager::install("coRdon")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
library(coRdon)
library(ggplot2)
library(dplyr)

# 1. CARGA DE ARCHIVOS FASTA (Secuencias Concatenadas)
cat("Selecciona el FASTA de ALTA expresión:\n")
high_file <- file.choose()
cat("Selecciona el FASTA de BAJA expresión:\n")
low_file <- file.choose()

seqs_high <- readSet(file = high_file)
seqs_low  <- readSet(file = low_file)

# 2. CARGA DE TABLA DE tRNAs (GCN)
cat("Selecciona el archivo CSV de tRNAs (columnas: anti_codon, GCN):\n")
trna_path <- file.choose()
trna_input <- read.csv(trna_path, stringsAsFactors = FALSE)

# 3. REFERENCIA Y FUNCIONES DE CONVERSIÓN
aa_ref <- c(AAA="K", AAC="N", AAG="K", AAT="N", ACA="T", ACC="T", ACG="T", ACT="T",
            AGA="R", AGC="S", AGG="R", AGT="S", ATA="I", ATC="I", ATG="M", ATT="I",
            CAA="Q", CAC="H", CAG="Q", CAT="H", CCA="P", CCC="P", CCG="P", CCT="P",
            CGA="R", CGC="R", CGG="R", CGT="R", CTA="L", CTC="L", CTG="L", CTT="L",
            GAA="E", GAC="D", GAG="E", GAT="D", GCA="A", GCC="A", GCG="A", GCT="A",
            GGA="G", GGC="G", GGG="G", GGT="G", GTA="V", GTC="V", GTG="V", GTT="V",
            TAA="*", TAC="Y", TAG="*", TAT="Y", TCA="S", TCC="S", TCG="S", TCT="S",
            TGA="*", TGC="C", TGG="W", TGT="C", TTA="L", TTC="F", TTG="L", TTT="F")

# Función para convertir Anticodón (3'-5') a Codón (5'-3')
anti_to_codon <- function(anti) {
  # 1. Complemento (A->T, T->A, G->C, C->G)
  comp <- chartr("ATGC", "TACG", anti)
  # 2. Reverso (para obtener dirección 5'-3')
  codon <- paste(rev(strsplit(comp, "")[[1]]), collapse = "")
  return(codon)
}

# 4. CÁLCULO DE RSCU
calc_rscu_fast <- function(seqs) {
  cc <- codonTable(seqs)
  counts <- colSums(cc@counts)
  names(counts) <- gsub("U", "T", names(counts))
  
  rscu_vals <- setNames(numeric(64), names(aa_ref))
  unique_aas <- unique(aa_ref[!(aa_ref %in% c("*", "M", "W"))])
  
  for (aa in unique_aas) {
    cods <- names(aa_ref[aa_ref == aa])
    obs <- counts[cods]
    obs[is.na(obs)] <- 0
    total_aa <- sum(obs)
    if (total_aa > 0) {
      rscu_vals[cods] <- obs / (total_aa / length(cods))
    }
  }
  return(rscu_vals)
}

rscu_high <- calc_rscu_fast(seqs_high)
rscu_low  <- calc_rscu_fast(seqs_low)

# 5. ENSAMBLAJE DE TABLA FINAL
trna_input$Codon <- sapply(trna_input$anti_codon, anti_to_codon)

results <- data.frame(
  Codon = names(aa_ref),
  AA = aa_ref,
  Delta_RSCU = as.numeric(rscu_high - rscu_low),
  stringsAsFactors = FALSE
) %>%
  left_join(trna_input[, c("Codon", "GCN")], by = "Codon") %>%
  mutate(GCN = ifelse(is.na(GCN), 0, GCN)) %>%
  filter(!(AA %in% c("*", "M", "W")))

# 6. VISUALIZACIÓN
# Gráfico 1: Correlación entre abundancia de tRNA y preferencia de uso
p1 <- ggplot(results, aes(x = GCN, y = Delta_RSCU)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", color = "blue", fill = "lightblue", alpha = 0.2) +
  geom_point(aes(color = AA), size = 3, alpha = 0.8) +
  geom_text(aes(label = Codon), vjust = -1, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Co-adaptación Codón-tRNA",
       subtitle = "Valores positivos en Delta RSCU indican codones optimizados",
       x = "Número de Copias de tRNA (GCN)", y = "Δ RSCU (Alta - Baja)")

print(p1)

# Mostrar resumen de correlación
cor_test <- cor.test(results$Delta_RSCU, results$GCN)
cat("\n--- ESTADÍSTICAS DE CO-ADAPTACIÓN ---\n")
cat("Correlación de Pearson:", round(cor_test$estimate, 3), "\n")
cat("P-valor:", cor_test$p.value, "\n")
