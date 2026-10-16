# ============================================================
# CÁLCULO DE ΔRSCU (VERSIÓN ROBUSTA)
# ============================================================

if (!require("coRdon")) BiocManager::install("coRdon")
if (!require("ggplot2")) install.packages("ggplot2")
library(coRdon)
library(ggplot2)

# 1. CARGA DE ARCHIVOS
cat("1. Selecciona el archivo FASTA de ALTA expresión:\n")
high_file <- file.choose()
cat("2. Selecciona el archivo FASTA de BAJA expresión:\n")
low_file <- file.choose()

seqs_high <- readSet(file = high_file)
seqs_low  <- readSet(file = low_file)

cc_high <- codonTable(seqs_high)
cc_low  <- codonTable(seqs_low)

# 2. REFERENCIA AA-CODÓN (Normalizada a T para evitar conflictos con U)
aa_ref <- c(AAA="K", AAC="N", AAG="K", AAT="N", ACA="T", ACC="T", ACG="T", ACT="T",
            AGA="R", AGC="S", AGG="R", AGT="S", ATA="I", ATC="I", ATG="M", ATT="I",
            CAA="Q", CAC="H", CAG="Q", CAT="H", CCA="P", CCC="P", CCG="P", CCT="P",
            CGA="R", CGC="R", CGG="R", CGT="R", CTA="L", CTC="L", CTG="L", CTT="L",
            GAA="E", GAC="D", GAG="E", GAT="D", GCA="A", GCC="A", GCG="A", GCT="A",
            GGA="G", GGC="G", GGG="G", GGT="G", GTA="V", GTC="V", GTG="V", GTT="V",
            TAA="*", TAC="Y", TAG="*", TAT="Y", TCA="S", TCC="S", TCG="S", TCT="S",
            TGA="*", TGC="C", TGG="W", TGT="C", TTA="L", TTC="F", TTG="L", TTT="F")

# 3. FUNCIÓN MEJORADA PARA CALCULAR RSCU
calc_group_rscu <- function(cc_obj) {
  counts <- colSums(cc_obj@counts)
  # Normalizar nombres de columnas: cambiar U por T si es necesario
  names(counts) <- gsub("U", "T", names(counts))
  
  rscu_vals <- setNames(numeric(64), names(aa_ref))
  unique_aas <- unique(aa_ref[!(aa_ref %in% c("*", "M", "W"))])
  
  for (aa in unique_aas) {
    cods <- names(aa_ref[aa_ref == aa])
    obs <- counts[cods]
    # Reemplazar NAs por 0
    obs[is.na(obs)] <- 0
    
    total_aa <- sum(obs)
    if (total_aa > 0) {
      expected <- total_aa / length(cods)
      rscu_vals[cods] <- obs / expected
    } else {
      rscu_vals[cods] <- 0
    }
  }
  return(rscu_vals)
}

# 4. CÁLCULO
rscu_high <- calc_group_rscu(cc_high)
rscu_low  <- calc_group_rscu(cc_low)
delta_rscu <- rscu_high - rscu_low

# 5. RESULTADOS
results <- data.frame(
  Codon = names(delta_rscu),
  AA = aa_ref[names(delta_rscu)],
  RSCU_High = as.numeric(rscu_high),
  RSCU_Low = as.numeric(rscu_low),
  Delta_RSCU = as.numeric(delta_rscu),
  stringsAsFactors = FALSE
)

# Filtrar aminoácidos no informativos
results <- results[!(results$AA %in% c("*", "M", "W")), ]

# 6. PLOT
p <- ggplot(results, aes(x = reorder(Codon, AA), y = Delta_RSCU, fill = Delta_RSCU > 0)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "#F8766D"), 
                    labels = c("TRUE" = "Favorecido en Alta", "FALSE" = "Favorecido en Baja"),
                    name = "Tendencia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, family = "mono")) +
  labs(title = "Análisis de Delta RSCU",
       subtitle = "Positivo = Codón Óptimo (Potencial)",
       x = "Codón", y = "Delta RSCU")

print(p)

# Mostrar tabla simplificada
cat("\n--- TOP 10 CODONES ÓPTIMOS (Mayor Delta) ---\n")
print(head(results[order(results$Delta_RSCU, decreasing = TRUE), c("Codon", "AA", "Delta_RSCU")], 10))