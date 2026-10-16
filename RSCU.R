# ============================================================
# ANÁLISIS DE USO DE CODONES EN CDS DE UN GENOMA (SIN GRUPOS)
# ============================================================

# ===============================
# CONFIGURACIÓN E INSTALACIÓN
# ===============================

install.packages("devtools", dependencies = TRUE)
library(devtools)

if (!require("coRdon")) BiocManager::install("coRdon")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mariodosreis/tai")

load.lib <- c('ggplot2', 'ggrepel', 'ggthemes', 
              'ggpubr', 'FactoMineR', 'factoextra', 
              'BiocManager', 'devtools', 'data.table')

install.lib <- load.lib[!load.lib %in% installed.packages()]
if (length(install.lib) > 0) {
  cat("Installing packages:", paste(install.lib, collapse = ", "), "\n")
  install.packages(install.lib, dependencies = TRUE)
}

suppressMessages(sapply(load.lib, require, character.only = TRUE))

library('coRdon')

# ===============================
# CARGA DE ARCHIVO
# ===============================

cat("Selecciona el archivo con las secuencias codificantes (CDS):\n")
cds <- file.choose()
Filter_smaller_than_80_codons <- TRUE

# ===============================
# LECTURA DE SECUENCIAS
# ===============================

coding_Sequences <- readSet(file = cds)
if (Filter_smaller_than_80_codons){
  coding_Sequences <- coding_Sequences[which(coding_Sequences@ranges@width >= 80*3)]
}

cc_cds <- codonTable(coding_Sequences)

# ===============================
# FUNCIONES GC3 y GC12
# ===============================

calcGC3 <- function(countsTable){
  counts <- countsTable@counts
  len <- countsTable@len
  GC3 <- counts[,grep("..G|..C",colnames(counts))]
  GC3.sum <- apply(as.matrix(GC3),1,sum)
  GC3.res <- GC3.sum/len*100
  return(GC3.res)
}

calcGC12 <- function(countsTable){
  counts <- countsTable@counts
  len <- countsTable@len
  GC1 <- counts[,grep("G..|C..",colnames(counts))]
  GC1.sum <- apply(as.matrix(GC1),1,sum)
  GC2 <- counts[,grep(".G.|.C.",colnames(counts))]
  GC2.sum <- apply(as.matrix(GC2),1,sum)
  GC12.res <- ((GC1.sum/len*100)+(GC2.sum/len*100))/2
  return(GC12.res)
}

# ===============================
# CÁLCULO DE ÍNDICES
# ===============================

milc <- MILC(cc_cds)[,1]
enc <- ENC(cc_cds)
mcb <- MCB(cc_cds)[,1]
encp <- ENCprime(cc_cds, subsets = list(cc_cds))[,2]
scuo <- SCUO(cc_cds)
gc12 <- calcGC12(cc_cds)
gc3 <- calcGC3(cc_cds)

table <- data.frame(cds = coding_Sequences@ranges@NAMES, 
                    MILC = milc, Enc = enc, MCB = mcb, 
                    EncPrime = encp, SCUO = scuo, 
                    GC3 = gc3, GC12 = gc12)

# ===============================
# PLOTS BÁSICOS (Boxplots)
# ===============================

library('data.table')
ltable <- melt(setDT(table), id.vars = 'cds', variable.name = "Index")

plot_box <- function(df, index_name) {
  subdf <- df[df$Index == index_name,]
  ggplot(subdf, aes(x="CDS", y=value)) +
    geom_boxplot(fill="#00BFC4", alpha=0.7) +
    theme_minimal() + 
    ylab(index_name) +
    ggtitle(paste(index_name, "en CDS"))
}

print(plot_box(ltable, "MILC"))
print(plot_box(ltable, "Enc"))
print(plot_box(ltable, "MCB"))
print(plot_box(ltable, "EncPrime"))
print(plot_box(ltable, "SCUO"))
print(plot_box(ltable, "GC3"))

# ===============================
# PLOT: Nc vs GC3 (Nc Plot)
# ===============================

fgc3s <- function (gc3) {
  a = -6; b = 34; c = 1.025
  x = a + gc3 + (b/(gc3^2 + (c - gc3)^2))
  return(x)
}

ncplot <- data.frame(x=seq(0,1,0.01)*100, y=fgc3s(seq(0,1,0.01)))

pNCPlot <- ggplot(data=table, aes(x=GC3, y=Enc)) +
  geom_point(color="#00BFC4", size=2) +
  geom_line(data=ncplot, aes(x=x,y=y), color='black', linetype="dashed")+
  theme_classic() +
  labs(x="GC3s (%)", y="Nc", title="ENc vs GC3s")

print(pNCPlot)

# ===============================
# NEUTRALITY PLOT
# ===============================

# ===============================
# NEUTRALITY PLOT - ANÁLISIS FORMAL
# ===============================

# Modelo lineal
model <- lm(GC12 ~ GC3, data = table)
model_summary <- summary(model)

# Parámetros
slope <- coef(model)[2]
intercept <- coef(model)[1]
r_value <- cor(table$GC3, table$GC12, method = "pearson")
r_squared <- model_summary$r.squared
p_value <- model_summary$coefficients[2,4]

# Error estándar de la pendiente
se_slope <- model_summary$coefficients[2,2]

# Test: pendiente distinta de 1
t_value_vs1 <- (slope - 1) / se_slope
df <- nrow(table) - 2
p_value_vs1 <- 2 * pt(-abs(t_value_vs1), df)

# Contribución mutacional (modelo clásico)
mutational_contribution <- slope * 100
selective_contribution <- (1 - slope) * 100

# Mostrar resultados
cat("\n========== NEUTRALITY PLOT RESULTS ==========\n")
cat("Intercept:", round(intercept,4), "\n")
cat("Slope:", round(slope,4), "\n")
cat("Pearson r:", round(r_value,4), "\n")
cat("R²:", round(r_squared,4), "\n")
cat("p-value (slope ≠ 0):", signif(p_value,4), "\n")
cat("p-value (slope ≠ 1):", signif(p_value_vs1,4), "\n")
cat("Mutational contribution (%):", round(mutational_contribution,2), "\n")
cat("Selective constraint (%):", round(selective_contribution,2), "\n")
cat("=============================================\n")

# ===============================
# COA - RAW CODON COUNTS
# ===============================

cc4COA <- cc_cds@counts
cc4COA = cc4COA[,-c(15,49,51,57,59)]  # eliminar stop codons
res.ca <- CA(cc4COA, graph = FALSE)
varpc <- get_eigenvalue(res.ca)[1:2,2]

cds_pcs <- data.frame(get_ca_row(res.ca)$coord[,1:2])
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

coaplot_raw <- ggplot(cds_pcs, aes(x=Dim.1, y=Dim.2)) +
  geom_point(size=2, color="#00BFC4") +
  theme_minimal() +
  labs(title="Raw codon counts", x=xlabtext, y=ylabtext)

print(coaplot_raw)

# ===============================
# COA - RSCU
# ===============================

counts_cds <- cc_cds@counts
aminoacids <- c('K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P',
                'R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
                '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F')
codons <- colnames(counts_cds)
dnacodons <- data.frame(codon = colnames(counts_cds), aa = aminoacids)

getRSCU <- function(codcounts) {
  rscu <- data.frame(matrix(ncol = 64, nrow = dim(codcounts)[1]))
  for (i in 1:64) {
    aa <- dnacodons[dnacodons[,1] == codons[i],2]
    cods <- which(dnacodons[,2] == aa)
    if (length(cods) == 1) {
      rscu[,i] <- 1
    } else {
      rscu[,i] <- codcounts[,i]/rowSums(codcounts[,cods])*length(cods)
    }
  }
  colnames(rscu) = colnames(codcounts)
  rscu[is.na(rscu)] <- 1
  return(rscu)
}

rscu <- getRSCU(counts_cds)
rscu = rscu[,-c(15,49,51,57,59)]
avg_rscu <- colMeans(rscu)
rscu4COA <- rbind.data.frame(rscu,avg_rscu)
ncds <- dim(rscu)[1]
totalrows <- dim(rscu4COA)[1]

res.ca <- CA(rscu4COA, graph = FALSE, row.sup = (ncds+1):totalrows)
varpc <- get_eigenvalue(res.ca)[1:2,2]

cds_pcs <- data.frame(get_ca_row(res.ca)$coord[,1:2])
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

coaplot_rscu <- ggplot(cds_pcs, aes(x=Dim.1, y=Dim.2)) +
  geom_point(size=2, color="#00BFC4") +
  theme_minimal() +
  labs(title="RSCU", x=xlabtext, y=ylabtext)

print(coaplot_rscu)

# ===============================
# COA - RAW CODON COUNTS (con etiquetas de codones)
# ===============================

cc4COA <- cc_cds@counts
cc4COA = cc4COA[,-c(15,49,51,57,59)]  # eliminar stop codons

# Análisis de correspondencia
res.ca.raw <- CA(cc4COA, graph = FALSE)
varpc.raw <- get_eigenvalue(res.ca.raw)[1:2,2]

# Coordenadas de los CDS
cds_coords.raw <- data.frame(get_ca_row(res.ca.raw)$coord[,1:2])
colnames(cds_coords.raw) <- c("Dim.1", "Dim.2")

# Coordenadas de los codones
codon_coords.raw <- data.frame(get_ca_col(res.ca.raw)$coord[,1:2])
colnames(codon_coords.raw) <- c("Dim.1", "Dim.2")
codon_coords.raw$codon <- rownames(get_ca_col(res.ca.raw)$coord)

# Etiquetas de ejes
xlabtext.raw <- paste("C1 (", round(varpc.raw[1], 2), " %)", sep = "")
ylabtext.raw <- paste("C2 (", round(varpc.raw[2], 2), " %)", sep = "")

# Plot extendido con etiquetas repelidas
pCOA_raw_ext <- ggplot() +
  geom_point(data = cds_coords.raw, aes(x = Dim.1, y = Dim.2), color = "#00BFC4", size = 2) +
  geom_point(data = codon_coords.raw, aes(x = Dim.1, y = Dim.2), color = "black", size = 1.5) +
  geom_text_repel(
    data = codon_coords.raw,
    aes(x = Dim.1, y = Dim.2, label = codon),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    segment.color = "gray60"
  ) +
  theme_minimal() +
  labs(title = "Raw codon counts - codons",
       x = xlabtext.raw, y = ylabtext.raw)

print(pCOA_raw_ext)


# ===============================
# COA - RSCU (con etiquetas de codones)
# ===============================

# Generar RSCU si no está
counts_cds <- cc_cds@counts
aminoacids <- c('K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P',
                'R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
                '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F')
codons <- colnames(counts_cds)
dnacodons <- data.frame(codon = colnames(counts_cds), aa = aminoacids)

getRSCU <- function(codcounts) {
  rscu <- data.frame(matrix(ncol = 64, nrow = dim(codcounts)[1]))
  for (i in 1:64) {
    aa <- dnacodons[dnacodons[,1] == codons[i],2]
    cods <- which(dnacodons[,2] == aa)
    if (length(cods) == 1) {
      rscu[,i] <- 1
    } else {
      rscu[,i] <- codcounts[,i]/rowSums(codcounts[,cods])*length(cods)
    }
  }
  colnames(rscu) = colnames(codcounts)
  rscu[is.na(rscu)] <- 1
  return(rscu)
}

rscu <- getRSCU(counts_cds)
rscu = rscu[,-c(15,49,51,57,59)]
avg_rscu <- colMeans(rscu)
rscu4COA <- rbind.data.frame(rscu, avg_rscu)
ncds <- dim(rscu)[1]
totalrows <- dim(rscu4COA)[1]

# Análisis de correspondencia RSCU
res.ca.rscu <- CA(rscu4COA, graph = FALSE, row.sup = (ncds+1):totalrows)
varpc.rscu <- get_eigenvalue(res.ca.rscu)[1:2,2]

# Coordenadas de CDS
cds_coords.rscu <- data.frame(get_ca_row(res.ca.rscu)$coord[,1:2])
colnames(cds_coords.rscu) <- c("Dim.1", "Dim.2")

# Coordenadas de codones
codon_coords.rscu <- data.frame(get_ca_col(res.ca.rscu)$coord[,1:2])
colnames(codon_coords.rscu) <- c("Dim.1", "Dim.2")
codon_coords.rscu$codon <- rownames(get_ca_col(res.ca.rscu)$coord)

# Etiquetas de ejes
xlabtext.rscu <- paste("C1 (", round(varpc.rscu[1], 2), " %)", sep = "")
ylabtext.rscu <- paste("C2 (", round(varpc.rscu[2], 2), " %)", sep = "")

# Plot extendido con etiquetas repelidas
pCOA_rscu_ext <- ggplot() +
  geom_point(data = cds_coords.rscu, aes(x = Dim.1, y = Dim.2), color = "#00BFC4", size = 2) +
  geom_point(data = codon_coords.rscu, aes(x = Dim.1, y = Dim.2), color = "black", size = 1.5) +
  geom_text_repel(
    data = codon_coords.rscu,
    aes(x = Dim.1, y = Dim.2, label = codon),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    segment.color = "gray60"
  ) +
  theme_minimal() +
  labs(title = "RSCU - codons",
       x = xlabtext.rscu, y = ylabtext.rscu)

print(pCOA_rscu_ext)


# ===============================
# SPEARMAN: C1 vs GC3 (una planta)
# ===============================

# Extraer variables
C1 <- cds_coords.rscu$Dim.1
GC3 <- table$GC3

# Test de correlación Spearman
spearman_test <- cor.test(C1, GC3,
                          method = "spearman",
                          use = "complete.obs",
                          exact = FALSE)

# Extraer resultados
rho <- as.numeric(spearman_test$estimate)
p_value <- spearman_test$p.value

# Mostrar resultados
cat("\n===== SPEARMAN C1 vs GC3 =====\n")
cat("Spearman rho:", round(rho, 4), "\n")
cat("p-value:", signif(p_value, 4), "\n")
cat("===============================\n")


# ===============================
# RSCU PROMEDIO - BARPLOT GLOBAL
# ===============================

# Crear dataframe ordenado
rscu_avg_df <- data.frame(
  Codon = names(avg_rscu),
  RSCU = avg_rscu
)

# Ordenar por aminoácido (opcional pero recomendable)
rscu_avg_df$AA <- aminoacids[match(rscu_avg_df$Codon, codons)]
rscu_avg_df <- rscu_avg_df[order(rscu_avg_df$AA),]

# Definir codones preferidos (RSCU > 1)
rscu_avg_df$Preferred <- ifelse(rscu_avg_df$RSCU > 1, "Preferred", "Non-preferred")

# Gráfico
pRSCU_bar <- ggplot(rscu_avg_df, aes(x = Codon, y = RSCU, fill = Preferred)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("Preferred" = "#00BFC4", 
                               "Non-preferred" = "grey70")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, face = "bold"),
        panel.grid.major.x = element_blank()) +
  labs(title = "Promedio RSCU (Todos los genes)",
       y = "RSCU",
       x = "Codon")

print(pRSCU_bar)



