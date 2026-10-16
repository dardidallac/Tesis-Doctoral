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

# ===============================
# CARGA DE ARCHIVOS
# ===============================

cat("Selecciona el archivo con las secuencias codificantes (CDS):\n")
cds <- file.choose()
Filter_smaller_than_80_codons <- TRUE

cat("Selecciona el archivo de referencia (CDS):\n")
ref_cds <- file.choose()

library('coRdon')

# ===============================
# LECTURA DE SECUENCIAS
# ===============================

coding_Sequences <- readSet(file = cds)
if (Filter_smaller_than_80_codons){
  coding_Sequences <- coding_Sequences[which(coding_Sequences@ranges@width >= 80*3)]
}

cc_cds <- codonTable(coding_Sequences)

ref_coding_Sequences <- readSet(file = ref_cds)
if (Filter_smaller_than_80_codons){
  ref_coding_Sequences <- ref_coding_Sequences[which(ref_coding_Sequences@ranges@width >= 80*3)]
}
ref_cc_cds <- codonTable(ref_coding_Sequences)

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

table <- data.frame(cds = coding_Sequences@ranges@NAMES, MILC = milc, Enc = enc, 
                    MCB = mcb, EncPrime = encp, SCUO = scuo, GC3 = gc3, GC12 = gc12)

# Extraer si es conserved o nonconserved del nombre
table$Group <- ifelse(grepl("nonconserved", table$cds, ignore.case = TRUE), 
                      "nonconserved", "conserved")

# ===============================
# PLOTS BÁSICOS (Boxplots)
# ===============================

library('data.table')
ltable <- melt(setDT(table), id.vars = c('cds', 'Group'), variable.name = "Index")

# función para boxplots rápidos
plot_box <- function(df, index_name) {
  subdf <- df[df$Index == index_name,]
  ggplot(subdf, aes(x=Group, y=value, fill=Group)) +
    geom_boxplot(alpha=0.7) +
    scale_fill_manual(values=c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
    theme_minimal() + 
    ylab(index_name) +
    ggtitle(paste(index_name, "por grupo"))
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

pNCPlot <- ggplot(data=table, aes(x=GC3, y=Enc, color=Group)) +
  geom_point(size=2) +
  geom_line(data=ncplot, aes(x=x,y=y), color='black', linetype="dashed")+
  scale_color_manual(values=c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
  theme_classic() +
  labs(x="GC3s (%)", y="Nc", title="Nc vs GC3s")

print(pNCPlot)

# ===============================
# NEUTRALITY PLOT
# ===============================

neutrality_plot <- ggscatter(table, x = "GC3", y = "GC12", color = "Group",
                             palette = c("#00BFC4", "#F8766D"),
                             add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = Group), label.x = 30) +
  theme_minimal() +
  labs(x="GC3 (%)", y="GC12 (%)", title="Neutrality plot")

print(neutrality_plot)

# ===============================
# COA - RAW CODON COUNTS (colores por grupo)
# ===============================

cc4COA <- cc_cds@counts
cc4COA = cc4COA[,-c(15,49,51,57,59)]  # eliminar stop codons
res.ca <- CA(cc4COA, graph = FALSE)
varpc <- get_eigenvalue(res.ca)[1:2,2]

cds_pcs <- data.frame(get_ca_row(res.ca)$coord[,1:2])
cds_pcs$Group <- table$Group
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

coaplot_raw <- ggplot(cds_pcs, aes(x=Dim.1, y=Dim.2, color=Group)) +
  geom_point(size=2) +
  scale_color_manual(values=c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
  theme_minimal() +
  labs(title="Raw codon counts", x=xlabtext, y=ylabtext)

print(coaplot_raw)

# ===============================
# COA - RSCU (colores por grupo)
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
cds_pcs$Group <- table$Group

xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

coaplot_rscu <- ggplot(cds_pcs, aes(x=Dim.1, y=Dim.2, color=Group)) +
  geom_point(size=2) +
  scale_color_manual(values=c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
  theme_minimal() +
  labs(title="RSCU", x=xlabtext, y=ylabtext)

print(coaplot_rscu)













library(ggrepel)

# ===============================
# COA - RAW CODON COUNTS (etiquetas repelidas)
# ===============================

# Análisis de correspondencia
res.ca <- CA(cc4COA, graph = FALSE)
varpc <- get_eigenvalue(res.ca)[1:2,2]

# Coordenadas de los CDS
cds_coords <- data.frame(get_ca_row(res.ca)$coord[,1:2])
colnames(cds_coords) <- c("Dim.1", "Dim.2")
cds_coords$Group <- table$Group

# Coordenadas de los codones
codon_coords <- data.frame(get_ca_col(res.ca)$coord[,1:2])
colnames(codon_coords) <- c("Dim.1", "Dim.2")
codon_coords$codon <- rownames(get_ca_col(res.ca)$coord)

# Etiquetas de ejes
xlabtext <- paste("C1 (", round(varpc[1], 2), " %)", sep = "")
ylabtext <- paste("C2 (", round(varpc[2], 2), " %)", sep = "")

# Plot extendido con etiquetas repelidas
pCOA_raw_ext <- ggplot() +
  geom_point(data = cds_coords, aes(x = Dim.1, y = Dim.2, color = Group), size = 2) +
  geom_point(data = codon_coords, aes(x = Dim.1, y = Dim.2), color = "black", size = 1.5) +
  geom_text_repel(
    data = codon_coords,
    aes(x = Dim.1, y = Dim.2, label = codon),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    segment.color = "gray60",
    min.segment.length = 0.1
  ) +
  scale_color_manual(values = c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
  theme_minimal() +
  labs(title = "Raw codon counts - codons",
       x = xlabtext, y = ylabtext)

print(pCOA_raw_ext)

# ===============================
# COA - RSCU (etiquetas repelidas)
# ===============================

# Análisis de correspondencia
res.ca <- CA(rscu4COA, graph = FALSE, row.sup = (ncds+1):totalrows)
varpc <- get_eigenvalue(res.ca)[1:2,2]

# Coordenadas de CDS
cds_coords <- data.frame(get_ca_row(res.ca)$coord[,1:2])
colnames(cds_coords) <- c("Dim.1", "Dim.2")
cds_coords$Group <- table$Group

# Coordenadas de codones
codon_coords <- data.frame(get_ca_col(res.ca)$coord[,1:2])
colnames(codon_coords) <- c("Dim.1", "Dim.2")
codon_coords$codon <- rownames(get_ca_col(res.ca)$coord)

# Etiquetas de ejes
xlabtext <- paste("C1 (", round(varpc[1], 2), " %)", sep = "")
ylabtext <- paste("C2 (", round(varpc[2], 2), " %)", sep = "")

# Plot extendido con etiquetas repelidas
pCOA_rscu_ext <- ggplot() +
  geom_point(data = cds_coords, aes(x = Dim.1, y = Dim.2, color = Group), size = 2) +
  geom_point(data = codon_coords, aes(x = Dim.1, y = Dim.2), color = "black", size = 1.5) +
  geom_text_repel(
    data = codon_coords,
    aes(x = Dim.1, y = Dim.2, label = codon),
    size = 3,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    segment.color = "gray60",
    min.segment.length = 0.1
  ) +
  scale_color_manual(values = c("conserved"="#00BFC4", "nonconserved"="#F8766D")) +
  theme_minimal() +
  labs(title = "RSCU - codons",
       x = xlabtext, y = ylabtext)

print(pCOA_rscu_ext)


