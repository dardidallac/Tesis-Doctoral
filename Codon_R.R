# Configuración inicial y carga de librerías

install.packages("devtools", dependencies = TRUE)
library(devtools)

devtools::install_github("mariodosreis/tai")

load.lib <- c('ggplot2', 'ggrepel', 'ggthemes', 
              'ggpubr', 'FactoMineR', 'factoextra', 'BiocManager', 'devtools')

# Instalación de librerías necesarias
install.lib <- load.lib[!load.lib %in% installed.packages()]
if (length(install.lib) > 0) {
  cat("Installing packages:", paste(install.lib, collapse = ", "), "\n")
  install.packages(install.lib, dependencies = TRUE)
}

# Cargar las librerías instaladas
suppressMessages(sapply(load.lib, require, character.only = TRUE))

# Instalación de paquetes adicionales
if (!"coRdon" %in% installed.packages()) {
  BiocManager::install("coRdon")
}
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools")
}
devtools::install_github("mariodosreis/tai")

# Selección del archivo con las secuencias CDS
cat("Selecciona el archivo con las secuencias codificantes (CDS):\n")
cds <- file.choose()
Filter_smaller_than_80_codons <- TRUE


#Subir archivos de referencia
cat("Selecciona el archivo con las secuencias codificantes ref (CDS):\n")
ref_cds <- file.choose()
Filter_smaller_than_80_codons <- TRUE


# Lectura y procesamiento de secuencias CDS
library('coRdon')

#Leer los cdss
# cds <- 'nombre_de_archivo'
coding_Sequences <- readSet(file = cds)

#Filter CDS with less than 80 codons
if (Filter_smaller_than_80_codons){
  coding_Sequences <- coding_Sequences[which(coding_Sequences@ranges@width >= 80*3)]
}

#count codons
cc_cds <- codonTable(coding_Sequences)

##########

# ref_cds <- 'nombre_de_archivo'
ref_coding_Sequences <- readSet(file = ref_cds)

#Filter CDS with less than 80 codons
if (Filter_smaller_than_80_codons){
  ref_coding_Sequences <- ref_coding_Sequences[which(ref_coding_Sequences@ranges@width >= 80*3)]
}

#count codons
ref_cc_cds <- codonTable(ref_coding_Sequences)
#############

#GC3
calcGC3 <- function(countsTable){
  counts <- countsTable@counts
  len <- countsTable@len
  GC3 <- counts[,grep("..G|..C",colnames(counts))]
  GC3.sum <- apply(as.matrix(GC3),1,sum)
  GC3.res <- GC3.sum/len*100
  return(GC3.res)
}

#GC12
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

#Calculate using coRdon and calcGC3 function
milc <- MILC(cc_cds)[,1]
enc <- ENC(cc_cds)
mcb <- MCB(cc_cds)[,1]
encp <- ENCprime(cc_cds, subsets = list(cc_cds))[,2]
scuo <- SCUO(cc_cds)

gc12 <- calcGC12(cc_cds)
gc3 <- calcGC3(cc_cds)

table <- data.frame(cds = coding_Sequences@ranges@NAMES, MILC = milc, Enc = enc, MCB = mcb, EncPrime = encp, SCUO = scuo, GC3 = gc3)
print(head(table))
write.csv(table,'table.csv')

#Make plot

library('ggplot2')
library('ggthemes')
library('data.table')

#wide to long table
ltable <- melt(setDT(table), id.vars = c('cds'), variable.name = "Index")
milc_table <- ltable[which(ltable$'Index' == 'MILC'), ]
Enc_table <- ltable[which(ltable$'Index' == 'Enc'), ]
MCB_table <- ltable[which(ltable$'Index' == 'MCB'), ]
EncPrime_table <- ltable[which(ltable$'Index' == 'EncPrime'), ]
SCUO_table <- ltable[which(ltable$'Index' == 'SCUO'), ]
GC3_table <- ltable[which(ltable$'Index' == 'GC3'), ]

box_plot_MILC <- ggplot(milc_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("MILC")

box_plot_Enc <- ggplot(Enc_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("Enc")

box_plot_MCB <- ggplot(MCB_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("MCB")

box_plot_EncPrime <- ggplot(EncPrime_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("Enc'")

box_plot_SCUO <- ggplot(SCUO_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("SCUO")

box_plot_GC3<- ggplot(GC3_table, aes(x=Index, y=value, color = Index)) +
  geom_boxplot() +
  theme(legend.position='none') +
  ylab("GC3")

print(box_plot_MILC)
print(box_plot_Enc)
print(box_plot_MCB)
print(box_plot_EncPrime)
print(box_plot_SCUO)
print(box_plot_GC3)

library('ggplot2')
library('ggthemes')

#NC plot
fgc3s <- function (gc3) {
  a = -6
  b = 34
  c = 1.025
  x = a + gc3 + (b/(gc3^2 + (c - gc3)^2))
  return(x)
}

ncplot <- data.frame(x=seq(0,1,0.01), y=fgc3s(seq(0,1,0.01)))
ncplot[,1] <- ncplot[,1]*100

print(head(ncplot))

pNCPlot <- ggplot(data=table, aes(x=GC3, y=Enc)) +
  geom_point(colour='black')+
  geom_line(data=ncplot, aes(x=x,y=y), color='red')+
  ylab("Nc") +
  xlab("GC3s") +
  theme_classic() +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=18, face="bold"), #element_blank()
        axis.title.y = element_text(size=18, face="bold"))

print(pNCPlot)

#Neutrality plot
#GC12 vs GC3, correlation
library("ggpubr")
gc12_gc3 <- data.frame(GC12 = gc12, GC3 = gc3)

neutrality_plot <- ggscatter(gc12_gc3, x = "GC3", y = "GC12", add = "reg.line",
                             conf.int = TRUE,
                             #cor.coef = TRUE,
                             color = "dodgerblue4", #fill = "dodgerblue",
                             cor.method = "pearson", xlab = "GC3 (%)", ylab = "GC12 (%)") +
  stat_cor(, label.x = min(gc3), label.y = max(gc12)+(max(gc12)-min(gc12))/20) +
  stat_regline_equation(label.x = min(gc3), label.y = max(gc12)+2*(max(gc12)-min(gc12))/20  )+
  theme(legend.position="none")


print(neutrality_plot)

library("FactoMineR")
library("factoextra")
library("ggrepel")
#COA with raw counts [Tiene mas sentido para ver diferentres genes agrupados por tejido, nivel de expresión, etc...] Falta calcular un centro para
#diferentes conjuntos de datos.
#cc_cds
#discard stop codons and M W
cc4COA <- cc_cds@counts
cc4COA = cc4COA[,-c(15,49,51,57,59)]
res.ca <- CA(cc4COA, graph = T)

#variance percent
varpc <- get_eigenvalue(res.ca)[1:2,2]
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")
cds_pcs <- data.frame(get_ca_row(res.ca)$coord[,1:2])
cds_pcs <- cbind(cds_pcs, SET=rep("CDS",dim(cds_pcs)[1]))
codons <- data.frame(get_ca_col(res.ca)$coord[,1:2])

#COA plot
coaplot <- ggplot (cds_pcs, aes( y = Dim.2, x = Dim.1)) +
  geom_point(data = cds_pcs, color = "grey50", size=1) +
  xlab(xlabtext) +
  ylab(ylabtext) +
  labs(title = "Raw codon counts",fill = "Legend")  +
  theme(title = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"),
        legend.position='none')

print(coaplot)

#RSCU
counts_cds <- cc_cds@counts
aminoacids <- c('K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P',
                'R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
                '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F')
aas <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
codons = colnames(counts_cds)
dnacodons <- data.frame(codon = colnames(counts_cds), aa = aminoacids)

#codones con cuenta 0 total para una aminoácido los cuento como 1, random
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
#Average codon usage
rscu = rscu[,-c(15,49,51,57,59)]
avg_rscu <- colMeans(rscu)
rscu4COA <- rbind.data.frame(rscu,avg_rscu)
ncds <- dim(rscu)[1]
totalrows <- dim(rscu4COA)[1]

res.ca <- CA(rscu4COA, graph = T, row.sup = (ncds+1):totalrows)
#variance percent
varpc <- get_eigenvalue(res.ca)[1:2,2]
xlabtext <- paste("C1 (",round(varpc[1],2)," %)", sep="")
ylabtext <- paste("C2 (",round(varpc[2],2)," %)", sep="")

cds_pcs <- data.frame(get_ca_row(res.ca)$coord[,1:2])
cds_pcs <- cbind(cds_pcs, SET=rep("CDS",dim(cds_pcs)[1]))
codons <- data.frame(get_ca_col(res.ca)$coord[,1:2])
#COA plot
coaplot <- ggplot (cds_pcs, aes( y = Dim.2, x = Dim.1)) +
  geom_point(data = cds_pcs, color = "grey50", size=1) +
  xlab(xlabtext) +
  ylab(ylabtext) +
  labs(title = "RSCU",fill = "Legend")  +
  theme(title = element_text(size=16, face="bold"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"),
        legend.position='none')

print(coaplot)


#CAI

cai <- CAI(cc_cds, subsets = list(ref_cc_cds))[,1]

table <- data.frame(cds = coding_Sequences@ranges@NAMES, MILC = milc, Enc = enc, MCB = mcb, EncPrime = encp, SCUO = scuo, GC3 = gc3, CAI = cai)
print(head(table))
write.csv(table,'table.csv')

#MELP
melp <- MELP(cc_cds, subsets = list(ref_cc_cds))[,1]
head(melp)

table <- data.frame(cds = coding_Sequences@ranges@NAMES, MILC = milc, Enc = enc, MCB = mcb, EncPrime = encp, SCUO = scuo, GC3 = gc3, CAI = cai, MELP = melp)
print(head(table))
write.csv(table,'table.csv')

install.packages("tAI") # Si el paquete está disponible en CRAN
library(tAI)


#tAI
#Tabla de trna: columnas: codon trna , donde codon es NNN y trna es número de trnas para ese codon

wCodonOrder <- c('TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG','CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG','ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG','GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG')

wCodonOrder2 <- wCodonOrder[-c(11,12,15)]

trna <- read.table("~/Escritorio/uso de codones tesis/tRNATable_medicago.csv", sep = " ", header = TRUE)
trna <- trna[match(wCodonOrder,trna$codon),]
ws <- get.ws(tRNA=trna$trna,s=NULL, sking=0)
print (ws)

#order codon counts
codonCountdf <- as.data.frame(cc_cds@counts)
codonCountdf <- codonCountdf[,match(wCodonOrder,colnames(codonCountdf))][,-c(11,12,15)]
#ws could be used as an indirect measure of translational speed
#33 is methionine once stop codons have been removed
tai <- get.tai(codonCountdf[,-33],ws) #calculate tAIs
print (tai)




# Cargar tu tabla 'wsTable.csv'
wsTable <- read.csv("~/Escritorio/uso de codones tesis/wsTable_medicago.csv")


ws <- wsTable$ws
tai <- get.tai(codonCountdf[,-33],ws)
print(tai)







table <- data.frame(cds = coding_Sequences@ranges@NAMES, MILC = milc, Enc = enc, MCB = mcb, EncPrime = encp, SCUO = scuo, GC3 = gc3, CAI = cai, MELP = melp, TAI = tai)
print(head(table))
write.csv(table,'table.csv')

write.csv(table, file = "~/Escritorio/table.csv", row.names = FALSE)


