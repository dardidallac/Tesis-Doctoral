install.packages("readxl")  # Si no lo tienes instalado
library(readxl)


tabla_expresion <- read_excel("~/Escritorio/uso de codones tesis/tabla_expresion.xlsx")





cor_test_result <- cor.test(tabla_expresion$root_symbiotic_condition, tabla_expresion$MILC)


print(cor_test_result)
