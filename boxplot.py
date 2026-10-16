import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from statannotations.Annotator import Annotator

# -------- CONFIGURACIÓN --------
archivos = [
    "Cicer_CAI_ENC.xlsx",
    "Lens_CAI_ENC.xlsx",
    "Lotus_CAI_ENC.xlsx",
    "Medicago_CAI_ENC.xlsx",
    "Phaseolus_CAI_ENC.xlsx"
]

carpeta = "/home/mario/Escritorio/Codons_orthologues/codon_bias_80"

# -------- LEER Y CONCATENAR DATOS --------
dfs = []
for archivo in archivos:
    path = os.path.join(carpeta, archivo)
    df = pd.read_excel(path)
    
    # Normalizar decimales
    df["CAI"] = df["CAI"].astype(str).str.replace(",", ".").astype(float)
    df["ENC"] = df["ENC"].astype(str).str.replace(",", ".").astype(float)
    
    # Agregar nombre del archivo como columna "Species"
    especie = archivo.split("_")[0]
    df["Species"] = especie
    dfs.append(df)

df_all = pd.concat(dfs, ignore_index=True)

# -------- GRAFICAR CON TESTS --------
sns.set(style="whitegrid")

# Definir pares de comparación: Conserved vs Nonconserved para cada especie
pairs = [( (sp, "Conserved"), (sp, "Nonconserved")) for sp in df_all["Species"].unique()]

# === CAI ===
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=df_all, x="Species", y="CAI", hue="Group", palette="Set2")

annotator = Annotator(ax, pairs, data=df_all, x="Species", y="CAI", hue="Group")
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

plt.title("Distribución de CAI por especie (Conserved vs Nonconserved)")
plt.legend(title="Group")
plt.tight_layout()
plt.savefig("boxplot_CAI_all_stats.png", dpi=300)
plt.show()

# === ENC ===
plt.figure(figsize=(10, 6))
ax = sns.boxplot(data=df_all, x="Species", y="ENC", hue="Group", palette="Set2")

annotator = Annotator(ax, pairs, data=df_all, x="Species", y="ENC", hue="Group")
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

plt.title("Distribución de ENC por especie (Conserved vs Nonconserved)")
plt.legend(title="Group")
plt.tight_layout()
plt.savefig("boxplot_ENC_all_stats.png", dpi=300)
plt.show()

