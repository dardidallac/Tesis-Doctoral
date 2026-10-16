import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from statannotations.Annotator import Annotator

# -------- CONFIGURACIÓN --------
especies = ["Cicer", "Lotus","Phaseolus"]
carpeta = "/home/mario/Escritorio/Codons_orthologues/codon_bias_80"

# -------- LEER Y CONCATENAR DATOS --------
dfs = []

for especie in especies:
    for group in ["conserved", "nonconserved"]:
        archivo = f"{especie}_TAI_{group}.xlsx"
        path = os.path.join(carpeta, archivo)
        df = pd.read_excel(path)
        
        # Normalizar decimales
        df["TAI"] = df["TAI"].astype(str).str.replace(",", ".").astype(float)
        
        # Agregar columnas necesarias
        df["Species"] = especie
        df["Group"] = "Conserved" if group == "conserved" else "Nonconserved"
        
        dfs.append(df)

df_all = pd.concat(dfs, ignore_index=True)

# -------- GRAFICAR CON TESTS --------
sns.set(style="whitegrid")

# Definir pares Conserved vs Nonconserved por especie
pairs = [((sp, "Conserved"), (sp, "Nonconserved")) for sp in df_all["Species"].unique()]

# === TAI ===
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=df_all,
    x="Species",
    y="TAI",
    hue="Group",
    palette="Set2",
    dodge=True
)

annotator = Annotator(ax, pairs, data=df_all, x="Species", y="TAI", hue="Group")
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

plt.title("Distribución de TAI por especie (Conserved vs Nonconserved)")
plt.legend(title="Group")
plt.tight_layout()
plt.savefig("boxplot_TAI_all_stats.png", dpi=300)
plt.show()

