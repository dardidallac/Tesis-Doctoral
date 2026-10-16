import pandas as pd

# Archivos
input_file = "trnas_lens.txt"
output_file = "trna_abundance_lens.csv"

# Definir nombres de columnas (según tu salida)
colnames = [
    "Sequence", "tRNA#", "Begin", "End", "Type",
    "Anticodon", "IntronBegin", "IntronEnd", "Score", "Note"
]

# Leer el archivo, saltando las 2 primeras líneas (header + guiones)
df = pd.read_csv(
    input_file,
    sep=r"\s+",
    skiprows=2,
    names=colnames,
    engine="python"
)

# Contar abundancia de anticodones
anticodon_counts = df["Anticodon"].value_counts().reset_index()
anticodon_counts.columns = ["Anticodon", "Count"]

# Guardar en CSV
anticodon_counts.to_csv(output_file, index=False)

print(f"✅ Tabla de abundancia guardada en {output_file}")

