import pandas as pd
from codonbias.scores import TrnaAdaptationIndex
from Bio import SeqIO

# -------- CONFIGURACIÓN --------
species = {
    "Cicer": {
        "conserved": "Cicer_nonconserved.fa",
        "tGCN": "Cicer_tRNA_table_clean.csv"
    }
}

WINDOW_SIZE = 1000  # codones

# -------- PARTE 1: LIMPIEZA DE TABLA --------
raw_table = pd.read_csv("trna_abundance_cicer.csv")

# Renombrar columnas
raw_table = raw_table.rename(columns={"Anticodon": "anti_codon", "Count": "GCN"})

# Filtrar anticodones válidos (solo tripletes ACGT)
clean_table = raw_table[raw_table["anti_codon"].str.fullmatch(r"[ACGT]{3}")].copy()

# Agrupar duplicados
clean_table = clean_table.groupby("anti_codon", as_index=False)["GCN"].sum()

# Guardar tabla limpia
clean_table.to_csv("Cicer_tRNA_table_clean.csv", index=False)
print("✅ Tabla limpia guardada en Cicer_tRNA_table_clean.csv")
print(clean_table.head())

# -------- PARTE 2: FUNCIONES --------
def split_windows(seq, window_size=1000):
    """Divide una secuencia en ventanas de tamaño fijo de codones."""
    codons = [str(seq[i:i+3]) for i in range(0, len(seq), 3)]
    windows = []
    for i in range(0, len(codons), window_size):
        window = "".join(codons[i:i+window_size])
        if len(window) % 3 == 0 and len(window) > 0:
            windows.append(window)
    return windows

def analyze_species(name, files):
    """Analiza las secuencias conserved para TAI."""
    # Leer tabla de tRNAs
    tGCN = pd.read_csv(files["tGCN"])
    
    # Crear modelo TAI
    tai_model = TrnaAdaptationIndex(tGCN=tGCN)

    results = []
    for group in ["conserved"]:
        for rec in SeqIO.parse(files[group], "fasta"):
            seq = rec.seq
            windows = split_windows(seq, window_size=WINDOW_SIZE)

            for w_i, window in enumerate(windows, start=1):
                try:
                    tai_val = tai_model.get_score(window)
                except Exception as e:
                    print(f"⚠️ Error en {rec.id} ventana {w_i}: {e}")
                    tai_val = None

                results.append({
                    "Species": name,
                    "Sequence_ID": rec.id,
                    "Group": group.capitalize(),
                    "Window": w_i,
                    "TAI": tai_val
                })

    df = pd.DataFrame(results)
    out_file = f"{name}_TAI.xlsx"
    df.to_excel(out_file, index=False)
    print(f"✔ {name}: resultados guardados en {out_file}")

# -------- PARTE 3: EJECUCIÓN --------
for sp, files in species.items():
    analyze_species(sp, files)

