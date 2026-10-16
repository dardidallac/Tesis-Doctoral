from codonbias.scores import CodonAdaptationIndex, EffectiveNumberOfCodons
from Bio import SeqIO
import pandas as pd

# -------- CONFIGURACIÓN --------
species = {
    "Phaseolus": {
        "conserved": "Phaseolus_conserved.fa",
        "nonconserved": "Phaseolus_nonconserved.fa",
        "ribosomal": "vulgaris_ribosomal_cds.fa"
    }
}

WINDOW_SIZE = 1000  # codones

# -------- FUNCIONES --------
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
    """Analiza conserved y nonconserved de una especie."""
    # Construir modelo de CAI usando secuencias ribosomales
    ribosomal_records = list(SeqIO.parse(files["ribosomal"], "fasta"))
    ribosomal_seqs = [str(rec.seq) for rec in ribosomal_records]
    cai_model = CodonAdaptationIndex(ribosomal_seqs)

    # Modelo de ENC (no necesita referencia)
    enc_model = EffectiveNumberOfCodons()

    results = []
    for group in ["conserved", "nonconserved"]:
        for rec in SeqIO.parse(files[group], "fasta"):
            seq = rec.seq
            windows = split_windows(seq, window_size=WINDOW_SIZE)

            for w_i, window in enumerate(windows, start=1):
                try:
                    cai_val = cai_model.get_score(window)
                    enc_val = enc_model.get_score(window)
                except Exception:
                    cai_val, enc_val = None, None

                results.append({
                    "Species": name,
                    "Sequence_ID": rec.id,
                    "Group": group.capitalize(),
                    "Window": w_i,
                    "CAI": cai_val,
                    "ENC": enc_val
                })

    df = pd.DataFrame(results)
    out_file = f"{name}_CAI_ENC.xlsx"
    df.to_excel(out_file, index=False)
    print(f"✔ {name}: resultados guardados en {out_file}")

# -------- EJECUCIÓN --------
for sp, files in species.items():
    analyze_species(sp, files)

