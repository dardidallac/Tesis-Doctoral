import os
from Bio import SeqIO

input_dir = "/home/mario/Escritorio/Codons_orthologues/Single_Copy_Orthologue_CDS"
output_dir = "Filtered_Orthogroups_cds"
os.makedirs(output_dir, exist_ok=True)

# Contador de archivos filtrados
filtered_count = 0

for file in os.listdir(input_dir):
    if not file.endswith(".fa"):
        continue
    input_path = os.path.join(input_dir, file)
    sequences = list(SeqIO.parse(input_path, "fasta"))
    
    if len(sequences) == 5:  # Exactamente 5 secuencias
        output_path = os.path.join(output_dir, file)
        SeqIO.write(sequences, output_path, "fasta")
        filtered_count += 1
        print(f"[FILTERED] {file} -> {len(sequences)} secuencias")

print(f"\nTotal archivos filtrados con exactamente 5 secuencias: {filtered_count}")

