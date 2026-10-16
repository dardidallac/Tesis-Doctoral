import os
import subprocess

# Directorios
protein_aligned_dir = "/home/mario/Escritorio/Codons_orthologues/Aligned_Orthogroups_proteinas"
cds_dir = "/home/mario/Escritorio/Codons_orthologues/Filtered_Orthogroups_cds"
output_dir = "Codon_Aligned_Orthogroups"
os.makedirs(output_dir, exist_ok=True)

# Opciones de PAL2NAL
pal2nal_options = "-output fasta"  # También puede ser clustal, codon, paml, etc.

for file in os.listdir(protein_aligned_dir):
    if not file.endswith(".fa"):
        continue
    
    protein_path = os.path.join(protein_aligned_dir, file)
    cds_path = os.path.join(cds_dir, file)
    
    if not os.path.exists(cds_path):
        print(f"[WARNING] No se encontró el CDS para {file}, se omite.")
        continue
    
    output_path = os.path.join(output_dir, file)
    
    print(f"[INFO] Procesando {file}...")
    
    # Comando PAL2NAL
    cmd = f"pal2nal.pl {protein_path} {cds_path} {pal2nal_options} > {output_path}"
    subprocess.run(cmd, shell=True)

print("\n✅ Alineamientos codon-aware completos en:", output_dir)

