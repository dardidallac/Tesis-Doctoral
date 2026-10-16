import os
import subprocess

# Directorios
input_dir = "/home/mario/Escritorio/Codons_orthologues/Single_Copy_Orthologue_Sequences"  # Ortogrupos de proteínas
output_dir = "Aligned_Orthogroups"
os.makedirs(output_dir, exist_ok=True)

# Opciones de MAFFT (puedes ajustarlas según necesidad)
mafft_options = "--auto"  # O puedes usar --maxiterate 1000 --localpair para más precisión

# Procesar cada archivo FASTA
for file in os.listdir(input_dir):
    if not file.endswith(".fa"):
        continue
    input_path = os.path.join(input_dir, file)
    output_path = os.path.join(output_dir, file)
    
    print(f"[INFO] Alineando {file}...")
    
    # Ejecutar MAFFT
    cmd = f"mafft {mafft_options} {input_path} > {output_path}"
    subprocess.run(cmd, shell=True)
    
print("\n✅ Alineamientos completos en la carpeta:", output_dir)

