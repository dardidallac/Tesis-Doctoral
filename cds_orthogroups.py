import os
import re
from Bio import SeqIO

# ==============================
# CONFIGURACIÓN DE RUTAS
# ==============================
cds_dir = "/home/mario/Escritorio/Codons_orthologues"  # Carpeta con los .fa de CDS
orthogroups_dir = "/home/mario/Escritorio/Codons_orthologues/Single_Copy_Orthologue_Sequences"
output_dir = "Single_Copy_Orthologue_CDS"
log_file = "missing_ids.log"

os.makedirs(output_dir, exist_ok=True)

# Archivos CDS por especie
cds_files = {
    "vulgaris": "cds_vulgaris.fa",
    "cicer": "cicar.CDCFrontier.gnm3.ann1.NPD7.cds.fa",
    "lens": "Lculinaris_718_v1.cds.fa",
    "lotus": "Ljaponicus_571_Lj1.0v1.cds.fa",
    "medicago": "Mtruncatula_cds.fa"
}

# ==============================
# FUNCIÓN PARA EXTRAER IDs
# ==============================
def extract_id(header, species):
    if species == "vulgaris":
        # Captura protein_id=ABC123.1 en cualquier formato
        match = re.search(r"protein_id=([\w\.\-]+)", header)
        return match.group(1) if match else header.split()[0]
    elif species == "cicer":
        return header.split()[0]
    elif species == "lotus":
        # Priorizar polypeptide= (coincide con Orthofinder)
        match = re.search(r"polypeptide=(\S+)", header)
        if match:
            return match.group(1)
        # Si no hay, usar ID=
        match = re.search(r"ID=([^\s]+)", header)
        return match.group(1) if match else header.split()[0]
    elif species == "lens":
        return header.split()[0]
    elif species == "medicago":
        return header.split()[0]
    return None

# ==============================
# CARGAR TODOS LOS CDS EN DICCIONARIO
# ==============================
cds_dict = {}
species_counts = {}

print("=== Cargando CDS ===")
for sp, filename in cds_files.items():
    path = os.path.join(cds_dir, filename)
    print(f"[INFO] Cargando {filename} ({sp})...")
    count = 0
    for record in SeqIO.parse(path, "fasta"):
        pid = extract_id(record.description, sp)
        if pid:
            cds_dict[pid] = str(record.seq)
            count += 1
    species_counts[sp] = count
    print(f"   -> {count} secuencias cargadas para {sp}")

print(f"TOTAL CDS cargados: {len(cds_dict)}")
print(species_counts)

# ==============================
# PROCESAR ORTHOGROUPS
# ==============================
missing_ids = []

print("\n=== Procesando ortogrupos ===")
for file in os.listdir(orthogroups_dir):
    if not file.endswith(".fa"):
        continue
    input_path = os.path.join(orthogroups_dir, file)
    output_path = os.path.join(output_dir, file)
    output_seqs = []
    
    for record in SeqIO.parse(input_path, "fasta"):
        prot_id = record.id  # ID del ortogrupo
        if prot_id in cds_dict:
            seq = cds_dict[prot_id]
            output_seqs.append(f">{prot_id}\n{seq}")
        else:
            missing_ids.append(f"{prot_id}\t{file}")
            print(f"[WARNING] No se encontró CDS para {prot_id} en {file}")

    if output_seqs:
        with open(output_path, "w") as out:
            out.write("\n".join(output_seqs))

# ==============================
# GUARDAR LOG DE FALTANTES
# ==============================
if missing_ids:
    with open(log_file, "w") as log:
        log.write("\n".join(missing_ids))
    print(f"\n[WARNING] {len(missing_ids)} IDs no encontrados. Guardado en {log_file}")
else:
    print("\nTodos los IDs fueron encontrados correctamente.")

