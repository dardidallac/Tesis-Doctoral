from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def count_snps(alignment):
    snp_count = 0
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        if len(set(column)) > 1:
            snp_count += 1
    return snp_count

# Cambia 'alineamiento.fasta' al nombre de tu archivo de alineamiento múltiple en formato FASTA
alignment_file = 'alineamiento.fasta'

try:
    alignment = AlignIO.read(alignment_file, 'fasta')
    snp_count = count_snps(alignment)
    print(f"La cantidad de SNPs en el alineamiento es: {snp_count}")
except FileNotFoundError:
    print(f"El archivo '{alignment_file}' no se encuentra.")
except ValueError:
    print("El archivo no parece ser un alineamiento múltiple en formato FASTA válido.")


