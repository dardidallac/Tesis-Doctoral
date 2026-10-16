import os
from Bio import AlignIO
import pandas as pd

def calcular_matriz_snps(andino_fasta, meso_fasta, ref_keyword="original_cds"):
    """
    Lee los alineamientos de ambos acervos, identifica la secuencia de referencia,
    cuenta los SNPs por posición y calcula los porcentajes de frecuencia.
    """
    # 1. Cargar los alineamientos hechos con MAFFT
    print("Cargando alineamientos...")
    align_andino = AlignIO.read(andino_fasta, "fasta")
    align_meso = AlignIO.read(meso_fasta, "fasta")
    
    # Obtener el largo total del alineamiento (asumiendo que coinciden por la referencia)
    align_len = align_andino.get_alignment_length()
    
    # 2. Determinar el total de muestras reales por grupo (restando la referencia si está dentro)
    total_andino = len(align_andino)
    total_meso = len(align_meso)
    
    # 3. Buscar la secuencia de referencia para usarla como molde de comparación
    ref_seq = None
    for record in align_andino:
        if ref_keyword in record.id.lower():
            ref_seq = str(record.seq)
            total_andino -= 1  # No contar la referencia en el porcentaje poblacional
            break
            
    if not ref_seq:
        for record in align_meso:
            if ref_keyword in record.id.lower():
                ref_seq = str(record.seq)
                total_meso -= 1
                break
                
    if not ref_seq:
        raise ValueError(f"No se encontró ninguna secuencia que contenga '{ref_keyword}' en el ID para usar como referencia.")

    print(f"Muestras a procesar -> Andinas: {total_andino} | Mesoamericanas: {total_meso}")
    
    snp_records = []
    
    # 4. Recorrer el alineamiento columna por columna
    for pos in range(align_len):
        ref_base = ref_seq[pos].upper()
        
        # Ignorar si la propia referencia tiene un gap en esa posición del alineamiento
        if ref_base == "-":
            continue
            
        # Extraer las bases de esa columna para cada población
        col_andino = [seq[pos].upper() for seq in align_andino if ref_keyword not in seq.id.lower()]
        col_meso = [seq[pos].upper() for seq in align_meso if ref_keyword not in seq.id.lower()]
        
        # Contar cuántas secuencias difieren de la referencia (ignorando gaps y ambigüedades 'N')
        snps_andino = sum(1 for base in col_andino if base != ref_base and base not in ["-", "N"])
        snps_meso = sum(1 for base in col_meso if base != ref_base and base not in ["-", "N"])
        
        total_snps = snps_andino + snps_meso
        
        # Si se detectó al menos una mutación en esta coordenada, guardamos los datos
        if total_snps > 0:
            # Calcular porcentajes relativos al tamaño de su propia población
            pct_andino = (snps_andino / total_andino) * 100 if total_andino > 0 else 0
            pct_meso = (snps_meso / total_meso) * 100 if total_meso > 0 else 0
            
            snp_records.append({
                "posicion": pos + 1,  # Coordenada 1-based (estilo Jalview/Genética)
                "ref_base": ref_base,
                "andino_snps": snps_andino,
                "andino_pct": round(pct_andino, 2),
                "meso_snps": snps_meso,
                "meso_pct": round(pct_meso, 2),
                "total_snps": total_snps
            })
            
    # Convertir a DataFrame de Pandas
    df_snps = pd.DataFrame(snp_records)
    return df_snps

# --- CONFIGURACIÓN DE RUTAS ---
# Reemplaza con los nombres exactos de tus archivos de alineamiento intermedios
andino_file = "/home/mario/Escritorio/SNPs/secuencias_andinas_nfr1_mesoref.fa"
meso_file = "/home/mario/Escritorio/SNPs/secuencias_mesoamericanas_nfr1_mesoref.fa"
output_csv = "reporte_conteo_y_porcentaje_snps_nfr1.csv"

if __name__ == "__main__":
    try:
        df_resultado = calcular_matriz_snps(andino_file, meso_file, ref_keyword="original_cds")
        
        # Guardar a un archivo CSV compatible con Excel
        df_resultado.to_csv(output_csv, index=False)
        print(f"\n¡Proceso completado con éxito!")
        print(f"Se encontraron {len(df_resultado)} posiciones polimórficas.")
        print(f"Archivo guardado en: {os.path.abspath(output_csv)}")
        
        # Mostrar las primeras filas como vista previa
        print("\nVista previa del resultado:")
        print(df_resultado.head(10).to_string(index=False))
        
    except Exception as e:
        print(f"\n❌ Error al procesar: {e}")
