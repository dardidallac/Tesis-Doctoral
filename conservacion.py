import pandas as pd
from Bio import AlignIO
from collections import Counter
import math
import os

def calcular_conservacion_identidad(alineamiento):
    """
    Calcula el grado de conservación por posición basado en la identidad simple.
    """
    conservacion_por_posicion = []
    num_secuencias = len(alineamiento)
    longitud_alineamiento = alineamiento.get_alignment_length()

    for i in range(longitud_alineamiento):
        columna = alineamiento[:, i]
        frecuencias = Counter(columna)

        # Ignorar los gaps para el cálculo de identidad
        if '-' in frecuencias:
            del frecuencias['-']

        if frecuencias:
            caracter_mas_frecuente, count = frecuencias.most_common(1)[0]
            grado_conservacion = (count / num_secuencias) * 100
        else:
            grado_conservacion = 0

        conservacion_por_posicion.append(grado_conservacion)
        
    return conservacion_por_posicion

def calcular_conservacion_shannon_entropy(alineamiento):
    """
    Calcula la entropía de Shannon para cada posición del alineamiento.
    Un valor bajo de entropía indica alta conservación.
    """
    entropia_por_posicion = []
    num_secuencias = len(alineamiento)
    longitud_alineamiento = alineamiento.get_alignment_length()
    
    for i in range(longitud_alineamiento):
        columna = alineamiento[:, i]
        frecuencias = Counter(columna)
        
        entropia = 0
        for caracter in frecuencias:
            probabilidad = frecuencias[caracter] / num_secuencias
            if probabilidad > 0:
                entropia += -probabilidad * math.log2(probabilidad)
        
        entropia_por_posicion.append(entropia)
        
    return entropia_por_posicion

# --- Main script ---

# Directorio que contiene los archivos de alineamiento
directorio_alineamientos = "/home/mario/Escritorio/Codons_orthologues/Aligned_Orthogroups_proteinas"

# Verificar si el directorio existe
if not os.path.isdir(directorio_alineamientos):
    print(f"Error: El directorio '{directorio_alineamientos}' no existe.")
else:
    # Recorrer todos los archivos en el directorio
    for nombre_archivo in os.listdir(directorio_alineamientos):
        # Procesar solo los archivos que terminan en .fa
        if nombre_archivo.endswith('.fa'):
            ruta_archivo_alineamiento = os.path.join(directorio_alineamientos, nombre_archivo)
            
            # Crear el nombre del archivo de salida de Excel
            # Cambia la extensión .fa por .xlsx
            nombre_excel = os.path.splitext(nombre_archivo)[0] + ".xlsx"
            ruta_archivo_salida_excel = os.path.join(directorio_alineamientos, nombre_excel)
            
            print(f"Procesando el archivo: {nombre_archivo}...")
            
            try:
                # Cargar el alineamiento
                alineamiento = AlignIO.read(ruta_archivo_alineamiento, "fasta")

                # Calcular ambos tipos de conservación
                conservacion_identidad = calcular_conservacion_identidad(alineamiento)
                entropia_shannon = calcular_conservacion_shannon_entropy(alineamiento)

                # Crear DataFrames de pandas con los resultados
                posiciones = list(range(1, alineamiento.get_alignment_length() + 1))
                
                df_identidad = pd.DataFrame({
                    'Posición': posiciones,
                    'Conservación (%)': conservacion_identidad
                })

                df_shannon = pd.DataFrame({
                    'Posición': posiciones,
                    'Entropía de Shannon': entropia_shannon
                })

                # Exportar los DataFrames a un archivo de Excel con dos hojas
                with pd.ExcelWriter(ruta_archivo_salida_excel, engine='openpyxl') as writer:
                    df_identidad.to_excel(writer, sheet_name='Identidad Simple', index=False)
                    df_shannon.to_excel(writer, sheet_name='Entropía de Shannon', index=False)
                
                print(f"-> Los resultados han sido exportados a '{ruta_archivo_salida_excel}' con éxito.")
            
            except Exception as e:
                print(f"Error al procesar '{nombre_archivo}': {e}")
                print("-" * 50)
    
    print("\nProceso de todos los archivos finalizado.")
