def leer_fasta(archivo):
    secuencias = {}
    with open(archivo, 'r') as f:
        nombre = None
        secuencia = ''
        for linea in f:
            if linea.startswith('>'):
                if nombre:
                    secuencias[nombre] = secuencia
                nombre = linea.strip()[1:]
                secuencia = ''
            else:
                secuencia += linea.strip()
        if nombre:
            secuencias[nombre] = secuencia
    return secuencias

def calcular_proporcion_snps(secuencias):
    if 'P.vulgarisG23939/1-1722' not in secuencias:
        raise ValueError("La secuencia 'original_cds' no se encuentra en las secuencias proporcionadas.")
    
    original_cds = secuencias['P.vulgarisG23939/1-1722']
    longitud_secuencia = len(original_cds)
    cantidad_secuencias = len(secuencias)
    propor_snps_posiciones = []

    for i in range(longitud_secuencia):
        secuencias_con_snps = 0
        for nombre, secuencia in secuencias.items():
            if nombre != 'P.vulgarisG23939/1-1722' and secuencia[i] != original_cds[i]:
                secuencias_con_snps += 1
        if secuencias_con_snps > 0:
            propor_snps = (secuencias_con_snps / (cantidad_secuencias - 1)) * 100
            propor_snps_posiciones.append((i + 1, propor_snps, secuencias_con_snps))

    return propor_snps_posiciones

def main():
    archivo_andinos = 'porotos_andinos.fasta'
    archivo_mesoamericanos = 'porotos_mesoamericanos.fasta'

    secuencias_andinos = leer_fasta(archivo_andinos)
    secuencias_mesoamericanos = leer_fasta(archivo_mesoamericanos)

    print("Proporción de SNPs en porotos andinos:")
    for posicion, propor_snps, num_secuencias_con_snps in calcular_proporcion_snps(secuencias_andinos):
        print(f"Posición {posicion}: Proporción {propor_snps:.2f}%, Secuencias con SNPs {num_secuencias_con_snps}")

    print("Proporción de SNPs en porotos mesoamericanos:")
    for posicion, propor_snps, num_secuencias_con_snps in calcular_proporcion_snps(secuencias_mesoamericanos):
        print(f"Posición {posicion}: Proporción {propor_snps:.2f}%, Secuencias con SNPs {num_secuencias_con_snps}")

if __name__ == "__main__":
    main()

