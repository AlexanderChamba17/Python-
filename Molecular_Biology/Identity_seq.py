# Calculo de porcentaje de identidad entre secuencias
title = """
La siguiente calculadora determina el porcentaje de identidad 
entre dos secuancias de ADN, si la secuencia tiene GAPS marcarlos 
con un guión (-). Ejemplo ATG-CATG- 
Nota. las dos secuencias deben tener la misma longitud de caracteres
para el análisis.
"""
print(title)

#Definición de variables

q1_seq1 = input("Ingresar la primera secuencia: ")
q2_seq2 = input("Ingresar la segunda secuencia: ")

#Función de % identidad

def ident_seq(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Las secuencias NO tienen la misma longitud")
    
    identidad = sum(1 for base1, base2 in zip(seq1, seq2) if base1 == base2)
    porcentage = (identidad /  len(seq1) * 100)
    return porcentage

try: 
    porcnt_result = ident_seq(q1_seq1, q2_seq2)
    print(f"El porcentaje de identidad entre las secuencias ingresadas es: {porcnt_result: .2f}%")

except ValueError as e:
    print(f"Error: {e}")

