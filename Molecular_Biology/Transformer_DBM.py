# Genetic code transformation
# DNA to RNA function
def DNA_to_RNA(DNAseq):
    return DNAseq.replace("T", "U")


# RNA to DNA function 
def RNA_to_DNA(RNAseq):
    return RNAseq.replace("U", "T")

# Genetic code

# DNA to protein 
# Codons_DNA to aminoacids (Dictonary)
genetic_code_DNA_to_Protein = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# DNA to protein function
def DNA_to_protein(DNAseq):
    proteinSeq = ""
    for codon in range(0, len(DNAseq), 3):
        codon = DNAseq[codon : codon + 3]
        amino_acid = genetic_code_DNA_to_Protein.get(codon, "X")
        proteinSeq += amino_acid
    return proteinSeq


# RNA to protein 
# Codons_RNA to aminoacids (Dictonary)
genetic_code_RNA_to_Protein = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# RNA to protein function 
def RNA_to_protein(RNAseq):
    proteinSeq = ""
    for codon in range(0, len(RNAseq), 3):
        codon = RNAseq[codon : codon + 3]
        amino_acid = genetic_code_RNA_to_Protein.get(codon, "X")
        proteinSeq += amino_acid
    return proteinSeq
    

# Main function
def main():
    q1 = input("Ingresar la secuencia de ADN, ARN o proteina: ").upper()

    while True:
        q2 = input("La secuencia ingresada es ADN (D), ARN (R) o proteina (P)?: ").upper()
        q3 = input("A que desea convertir la secuencia? ADN (D), ARN (R) o proteina (P): ").upper()

        if q2 == "D" and q3 == "R":
            converted_sequence = DNA_to_RNA(q1)
            break
        elif q2 == "R"  and q3 == "D":
            converted_sequence = RNA_to_DNA(q1)
            break
        elif q2 == "D" and q3 == "P":
            converted_sequence = DNA_to_protein(q1)
            break
        elif q2 == "R" and q3 == "P":
            converted_sequence = RNA_to_protein(q1)
            break 
        else: 
            print("Datos incorrectos")
    print(f"La secuencia convertida es: {converted_sequence}")


if __name__ == "__main__":
    main()
    