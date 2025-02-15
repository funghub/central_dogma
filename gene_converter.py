# convert dna to rna to protein and dna to rna and back 
# examples TAUGCACCGTCACAC  AUGUCUUGAGG   TAGCATGATGAAATAAA

import re

rna_letters = '^[AGCU]*$'
dna_letters = '^[AGCT]*$'

# ask user for a sequence
def user_start():
    user_input = input('Type a RNA/DNA/amino acid sequence for a gene assuming coding sequence\n')
    user_input = user_input.upper()
    return user_input

def protein_converter(strand):
    import pandas as pd
    import math

    codons = pd.read_csv('codon_chart.csv')
    codons_dict = codons.to_dict(orient='list')

    for key, codon_list in codons_dict.items():
        indicies = [i for i,codon in enumerate(codon_list) if isinstance(codon, float) and math.isnan(codon)]
        # enumerate the location of nan in the list so that we can remove it in the codon_list for the specific key

        # create that list of where nan appears in the codon_list for the specific key
        [codon_list.pop(index) for index in reversed(indicies)]
        # removing an item from the list shifts all subsequent elements 
        # left by one position, and you may unintentionally skip over the element 
        # immediately following a removed one due to nature of index

        codons_dict[key] = codon_list


    ask_protein = input(
        'Do you want to convert this strand to protein with/without Start and Stop detection?\n'
        'A: Yes, but do not detect Start or Stop\n'
        'B: Yes, please detect both Start and Stop\n'
        'C: No\n')
    
    # DOES NOT DETECT START OR STOP
    if ask_protein in ("A", "a"):
        status = 'A'
        list_strand = re.findall('...', strand)

    # DOES DETECT START OR STOP
    elif ask_protein in ("B", "b"):
        status = 'B'
        start_location = strand.find('AUG')
        strand = strand[start_location:]
        list_strand = re.findall('...', strand)

    else:
        return 

    # print('list_strand is:', list_strand)
    # print('list codons_dict is:', codons_dict)

    amino_acids = [key for codon in list_strand for key, codon_list in codons_dict.items() if codon in codon_list]
    amino_acids_letters = [codon_list[0] for codon in list_strand for key, codon_list in codons_dict.items() if codon in codon_list]

    # recognize STOP codon in gene and remove proceeding codons
    if status == 'B' and 'STOP' in amino_acids:
        location_stop = amino_acids.index('STOP')
        amino_acids = amino_acids[:location_stop + 1]
    if status == 'B' and '*' in amino_acids_letters:
        location_stop = amino_acids_letters.index('*')
        amino_acids_letters = amino_acids_letters[:location_stop + 1]

    # display the amino acids
    print(f'Your amino acids are: {amino_acids}\nas letters: ', *amino_acids_letters, sep="")
    print()
    return


def test_seq(user_input):
    if re.fullmatch(rna_letters, user_input):
        rna_to_dna = user_input.maketrans('U','T')
        converted = user_input.translate(rna_to_dna)
        print('\nYour converted rna to dna strand is', converted, '\n')
        protein = protein_converter(user_input)
        return
    elif re.fullmatch(dna_letters, user_input):
        dna_to_rna = user_input.maketrans('T','U')
        converted = user_input.translate(dna_to_rna)
        print('\nYour converted dna to rna strand is', converted, '\n')
        protein = protein_converter(converted)
        return
    else:
        print('\nPlease enter a valid sequence\n')


while True:
    user_input = user_start()
    test = test_seq(user_input)