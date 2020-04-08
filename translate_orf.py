#! /usr/bin/env/ python3
# the goal of this script is to search for an open reading frame and
#translate that open reading frame. The open reading frame that gets
#translated will be the first one that is read by the script.

import sys
import re

#pulling the first ORF
#first_orf = r'('+r'|'.join(starts)+r')([AUCG]{3})*('+r'|'.join(stops)+r')'
#protein = []
#translating the first ORF
#for i in range(0, len(first_orf), 3):
#    codon = first_orf[i: i+3]
#    amino_acid = genetic_code.get(codon, '*')
#    if codon == "AUG":
#        aa_seq = translate_sequence(
#            first_orf = first_orf[i:],
#            genetic_code = genetic_code)
#        if aa_seq:
#            protein.append(aa_seq)
#return protein

def main():
    import argparse
    import find_orf
    import translate

    # Create a command-line parser object
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codons',
            type = str,
            nargs = '+', #one or more arguments
            default = ['UAA', 'UAG', 'UGA'],
            help = ('One or more possible stop codons.'))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    #Check to see if the path option was set to True by the caller. If so, parse
    #The sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codons,
            stop_codons = args.stop_codons)
    sys.stdout.write('{}\n'.format(orf))

genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}

translation = translate.translate_sequence(rna_sequence = orf, genetic_code = genetic_code)
sys.stdout.write('{}\n'.format(translation))

if __name__ == '__main__':
    main()

