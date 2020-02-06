# Sequence_Alignment
A sequence aligner, that given: two FASTA files, type of alignment, and a scoring matrix, will print the optimal alignment and its score. 
The algorithm is using dynamic programming.  

FASTA is a common format for biological sequence:  
">name of seq1  
ACACGGTGGACCGGAT  
AACACGGTAATACCAG"  

## Input

### FASTA files
For simplicity we will assume the aligner should handle only FASTA sequences from the nucleotide alphabet , e.g = A; T; G;C.  

### Score matrix
The first row and first column describe the characters in the alphabet or gap (-), and the rest of the cells describe the score for substitution or deletion.  

## Alignments
Given a pair of input sequences X and Y, the aligner should support the following types of alignments:
### Global alignment ("global")
In this type of alignment, the algorithm seeks for the best match between X and Y , such that all the characters of X and Y are aligned either one to another , or to a gap.
### Local alignment ("local")
In this type of alignment, the algorithm seeks for the best substring match of X and Y.  

## Example
The program prints the optimal alignment and bellow the alignment "alignment_type: score".  

TCGAATCG--CACGCGCGGCTCTCCTTAGAACCGGCCGGCTCCCGAATAA  
TTGGGTCGGTTTCACCCGGTCTTCATCCG--CCGACTGTTTAAAAACCAA  

TGTTTCAGTGTTTGACAAACTCAATCGGAGGTCTCGGAAGA---AGTATC  
CAAGGTAAGAGGAGGGGAGCTTTGTTGTTGTTTTAACGTGTGTTAGTGAC  

AAAAAAAAAAAA  
AAAAAAAAAAAA  
global:20  
