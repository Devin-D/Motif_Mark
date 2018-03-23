# Motif_Mark
A program to search intronic and exonic sequences for motifs and plot gene sequence introns and exons with color coded motifs  and postion information

Requires (python 3+)

Requires (pycairo)

How to use this tool:
- (Standard): User supplies Two input files,
 1) A gene sequence  file in fasta format with exonic regions capitalized (see example fasta header for header format)
 2) A motif file with 1 motif per line no spaces
 
    *example files can be found on github*
    
  specify -o flag  for the name of the output graph *will be in .svg format*
 
Example fasta header:

 ">INSR chr19:7149896-7151209 (reverse complement)"

Example of how to run this script from terminal:

./MotifMark.py sequence.fasta motifs.txt -o motif_plot 
