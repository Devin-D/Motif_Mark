#!/usr/bin/env python3
"""A program to search intronic and exonic sequences for motifs and plot gene sequence introns and exons with color coded motifs 
and postion information

Requires (python 3+)
Requires (pycairo)

How to use this tool:
- (Standard): User supplies Two input files,
 1) A gene sequence  file in fasta format with exonic regions capitalized (see example fasta header for header format)
 2) A motif file with 1 motif per line no spaces
    *example files can be found on github*
 3) specify -o flag with string for the name of the output graph *will be in .svg format*
Example fasta header:
 >INSR chr19:7149896-7151209 (reverse complement)

Example of how to run this script from terminal:
./MotifMark.py sequence.fasta motifs.txt -o motif_plot 
    
"""
__author__ = 'Devin Dinwiddie'
__contact__ = "devin-d.github.io"
__version__ = '1.0'
__copyright__ = """Copyright (C) 2018 Devin Dinwiddie.
   
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>."""


class motif_marker(object):
    
    """ Function to filters a fasta file and a motif file grabbing sequence info and storing in dictionarys. as well 
    as get postion information for all motif locations on a sequence
    
    Attributes:
        
    """
    
    def __init__(self):
        
        self.self=self
        """ init
        
        arguments
        
        """
        
        
    def get_inex(line):
        ''' a function to split a fasta sequence into intron, exon, intron parts.
            exonic regions must be uppercase in fasta file
            aaaccttgAACCTTGccttgga --> aaaccttg,AACCTTG,ccttgga'''

        inexin = ["","",""] #list to store parts
        part = re.search('[A-Z]',line) #search fot the exon (uppercase)
        place = part.start() #mark the start of the exon
        inexin[0] = line[0:place] #1st intron was everything up to start of exon
        line=line[place:] #move to start of exon seq
        part = re.search('[a-z]',line) #find 2nd intron
        place = part.start() #mark start of 2nd intron
        inexin[1],inexin[2] = line[0:place], line[place:] # store exon and 2nd intron
        return inexin

    def conv_file_to_dict(file):
        '''Takes a fasta file and makes a dict with header line as key and list get_index as value
            {MBNL1500 chroX rev : [intron1,exon,intron2]'''
        with open(infasta)as fh:
            fasta_dict = {}
            head = ""
            seq = ""
            first = True #bool for recognizing first line

            while True:
                line = fh.readline().strip() #read lines + strip /n char
                if line =="": #eof
                    break

                if line[0] == ">": # if on a header line
                    # if  not on the first header line thenwe have stored the seq lines in seq (see else statement) 
                    #now split sequences into in,ex,in
                    #store in a dict with header as key
                    #make header the current line clear seq var <-- just this if on first header line
                    if not first: 
                        fasta_dict[head] = motif_marker.get_inex(seq)   
                    head = line
                    seq = ""
                else: 
                    seq += line #store seq info if not on header line
                    first = False

        fasta_dict[head] = motif_marker.get_inex(seq) #last seq info
        return fasta_dict


    def motif_char_find(motif_file):
        '''Takes a list of motifs and creates a dict of possible characters for each motif with motif as key
            and all IUPAC characters associated with each character (upper and lower) as value.
            {YGCY: [TtCcUu,Gg,Cc,TtCcUu]}'''

        iup = {'G':'[Gg]','C':'[Cc]','A':'[Aa]','T':'[TtUu]','U':'[TtUu]',
                'W':'[AaTtUu]','S':'[GgCc]','M':'[AaCc]','K':'[GgTtUu]','R':'[AaGg]','Y':'[CcTtUu]',
                'B':'[CcGgTtUu]','D':'[AaGgTtUu]','H':'[AaCcTtUu]','V':'[AaCcGg]',
                'N':'[GgCcAaTtUu]'}           
        motif_dict = {}
        with open(motifs_file)as fh:
            for line in fh:
                line = line.strip().upper()
                chars = []
                for c in line:
                    chars.append(iup[c]) # append list with chars from iup{}
                motif_dict[line] = chars #put in dict 
            #print(motif_dict)
        return motif_dict

    def motif_mark(seqs, motif_dict):
        '''Main function takes converted fasta file and motifs makes a list of dictionarys
            each dictionary contains intron, exon, intron with lengths {intron1:120}
            and all postions for each motif {MBNLChr19..._YGCY: 12,42,152,200}'''
        pos_list=[]
        for head in seqs:
            positions = {}
            #place intron and exon length in dict
            positions["intron1"] = len(seqs[head][0])
            positions["exon"] = len(seqs[head][1])
            positions["intron2"] = len(seqs[head][2])
            seq = "".join(seqs[head]) #join the sequences(values) to be able to interate

            for motif in motif_dict.keys():
                lookup = "".join(motif_dict[motif]) #make the motifs of interest (values in motif dict) in a list
                finder = re.finditer(lookup, seq) # #iterate motifs over the seq to Yeild matches 
                positions[head+"_"+motif]=[] # key to a list to hold postion info
                for pos in finder:
                    positions[head+"_"+motif].append(pos.start()) #place the postition info in the list {MBNL.._ygcy:[12,152,200]}

            pos_list.append(positions) # append each dict made to a list *postions dict is cleared w each key in fasta dict

        return pos_list
    
    
class motif_Draw(object):
    
    """ functions to draw gene with intron exon intron and draws motifs on the gene with postion information, legend and gene ID.
    Attributes:
       
    """
    
    def __init__(self):
        
        self.self=self
        """ init
        
        arguments
        
        """    
        
    def draw_inex(intron1 ,ex_len, intron2, context, start):
        '''Draws the intron as a skinny line based on intron lengths from pos_list
            exon is a fat line also based on length from pos_list'''
        context.set_line_width(1) #skinny line
        context.set_source_rgb(0,0,0) #blk

        context.move_to(start[0], start[1]) #(x,y) start
        context.line_to(start[0] + intron1 + ex_len + intron2, start[1]) #how long of a line to make
        context.stroke() # draw it


        context.set_line_width(10) #fat line (exon)
        context.move_to(start[0] + intron1, start[1]) #move to start of exon
        context.line_to(start[0] + intron1 + ex_len, start[1]) #how long of a line
        context.stroke() #draw it

    def draw_motifs(pos_list, context, start, motif):
        '''draws the position of the motifs on the introns/exon. Returns the color values for motif'''
        r = random.random() #generate three random floats for color
        g = random.random()
        b = random.random()
        global OE # for postion of motif location text (odd even...up down)
        context.set_line_width(10) #line size
        context.set_source_rgb(r, g, b) #color


        for p in pos_list:
            #print(pos)
            context.move_to(start[0] + p, start[1]) #move to postion of motif
            context.line_to(start[0] + p + len(motif), start[1]) # line it
            context.stroke() #draw it
            if OE%2==0: #if even draw postion on top
                context.move_to(start[0] + p, start[1] - 20)
                context.show_text(str(p))

            else: #if odd draw postion on bottom *not perfect some overlap
                context.move_to(start[0] + p, start[1] + 20)
                context.show_text(str(p))
            OE+=1 #increase counter
            
        return [r, g, b]

    def draw_legend(r, g, b, start, lpos, motif):
        '''draws the legend  with each motif in color'''
        context.set_source_rgb(r, g, b) #color
        context.move_to(start[0] - 50, start[1] + lpos) # where to put legend
        context.show_text(motif) #draw the text (motif)

    def draw_gene(start, name):
        '''draws the name of the gene of interest next to image'''
        context.set_source_rgb(0,0,0) #color
        context.move_to(start[0] - 100, start[1]) #where to put gene
        context.show_text(name)  #draw gene name
    
    def draw_draw(postionList):
        global start
        for z in out:
            #print(z["intron1"])
            motif_Draw.draw_inex(z["intron1"], z["exon"], z["intron2"], context, start) #draw the intron exon intron

            lpos = -20 #postion of legend 
            for x in z.items():

                if x[0][0] == ">":
                    #get motif 
                    motif = x[0].split("_")[1]
                    
                    #get gene name 
                    gene = x[0].split(" ")[0][1:]
                    
                    if x[1]!=[]:
                        color = motif_Draw.draw_motifs(x[1], context, start, motif) #draw motif
                        motif_Draw.draw_gene(start, gene) #put in gene name
                        motif_Draw.draw_legend(color[0], color[1], color[2], start, lpos, motif) #put in legend (motif name colored same as motif)
                        lpos = lpos + 10 #move down 

            start = [start[0], start[1] + 100] # adjust stat pos
            
def file_check(parser, arg):
    """ checks if input files exist"""
    if not os.path.exists(arg):
        parser.error("The file {0} does not exist!".format(arg))
    else:
        return str(arg)
            
if __name__ == '__main__':
    
    import argparse,re,cairo,random,os
    

    parser = argparse.ArgumentParser(description=__doc__.format(author=__author__, contact=__contact__), formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

    pgroup = parser.add_argument_group("Input Files")
    pgroup.add_argument('Fasta', metavar='IN.fasta', type=lambda x:file_check(parser, x), help='input fasta file with lowercase introns upper case Exons')
    pgroup.add_argument('Motifs', metavar='IN.motif', type=lambda x:file_check(parser, x), help='input motif file with with 1 motif per line no spaces')
    
    ogroup = parser.add_argument_group("Options")
    ogroup.add_argument('-o','--out', dest='out_prefix', help='prefix of output file for graph will be in .svg format',type=str)
    ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
    ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')
    
    args = parser.parse_args()
    #set up the drawing surface
    surface = cairo.SVGSurface(args.out_prefix+".svg", 1600, 1000)
    context = cairo.Context(surface)
    context.set_line_width(1)
    start = [200,100]
    OE=0 # counter for motif postion text 

    #draw introns/exons and their associated motif positions

    infasta=args.Fasta
    motifs_file=args.Motifs
    #seqs=motif_marker.conv_file_to_dict(infasta)
    # motifs=motif_char_find(motifs_file)
    out=motif_marker.motif_mark(motif_marker.conv_file_to_dict(infasta),motif_marker.motif_char_find(motifs_file))
    #print(out)
    motif_Draw.draw_draw(out)
