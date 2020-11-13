from matplotlib import pyplot as plt   
import forgi.visual.mplotlib as fvm  
import forgi.graph.bulge_graph as fgb
from time import time as t
import typing
import re
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import motifs
import seqlogo
from Bio.Alphabet import IUPAC

start = t()


#opening a dot or other FASTA-like file
lines = open("text.txt").readlines()


dot_brackets = []       #will be list of all loops in dot bracket notation
sequences = []          #will be list of all sequeces relating to loops
big_list_dots = [] 
big_list_seq = []
doc = ""
emp = ""
stem_lengths = []
big_list_stems = []
stem_count = 0 
GC_count = 0   #number of times a nucleotides is either guanine or cytosine
GC_list = [] 
RNA_names = [i.strip(">").strip("\n").strip(" ") for i in lines[::3]]
structures = [i.strip("\n").strip(" ") for i in lines[2::3]]
seqs = [i.strip() for i in lines[1::3]]


for dots,seq in zip(lines[2::3],lines[1::3]):
    for j,x in zip(dots,seq):

        if x == "G" or x == "C":
            GC_count += 1
            len_seq = len(seq)

        if j == "(":
            stem_count+=1


        elif j != "(":

            stem_lengths.append(stem_count) 
            stem_lengths = [i for i in stem_lengths if i > 0]
            stem_count = 0


        if j == ".":

            doc += j
            emp += x

        elif j != ".":

            dot_brackets.append(doc)
            dot_brackets = [ i for i in dot_brackets if len(i) >= 4 ]
            sequences.append(emp)
            sequences = [i for i in sequences if len(i) >= 4]
            
            doc = ""
            emp = ""

    big_list_dots.append(dot_brackets)
    dot_brackets = []

    big_list_seq.append(sequences)
    sequences = []

    big_list_stems.append(stem_lengths)
    stem_lengths =[]
    
    GC_list.append((GC_count/len_seq)*100)
    GC_count = 0   

#print(big_list_dots)
#print(big_list_seq)
    
def GC():

    """

    Calculates Guanine - Cytosine amount in a specific RNA sequence as a percentage
    and displays it as a histogram

    """

    plt.style.use('ggplot')

    names = [i for i, _ in enumerate(RNA_names)]

    plt.bar(names, GC_list, color='green')
    plt.xlabel("RNA names")
    plt.ylabel("GC Content (%)")

    plt.title("GC Contents of RNAs from Dot File")
    plt.xticks(names, RNA_names)
    plt.show()


A,U,G,C = 0,0,0,0       #counts of each nucleotide appearance
large = ""              #represents sequence of all loop sequences combined together
for i in big_list_seq:
    for j in i:
        large += j

for nucleotide in large:

    if nucleotide =="A":
        A+=1
    if nucleotide =="U":
        U+=1
    if nucleotide =="G":
        G+=1
    if nucleotide =="C":
        C+=1    


def nucleotide_distribution():       

    """
    Displays a pie chart of the distribution of nucleotides found in 
    RNA loop motifs

    """

    labels = "A","U","G","C"
    sizes = [A,U,G,C]
    explode = (0.1,0,0.1,0)
    fig1, ax1 = plt.subplots()
    ax1.pie(
        sizes,labels=labels, autopct='%1.1f%%',
                shadow=False,explode=explode, startangle=90
                )
    ax1.axis('equal')  
    plt.title("Frequency of Nucleotides in Loops of Length 4 or Greater")
    plt.show()


def draw_structure():         #Spits out a 2D image representing the structure of the RNA sequence

    for i in lines[1::3]:
        
        cg = str(i)

        fvm.plot_rna(
            cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
                backbone_kwargs={"linewidth":3}
                )

        plt.show()



def group_hairpins()->typing.List:          #classifies various hairpin loops by size and groups similar terms 

    for i in structures:

        bg = fgb.BulgeGraph.from_dotbracket(i)
        a = (bg.to_bg_string())
        print(a)
    
    

def base_pair_distribution():    #Distribution of certain nucleotide pairs (A-U)% vs. (A-G)% for example
    
    pass


def seq_freq_hairpins()->typing.Dict:


    n = int(input("Enter the length of nucleotide sequence you wish to see: "))
    new_seqs = []
    for i in seqs:
        for j in i:
            j = str(j)
            if j =="U":
                i = i.replace("U","T")
                new_seqs.append(i)
    
    new = []
    for i in structures:
        result = re.finditer(r'\(\.*\)',i)
        new.append(result)
    
    hairpin_coordinates = []
    for match in new:
        for q in match:
            hairpin_coordinates.append(q.span())
    #print(hairpin_coordinates)
    #print("Coordinates of hairpin sequences: " + str(hairpin_coordinates))
    newer=[]
    for sequence,hairpin_coordinate in zip(new_seqs,hairpin_coordinates):
        newer.append(sequence[hairpin_coordinate[0]:hairpin_coordinate[1]])
    #print("Sequences of hairpin loops: " + str(newer))
    
    frequency_map = {}

    for i in newer:
        for j in range(len(i)-n+1):

            window = i[j:j+n]
            if not window in frequency_map:
                frequency_map[window] =1

            else:
                frequency_map[window] +=1

    # print(frequency_map)

    #Gets 5 highest motifs 
    top_5 = sorted(frequency_map, key=frequency_map.get, reverse=True)[:5]
    
    #Gets 5 highest motif frequencies
    top_5_vals = [frequency_map[i] for i in top_5]

    plt.bar(top_5,top_5_vals)
    plt.title("Top 5 Most Frequent Motifs of length " + str(n))
    plt.show()

    hairs = []
    for key in frequency_map:
        for idx in range(frequency_map[key]):
            hairs.append(key)

    new_seqs = [x for x in hairs if not 'N' in x]
    new_seqs = [Seq(i) for i in new_seqs]

    m = motifs.create(new_seqs)

    #Converting the Position Frequency Matrix into a pandas DataFrame 
    pfm = pd.DataFrame(m.counts)
    print("PFM: ")
    print(pfm)

    #Converting a Position Frequency Matrix into a Position Weight Matrix (PWM)
    pwm = seqlogo.pfm2pwm(pfm)
    print("PWM: ")
    print(pwm)
    

    #Converting a PWM into a Position Probability Matrix (PPM)
    ppm = seqlogo.pwm2ppm(pwm)
    print("PPM: ")
    print(ppm)
    
    #seqlogo.seqlogo(ppm, ic_scale = False, format = 'svg', size = 'medium')


    return ("Consensus sequence: " + m.consensus)

if __name__ == "__main__":

    #GC()
    #nucleotide_distribution()
    #draw_structure()
    #print(RNA_names)
    #print(big_list_seq)
    #print(large)
    #group_hairpins()
    print(seq_freq_hairpins())
    print(str(t()-start) + " seconds")    
    #print(seqs)