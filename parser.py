from matplotlib import pyplot as plt
from collections import Counter
import Bio   #Biopython
import forgi.visual.mplotlib as fvm  
from time import time as t
import statistics as stats
start = t()

#opening a dot file
lines = open("file.txt").readlines()

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
    
def GC():
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
    # Pie chart of nucleotide frequency 
    labels = "A","U","G","C"
    sizes = [A,U,G,C]
    explode = (0.1,0,0.1,0)
    fig1, ax1 = plt.subplots()
    ax1.pie(sizes,labels=labels, autopct='%1.1f%%',
                shadow=False,explode=explode, startangle=90)
    ax1.axis('equal')  
    plt.title("Frequency of Nucleotides in Loops of Length 4 or Greater")
    plt.show()

def draw_structure():
    pass

   
if __name__ == "__main__":
    GC()
    nucleotide_distribution()
    #print(RNA_names)
    #print(big_list_seq)
    #print(large)
    print(str(t()-start) + " seconds")    