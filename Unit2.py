# This code was implemented for Genome Sequencing (Bioinformatics II) 


####Helper Methods from Previous Unit

# This method prints a list, each element in the list is printed in a new line
def PrintList(output):
    for i in range(len(output)):
        print(output[i])

# This method combines all elements in a list into a single string
# and prints the string on a single line
def PrintListFlat(list):
    out=''
    for  i in range(len(list)):
        out+= str(list[i])+" "
    print(out)

# This method takes a DNA string that is deliminated by
# a space and stores each DNA sequence as a separate element
# in a list
# Input: DNA string that has sequences separated by " "
# Output: A list that contains each sequence
def DNAStringtoList(DNA_string):
    curr_seq=''
    DNA_list=[]
    for i in range(len(DNA_string)):
        if DNA_string[i]!=" ":
            curr_seq+=DNA_string[i]
        else:
            DNA_list.append(curr_seq)
            curr_seq=''
    DNA_list.append(curr_seq)
    return DNA_list

#This method calculates the hamming distance between two strings
#Input: text1 and text2 are two strings
#Output: distance which is an int that counts the number of instances the two stings
#do not match
def Hamming(text1,text2):
    n=len(text1)
    l=len(text2)
    text1 = text1.upper()
    text2 = text2.upper()
    distance = 0

    if n>l:
        for i in range(l):
            if text1[i]!=text2[i]:
                distance+=1
        distance+=n-l
    else:
        for i in range(n):
            if text1[i]!=text2[i]:
                distance+=1
        distance+=l-n
    #print(distance)
    return distance

#This method finds the reverse complement of a DNA sequence
#Input: DNA is a string containing only A,C,G,or T
#Output: A string that is the reverse complement of DNA
def RevComp(DNA):
    DNA=DNA.upper()
    cDNA=""
    for i in range(len(DNA)):
        if DNA[i]=="A":
            cDNA+="T"
        elif DNA[i]=="C":
            cDNA+="G"
        elif DNA[i]=="G":
            cDNA+="C"
        elif DNA[i]=="T":
            cDNA+="A"
    #print(cDNA[::-1])
    return cDNA[::-1]
###########################################################

import random
### Unit 2
### Week 1

# Input: A string Text and a positive int k
# Output: A list of all substrings of length k 
# from the text
def Composition(Text,k):
    comp=[]
    for i in range(len(Text)-k+1):
        comp.append(Text[i:i+k])
    #comp.sort()
    return comp

# This method stitches together a path list to complete
# a path. This is used to stitch together a list of kmers
# to reconstruct a genome
# Input: A list of strings pathList
# Output: The 
def PathToGenome(pathList):
    genome=""
    #pathList=DNAStringtoList(path)
    n= len(pathList)
    k=len(pathList[0])
    for i in range(n):
        if i == 0:
            genome+=pathList[i][0:k]
        else:
            genome+=pathList[i][-1]
    return genome

# This method returns the prefix defined as substring of a 
# string consisting of the entire string except the last element
# Input: A string pattern
# Output: A prefix substring
def Prefix(pattern):
    k=len(pattern)
    return pattern[0:k-1]

# This method returns the suffic defined as substring of a 
# string consisting of the entire string except the first element
# Input: A string pattern
# Output: A suffix substring
def Suffix(pattern):
    return pattern[1::]

# This method takes a string of DNA sequences
# Input: A string of DNA sequences delimited by " "
# Output: A dictionary graph with keys that are DNA 
# sequences from patterns and the corresponding values are 
# sequences that overlap with the key such that the
# Suffix(key)==Prefix(pattern)
def OverlapGraph(patterns):
    graph={}
    patternsList=DNAStringtoList(patterns)

    #Adds all patterns into a dict
    for pattern in patternsList:
        if pattern not in graph:
            graph[pattern]=[]

    #Adds the nodes
    for keys in graph:
        for pattern in patternsList:
            if Prefix(pattern)==Suffix(keys):
                graph[keys].append(pattern)

    # Removes keys with no entries
    keysToRemove=[]
    for keys in graph:
        if len(graph[keys])<1:
            keysToRemove.append(keys)
    for toRemove in keysToRemove:
        graph.pop(toRemove)
    return graph

# This method takes a dictionary graph and prints each
# key and its corresponding values one key per line
def PrintGraph(graph):
    keys = list(graph.keys())
    keys.sort()
    for key  in keys:
        out=''
        for i in range(len(graph[key])):
            out+=" "+str(graph[key][i])
        print(str(key)+":"+out)

# This method generates combinations of either 
# 0's or 1's with a mismatch d
# Input: pattern is a string consisting of either 0's or 1's
# and d is a positve int
# Output: A list of all combinations of 0's and 1's with at most
# d mismatches
def BinaryNeightbors(pattern,d):
    #Base cases of the recursive alg
    if d == 0:
        return pattern
    if len(pattern)==1:
        return ['0','1']

    Neightborhood=[]
    SuffixNeighbors=BinaryNeightbors(pattern[1::],d)
    for i in range(len(SuffixNeighbors)):
        text=SuffixNeighbors[i]
        if Hamming(text,pattern[1::]) < d:
            for binary in ['0','1']:
                Neightborhood.append(binary+text)
        else:
            Neightborhood.append(pattern[0]+text)
    return Neightborhood

# This method takes a string text and breaks it into 
# edges of length k and stores the corresponding nodes
# Input: A positive int k and a string text
# Output: A list of lists nodes and edges
def PathGraph(k,text):
    edge=[]
    node=[]
    for i in range(len(text)-k+1):
        edge.append(text[i:i+k])
        node.append(text[i:i+k-1])
    return [node,edge]

# This method generates a DeBruijn Graph
# from a complete genome
# Input: A string text and a positive int k
# Output: A dictionary with a DeBruijn graph
# with fragments of length k
def DeBruijnGraph(k,text):
    graph={}
    path=PathGraph(k,text)

    node=path[0]
    edge=path[1]

    for i in range(len(node)):
        if node[i] not in graph:
            graph[node[i]]=[]
            graph[node[i]].append(Suffix(edge[i]))
        else:
            graph[node[i]].append(Suffix(edge[i]))

    return graph

# This method generates a DeBruijn Graph
# from a list of strings pattern
# Input: A list of strings patterns and a positive int k
# Output: A dictionary with a DeBruijn graph
# with fragments of length k
def DeBruijn(patterns):
    #patternList=DNAStringtoList(patterns)
    graph={}
    for pattern in patterns:
        pre=Prefix(pattern)
        suf=Suffix(pattern)
        if pre not in graph:
            graph[pre]=[]
            graph[pre].append(suf)
            graph[pre].sort()
        else:
            graph[pre].append((suf))
            graph[pre].sort()
    return graph

##Week 2

# This is a helper method to read text files
# and store the output as a dictionary
def dictFromText(filename, dest = 0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    dict = {}
    for line in f:
        main_node,adj_nodes=line.rsplit(': ')
        adj_nodesList=adj_nodes.rsplit(' ')
        dict[main_node]=[]
        for adjNode in adj_nodesList:
            dict[main_node].append(adjNode.strip())
    return dict

# This method solves a EulerianCycle in a given graph
# Input: A dictionary graph and an optional argument
# for a set starting node. Default is to choose a random
# starting node
# Output: A list cycle that contains the order of nodes
# that complete a EulerianCycle from the given graph
def EulerianCycle(graph,curr_node=0):
    stack=[]
    cycle=[]

    if curr_node==0:
        keysList = list(graph.keys())
        n = len(keysList)
        dice = random.randint(0, n - 1)
        curr_node=keysList[dice]

    while len(list(graph.values()))>0:
        if len(graph[curr_node]) == 0:
            cycle.append(curr_node)
            if len(stack)>0:
                curr_node=stack[-1]
                stack.pop(-1)
            else:
                return cycle[::-1]
        else:
            stack.append(curr_node)
            out_degrees=graph[curr_node]
            if len(out_degrees)>1:
                dice=random.randint(0,len(out_degrees)-1)
            else:
                dice=0
            temp_node=graph[curr_node][dice]
            graph[curr_node].pop(dice)
            graph={k:v for k,v in graph.items() if v is not None}
            curr_node=temp_node[:]
    return cycle[::-1]


# This method solves a Eulerian Path in a given graph
# Input: A dictionary graph 
# Output: A list path that contains the order of nodes
# that complete a Eulerian Path
import copy
def EulerianPath(graph):
    path=[]
    #Get all Nodes
    allNodes={}
    edge_List=[]
    for node in graph:
        outbound=graph[node]
        allNodes[node]=0
        for out_node in outbound:
            allNodes[out_node]=0
            edge_List.append(out_node)
    allNodes=list(allNodes.keys())
    start_node=[]
    end_node=[]


    #Searches for Start and Ending Nodes
    for node in allNodes:
        count=0
        for edge in edge_List:
            if edge==node:
                count+=1

        if node not in graph:
            end_node=node
        elif (-count + len(graph[node])==1):
            start_node=node
        elif (count-len(graph[node]))==1:
            end_node=node

    #print('start node is: '+str(start_node))
    #print('end node is: ' +str(end_node))
    #add the fake edge
    if end_node not in graph:
        graph[end_node]=[start_node]
    else:
        graph[end_node].append(start_node)

    temp_graph=copy.deepcopy(graph)
    path=EulerianCycle(temp_graph)
    while (path[-2] != end_node):
        temp_graph=copy.deepcopy(graph)
        path=EulerianCycle(temp_graph,start_node)
    path.pop(-1)
    return path

# This method reconstructs a genome from a collection 
# of DNA strings
# Input: A list of strings patterns of equal length
# Output: A string that represents a reconstructed genome
def StringReconstruction(patterns):
    # Generates a DeBruijn graph from the given patterns
    dB=DeBruijn(patterns)
    # PrintGraph(dB)
    print("de Bruijn graph is done")
    print('Searching for Path...')
    # Solving a Eulerian Path through the graph 
    path=EulerianPath(dB)
    print('Path Done')
    # Concatenates the path to form a reconstructed genome
    genome=PathToGenome(path)
    return genome

# This method generates a universal binary circular string
# consisting of fragments of length k
# Input: A positive int k
# Output: A string out that contains the sequence 
def UniversalCircularString(k):
    first_kmer=""
    for i in range(k):
        first_kmer+='0'
    all_binary_kmers=BinaryNeightbors(first_kmer,k)
    db=DeBruijn(all_binary_kmers)
    path=EulerianCycle(db)
    path.pop(-1)
    genome=PathToGenome(path)
    out=''
    for i in range(len(genome)-(k-1)+1):
        out+=genome[i]
    return out

# This method breaks a genome into paired reads
# Input: A string genome, a positive int k, and 
# a positive int d that is the distance between
# the read pairs
# Output: a list of paired reads from the genome
def PairedComposition(genome,k,d):
    pairedRead=[]
    for i in range(len(genome)-2*k-d+1):
        pattern1=genome[i:i+k]
        pattern2=genome[i+k+d:i+d+k+k]
        out.append([pattern1,pattern2])
    pairedRead.sort()
    outString=""
    for i in range(len(pairedRead)):
        outString+="("+str(pairedRead[i][0])+"|"+str(pairedRead[i][1])+") "
    print(outString)
    return pairedRead

# This method stores paired reads from a text file and
# stores them in a list
def PairedReadsFromText(filename, dest = 0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    out=[]
    for line in f:
        allPairs=line.rsplit(' ')
        for pair in allPairs:
            pattern1,pattern2=pair.split("|")
            out.append([pattern1.strip(),pattern2.strip()])
    return out

# This method returns the prefix from a paried read
# Input: A list with two elements that are a paired read
# Output: A string with the prefix of each pattern from read
# separated by a "|"
def PairedPrefix(read):
    temp_pattern1=read[0]
    temp_pattern2=read[1]
    k=len(temp_pattern1)
    pattern1=''
    pattern2=''
    for i in range(k-1):
        pattern1+=temp_pattern1[i]
        pattern2+=temp_pattern2[i]
    return str(pattern1)+"|"+str(pattern2)

# This method returns the suffix from a paried read
# Input: A list with two elements that are a paired read
# Output: A string with the suffix of each pattern from read
# separated by a "|"
def PairedSuffix(read):
    pattern1=read[0]
    pattern2=read[1]

    return str(pattern1[1::])+"|"+str(pattern2[1::])

# This method generates a DeBruijn graph from paired reads
# Input: A list of strings patterns of equal length
# Output: A dictionary graph that contains a DeBruijn graph
# from the paired reads
def DeBruijnPairedReads(readsList):
    #patternList=DNAStringtoList(patterns)
    graph={}
    for read in readsList:
        #print("read is: "+str(read))
        pre=PairedPrefix(read)
        suf=PairedSuffix(read)
        #print("Prefix is: "+ str(pre))
        #print("Suffix is: "+str(suf))
        if pre not in graph:
            graph[pre]=[]
            graph[pre].append(suf)
            graph[pre].sort()
        else:
            graph[pre].append((suf))
            graph[pre].sort()
    return graph

# This method returns a genome from a collection of paired reads
# Input: A list of an ordered set of paired reads GappedPatterns, a 
# positive int k that is the length of a read, and d a positive int 
# that separates the paired reads
# Output: A string that represents the reconstructed string from the paired 
# reads. Returns none if no string exists
def StringSpelledByGappedPatterns(GappedPatterns,k,d):
    firstPatterns=[]
    secondPatterns=[]

    for i in range(len(GappedPatterns)):
        first, second = GappedPatterns[i].split('|')
        firstPatterns.append(first)
        secondPatterns.append(second)

    prefixString=PathToGenome(firstPatterns)
    suffixString=PathToGenome(secondPatterns)

    for i in range(k+d,len(prefixString)):
        j=i-k-d
        if prefixString[i]!=suffixString[j]:
            print("There is no string spelled by the gapped patterns")
            return None
    return prefixString+suffixString[-1*(k+d):]

# This method returns a genome from a collection of paired reads
# Input: A list of a set of paired reads patterns, a positive int 
# k that is the length of a read, and d a positive int that 
# separates the paired reads
# Output: A string that represents the reconstructed string from the paired 
# reads. Returns none if no string exists
def PairedStringReconstruction(k,d,patterns):
    dB=DeBruijnPairedReads(patterns)
    #print(dB)
    print("Graph Done")
    print('Searching for Path...')
    path=EulerianPath(dB)
    print('Path Done')
    print(path)
    genome=StringSpelledByGappedPatterns(path,k,d)
    return genome

# This method takes a dictionary graph and returns a maximal
# non-branching path from the graph
# Input: A dictionary graph
# Output: A list Paths that contain all non-branching paths
def MaximalNonBranchingPaths(graph):
    Paths=[]
    
    # Get all possible nodes in the graph
    allNodes = {}
    edge_List = []
    for node in graph:
        outbound = graph[node]
        allNodes[node] = 0
        for out_node in outbound:
            allNodes[out_node] = 0
            edge_List.append(out_node)
    allNodes = list(allNodes.keys())

    # Determines the in and out degrees
    # outputs coundDict, each key is a node in the graph and the value is a list. the first value is in(v) and the
    # second value is out(v)

    countDict={}
    for node in allNodes:
        in_count = 0
        for edge in edge_List:
            if edge == node:
                in_count += 1
        if node in graph:
            countDict[node]=[in_count,len(graph[node])]
        else:
            countDict[node]=[in_count,0]

    for node in countDict:
        in_v= countDict[node][0]
        out_v = countDict[node][1]
        if in_v!=1 or out_v!=1:
            if out_v>0:
                for edge in graph[node]:
                    edge_in=countDict[str(edge)][0]
                    edge_out=countDict[str(edge)][1]
                    NonBranch=[node]

                    while edge_in==1 and edge_out==1:
                        NonBranch.append(edge)
                        edge=graph[edge][0]
                        edge_in = countDict[str(edge)][0]
                        edge_out = countDict[str(edge)][1]
                    NonBranch.append(edge)
                    Paths.append(NonBranch)

    #find cycles
    cycleList=[]
    for node in graph:
        cycle=[node]
        curr_node=str(graph[node][0])
        while str(curr_node)!=str(node):
            cycleKeys = dict.fromkeys(cycle)

            if curr_node in cycleKeys:
                break
            cycle.append(curr_node)
            if curr_node not in graph:
                break
            elif len(graph[curr_node])>0:
                curr_node=str(graph[curr_node][0])
            else:
                break
        if curr_node==node:
            cycle.append(curr_node)
            cycleList.append(cycle)

    # Remove duplicate cycles

    to_remove=[]
    for i in range(len(cycleList)-1):
        keys=dict.fromkeys(cycleList[i])
        keys=list(keys.keys())
        keys.sort()

        for j in range(i+1,len(cycleList)):
            other_keys=dict.fromkeys(cycleList[j])
            other_keys=list(other_keys.keys())
            other_keys.sort()
            remove= True
            if len(keys)!=len(other_keys):
                remove=False
            else:
                for k in range(len(keys)):
                    if other_keys[k]!=keys[k]:
                        remove=False
            if remove:
                to_remove.append(i)
                to_remove=dict.fromkeys(to_remove)
                to_remove=list(to_remove.keys())

    to_remove=to_remove[::-1]

    for key in to_remove:
        cycleList.pop(key)

    to_remove=[]
    #Remove if duplicate in path
    for i in range(len(cycleList)):
        cycleNodes=dict.fromkeys(cycleList[i])
        cycleNodes=list(cycleNodes.keys())
        cycleNodes.sort()
        for path in Paths:
            pathNodes=dict.fromkeys(path)
            pathNodes=list(pathNodes.keys())
            pathNodes.sort()
            remove=False
            if len(cycleNodes)==len(pathNodes):
                remove=True
                for j in range(len(cycleNodes)):
                    if cycleNodes[j]!=pathNodes[j]:
                        remove=False
                        break
            if remove:
                to_remove.append(i)
                to_remove = dict.fromkeys(to_remove)
                to_remove = list(to_remove.keys())


    to_remove = to_remove[::-1]
    for key in to_remove:
        cycleList.pop(key)

    for key in cycleList:
        Paths.append(key)
    # To print
    #for path in Paths:
    #    PrintListFlat(path)
    return Paths

# This method generates contigs from a collection of DNA strings
# Input: A list of strings patterns
# Outputs: A list of strings of the reconstructed contigs
def Contig(patterns):
    dB=DeBruijn(patterns)
    maxPaths=MaximalNonBranchingPaths(dB)
    contigs=[]
    for path in maxPaths:
        contigs.append(PathToGenome(path))
    contigs.sort()
    return contigs

def Unit2Week2InputReader(filename,problem):
    location = "/Users/nathanng/PycharmProjects/bioinfo/"
    file = location + filename
    f = open(file, "r")
    out=[]
    for line in f:
        if problem == 0:
            out.append(line.strip())
        elif problem==1:
            pattern1,pattern2=line.split("|")
            pattern1=pattern1.strip()
            pattern2=pattern2.strip()
            pattern1=pattern1[1:]
            pattern2=pattern2[:-1]
            out.append([pattern1.strip(),pattern2.strip()])
    return out

#Week 3

# This method translates a string RNA into an amino acid
# Input: A string RNA
# Output: A string seq that is the peptide sequence
def Translation(RNA):
    codon={
        "UUU":"F",
        "UUC":"F",
        "UUA":"L",
        "UUG": "L",
        "CUG": "L",
        "AUG": 'M',
        "GUG": 'V',
        'UCU': 'S',
        'CCU': 'P',
        'ACU': 'T',
        'GCU': 'A',
        'UCC': 'S',
        'CCC':'P',
        'ACC': 'T',
        'GCC': 'A',
        'UCA': 'S',
        'CCA': 'P',
        "ACA": "T",
        'GCA': 'A',
        'UCG':'S' ,
        'CCG': 'P',
        'ACG': 'T',
        "GCG": 'A',
        'UAU': 'Y',      
        'CAU': 'H',     
        'AAU': 'N',      
        'GAU': 'D',
        'UAC': 'Y',      
        'CAC': 'H',      
        'AAC': 'N',      
        'GAC': 'D',
        'UAA': 'Stop',
        'CAA':'Q',
        "AAA": 'K',      
        "GAA": 'E',
        'UAG':'Stop',   
        'CAG': 'Q',      
        'AAG': 'K',      
        'GAG': 'E',
        'UGU': 'C',  
        "CGU": 'R', 
        "AGU": 'S',      
        'GGU': 'G',
        'UGC': 'C',      
        "CGC": 'R',      
        "AGC": 'S',      
        'GGC': 'G',
        'UGA':'Stop',  
        'CGA':'R',      
        'AGA':'R',      
        "GGA" :'G',
        'UGG': 'W' ,     
        'CGG': 'R',      
        'AGG': 'R',      
        'GGG': 'G',
        'CUU': 'L',      
        'AUU': 'I',     
        'GUU': 'V', 
        "CUC": 'L', 
        'AUC': 'I',      
        'GUC': 'V',
        'CUA': 'L',      
        'AUA': 'I',      
        'GUA': 'V'
    }
    seq=''
    for i in range(int(len(RNA)/3)):
        start = 3*i
        currcodon=RNA[start:start+3]
        if codon[currcodon]=='Stop':
            #print(seq)
            return seq
        seq+=codon[currcodon]
    #print(seq)
    return seq

# This method transcribes DNA into an RNA string
# Input: A string of DNA
# Output: A string of RNA
def DNA2RNA(dna):
    #dna=dna.upper()
    rna=''
    for i in range(len(dna)):
        if dna[i]== "T":
            rna+="U"
        else:
            rna+=dna[i]
    return rna

# This method determines if a peptide is found in a DNA 
# string
# Input: A string of DNA text and a string peptide
# Output: A list of the locations indicating where the
# peptide can be found in the string Text
def PeptideEncoding(Text,peptide):
    k=len(peptide)*3
    rna=DNA2RNA(Text)
    peptide_location=[]

    for i in range(len(rna)-1-k):
        dna_seq=Text[i:i+k]
        rev=RevComp(dna_seq)
        rna_seq=DNA2RNA(dna_seq)
        rev_rna_seq=DNA2RNA(rev)
        currpeptide=Translation(rna_seq)
        rev_curr_peptide=Translation(rev_rna_seq)
        if currpeptide==peptide:
            peptide_location.append(dna_seq)
        if rev_curr_peptide==peptide:
            peptide_location.append(dna_seq)
    return peptide_location

# This method reads a genome from a text file
# and returns a string
def GenomeFromText(filename,dest=0):
    genome=""
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"

    file = location + filename
    f = open(file, "r")
    for line in f:
        genome+=line.strip()
    return genome

# This method calculates the approximate mass
# of a peptide
# Input: A string peptide
# Output: An int mass
def PeptideMass(peptide):
    mass=0
    massTable={
        'G':57,
        'A':71,
        'S':87,
        'P':97,
        'V':99,
        'T':101,
        'C':103,
        'I':113,
        'L':113,
        'N':114,
        'D':115,
        'K':128,
        'Q':128,
        'E':129,
        'M':131,
        'H':137,
        'F':147,
        'R':156,
        'Y':163,
        'W':186
    }
    if peptide=="":
        return 0

    for i in range(len(peptide)):
        mass+=massTable[peptide[i]]
    return mass

# This method takes an int mass and returns
# the peptide of that mass. For masses 113 and 128
# only I and K are given
# Input: An int mass
# Output: A string corresponding to the amino acid with the given mass

def MassToPeptide(mass):
    if mass == 57:
        return "G"
    elif mass== 71:
        return 'A'
    elif mass == 87:
        return 'S'
    elif mass == 97:
        return 'P'
    elif mass == 99:
        return 'V'
    elif mass == 101:
        return 'T'
    elif mass == 103:
        return 'C'
    elif mass == 113:
        return 'I'
    elif mass == 114:
        return 'N'
    elif mass == 115:
        return 'D'
    elif mass == 128:
        return 'K'
    elif mass == 129:
        return 'E'
    elif mass == 131:
        return 'M'
    elif mass == 137:
        return 'H'
    elif mass == 147:
        return 'F'
    elif mass == 156:
        return 'R'
    elif mass == 163:
        return 'Y'
    elif mass == 186:
        return 'W'
    else:
        return ""

# This method takes a string peptide and returns
# a string comprising the masses of each amino acid
# in the peptide separated by "-"
# Input: A string peptide
# Output: A string of amino acid masses separated by
# "-"
def PeptideToMass(peptide):
    mass=""
    massTable={
        'G':57,
        'A':71,
        'S':87,
        'P':97,
        'V':99,
        'T':101,
        'C':103,
        'I':113,
        'L':113,
        'N':114,
        'D':115,
        'K':128,
        'Q':128,
        'E':129,
        'M':131,
        'H':137,
        'F':147,
        'R':156,
        'Y':163,
        'W':186
    }
    if peptide=="":
        return 0

    for i in range(len(peptide)):
        mass+=str(massTable[peptide[i]])+"-"
    return mass[:-1]

# This method generates the linear mass spectrum of a peptide
# Input: A string peptide
# Output: A list spectrum that contains the mass spectrum of the peptide
def LinearSpectrum(Peptide):
    PrefixMass=[0]
    massTable = {
        'G': 57,
        'A': 71,
        'S': 87,
        'P': 97,
        'V': 99,
        'T': 101,
        'C': 103,
        'I': 113,
        'L': 113,
        'N': 114,
        'D': 115,
        'K': 128,
        'Q': 128,
        'E': 129,
        'M': 131,
        'H': 137,
        'F': 147,
        'R': 156,
        'Y': 163,
        'W': 186
    }
    for i in range(1,len(Peptide)+1):
        for key in massTable:
            if key == Peptide[i-1]:
                PrefixMass.append(PrefixMass[i-1]+massTable[key])
    spectrum=[0]
    for i in range(len(Peptide)):
        for j in range(i+1,len(Peptide)+1):
            spectrum.append(PrefixMass[j]-PrefixMass[i])
    spectrum.sort()
    return spectrum

# This method generates the cyclic mass spectrum of a peptide
# Input: A string peptide
# Output: A list spectrum that contains the mass spectrum of the peptide
def CyclicSpectrum(Peptide):
    spectrum=LinearSpectrum(Peptide[:-1])
    TotalMass=PeptideMass(Peptide)
    diff=[]
    for mass in spectrum:
        diff.append(TotalMass-mass)
    for value in diff:
        spectrum.append(value)
    spectrum.sort()
    return spectrum

# This method calculates the number of subpeptides
# that can be generated from a peptide of length n
# Input: An int n
# Output: A int sum
def Subpeptides(n):
    sum=0
    for i in range(n):
        sum+=n-i
    return sum+1

# This method expands a set of peptides by every possible
# amino acid
# Input: A set of strings that represent peptides
# Outout: A set of strings containing every original string
# expanded by a single amino acid
def PeptideExpander(peptideSet):
    expanded=set()
    for peptide in peptideSet:
        for ammino_acid in ['G', 'A','S','P','V','T','C','I','N','D','K','E','M','H','F','R','Y','W']:
            expanded.add(str(peptide)+str(ammino_acid))
    return expanded

# This method expands a set of peptides by the amino acids found
# in a given spectrum
# Input: A set of strings that represent peptides
# Outout: A set of strings containing every original string
# expanded by a single amino acid
def LimitedPeptideExpander(peptideSet,convoluted_spec):
    limitedAminoAcids=[]
    for mass in convoluted_spec:
        peptide=MassToPeptide(mass)
        if peptide!="":
            limitedAminoAcids.append(peptide)
    expanded=set()
    for peptide in peptideSet:
        for ammino_acid in limitedAminoAcids:
            expanded.add(str(peptide)+str(ammino_acid))
    return expanded


# This method is a branch and bound algorithm
# Inout: A list that represents the mass spectrum of a peptide
# Output: A set of strings containing the masses of a linear 
# peptide sequence that match the spectrum
import copy
def CyclopeptideSequencing(Spectrum):
    CandidatePeptide={''}
    FinalPeptide=set()
    while len(CandidatePeptide)>0:
        CandidatePeptide=PeptideExpander(CandidatePeptide)
        CandidatePeptide.discard("")
        temp=copy.deepcopy(CandidatePeptide)
        for peptide in CandidatePeptide:
            print()
            print('peptide is: '+str(peptide))
            if PeptideMass(peptide) == Spectrum[-1]:
                if CyclicSpectrum(peptide) == Spectrum:
                    print('Final Peptide Found!')
                    FinalPeptide.add(PeptideToMass(peptide))
                    print(FinalPeptide)
                temp.discard(peptide)
            else:
                print('Checking for consistency')
                partialSpec=LinearSpectrum(peptide)
                print(partialSpec)
                for i in range(len(partialSpec)):
                    if partialSpec[i] not in Spectrum:
                        print('Discarded')
                        temp.discard(peptide)
                        break
        print(CandidatePeptide)
        CandidatePeptide=copy.deepcopy(temp)
        print('reitterating')
        print(CandidatePeptide)
    FinalPeptide=sorted(FinalPeptide)
    return FinalPeptide[::-1]

# This is a helper method to read a spectrum
# from a text file
def SpecFromText(filename,dest=0):
    spec=[]
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    for line in f:
        entryList=list(line.split(' '))
        for entry in entryList:
            spec.append(int(entry.strip()))
    return spec

# Thie method scores a given cyclic peptide and a mass spectrum
# the higher the score, the more the two match
# Input: A string pep and a list of int masses spec
# Output: A int count that represents the score
def PeptideScore(pep,spec):
    count=0
    pep_spec=CyclicSpectrum(pep)
    specDict={}
    #Generates a count of all masses in the spectrum
    for mass in spec:
        if mass not in specDict:
            specDict[mass]=1
        else:
            specDict[mass]+=1
    # Generates a count of all masses in the peptide
    pepDict={}
    for mass in pep_spec:
        if mass not in pepDict:
            pepDict[mass] = 1
        else:
            pepDict[mass] += 1

    # Begin Score
    for keys in pepDict:
        if keys in specDict:
            count+=min(pepDict[keys],specDict[keys])
    return count

# Thie method scores a given linear peptide and a mass spectrum
# the higher the score, the more the two match
# Input: A string pep and a list of int masses spec
# Output: A int count that represents the score
def LinearPeptideScore(pep,spec):
    count=0
    pep_spec=LinearSpectrum(pep)
    specDict={}
    #Generates a count of all masses in the spectrum
    for mass in spec:
        if mass not in specDict:
            specDict[mass]=1
        else:
            specDict[mass]+=1
    # Generates a count of all masses in the peptide
    pepDict={}
    for mass in pep_spec:
        if mass not in pepDict:
            pepDict[mass] = 1
        else:
            pepDict[mass] += 1

    # Begin Score
    for keys in pepDict:
        if keys in specDict:
            count+=min(pepDict[keys],specDict[keys])
    return count

# This method trims a list so that only the top n peptides
# are kept
# Input: A set of strings Leaderboard containing peptides, a
# list Spectrum that contains the masses of peptide fragments,
# and a positive int n that represents the top ranking
# Output A list trim_Leaderboard
import pandas as pd
def Trim(Leaderboard,Spectrum,n):
    print("Making Dataframe...")
    d={'peptide':[],'score':[]}
    index_List=[]
    for i in range(len(Leaderboard)):
        d['peptide'].append(Leaderboard[i])
        d['score'].append(PeptideScore(Leaderboard[i],Spectrum))
        index_List.append(i)
    df=pd.DataFrame(data=d,index=index_List)
    print("Scoring...")
    df['score_rank']=df['score'].rank(ascending=0)
    print("Eliminating...")
    trim_Leaderboard=[]
    for index in index_List:
        if df.loc[index,'score_rank']<=n:
            print(df.loc[index,'score_rank'])
            trim_Leaderboard.append(str(df.loc[index,'peptide']))
    return trim_Leaderboard

# This method sequences a cyclopeptide with a noisy spectrum
# Input: A list of int masses Spectrum and a positive int N
# Output: A string with the best scoring peptide with the given 
# mass spectrum
def LeaderboardCyclopeptideSequencing(Spectrum,N):
    LeaderBoard={''}
    LeaderPeptide=""
    while len(LeaderBoard)>0:
        print('Current Len of Leaderboard: '+str(len(LeaderBoard)))
        print("Expanding...")
        LeaderBoard=PeptideExpander(LeaderBoard)
        LeaderBoard.discard("")
        print("Deep copying...")
        temp=copy.deepcopy(LeaderBoard)
        print("Done copying...")
        for peptide in LeaderBoard:
            if PeptideMass(peptide) == Spectrum[-1]:
                if PeptideScore(peptide,Spectrum) > PeptideScore(LeaderPeptide,Spectrum):
                    LeaderPeptide=peptide
            elif PeptideMass(peptide) > Spectrum[-1]:
                 print("Discarded")
                 temp.discard(peptide)
        print("Trimming")
        LeaderBoard=Trim(list(temp),Spectrum,N)
        print('Current Len of Leaderboard: '+str(len(LeaderBoard)))
        print()
    print(LeaderPeptide[::-1])
    print(PeptideToMass(LeaderPeptide[::-1]))
    return LeaderPeptide[::-1]

# This method generates a convoluted spectrum
# from a given mass spectrum. Convolution of a spectrum is done by taking 
# all positive differences of masses in the spectrum and storing them
# Input: A list of int masses spectrum
# Output: The list of elements in the convolution of spectrum 
# If an element has multiplicity k, it should appear exactly k times
def SpectrumConvolution(spectrum):
    convolution={}
    for i in range(1,len(spectrum)):
        for j in range(0,i):
            diff=spectrum[i]-spectrum[j]
            if diff>=57 and diff <=200:
                if diff not in convolution:
                    convolution[diff]=1
                else:
                    convolution[diff]+=1
    # out=""
    # for keys in convolution:
    #     for i in range(convolution[keys]):
    #         out+=str(keys)+' '
    # print(out)
    return convolution

# This method sequences a cyclopeptide using convolution
# This accounts for missing masses in a noisy experimental 
# spectrum
# Input: A positive  int m that selects the top m masses in a
# convoluted spectrum, a positive int n that selects the top n
# scoring peptides, and a list of ints Spectrum that is the mass
# spectrum of a protein
# Output: A string LeaderPeptide that is the best scoring peptide
# of a given spectrum
def ConvolutionCyclopeptideSequencing(M,N,Spectrum):
    convolution = SpectrumConvolution(Spectrum)
    #Eliminating infrequent Convolutions
    trim_convolutions=[]
    print("Making Dataframe...")
    d = {'mass': [], 'freq': []}
    index_List = []
    i=0
    for keys in convolution:
        d['mass'].append(keys)
        d['freq'].append(convolution[keys])
        index_List.append(i)
        i+=1
        print(str(keys)+":"+str(convolution[keys]))

    df = pd.DataFrame(data=d, index=index_List)
    print("Scoring...")
    df['score_rank'] = df['freq'].rank(ascending=0)
    print("Eliminating...")
    trim_Leaderboard = []
    for index in index_List:
        if df.loc[index, 'score_rank'] <= M:
            trim_convolutions.append(int(df.loc[index, 'mass']))

    LeaderBoard = {''}
    LeaderPeptide = ""
    while len(LeaderBoard) > 0:
        print('Current Len of Leaderboard: ' + str(len(LeaderBoard)))
        print("Expanding...")
        LeaderBoard = LimitedPeptideExpander(LeaderBoard,trim_convolutions)
        LeaderBoard.discard("")
        print("Deep copying...")
        temp = copy.deepcopy(LeaderBoard)
        print("Done copying...")
        for peptide in LeaderBoard:
            print(peptide)
            if PeptideMass(peptide) == Spectrum[-1]:
                if PeptideScore(peptide, Spectrum) > PeptideScore(LeaderPeptide, Spectrum):
                    LeaderPeptide = peptide
            elif PeptideMass(peptide) > Spectrum[-1]:
                print("Discarded")
                temp.discard(peptide)
        print("Trimming")
        LeaderBoard = Trim(list(temp), Spectrum, N)
        print('Current Len of Leaderboard: ' + str(len(LeaderBoard)))
        print()
    print(LeaderPeptide[::-1])
    print(PeptideToMass(LeaderPeptide[::-1]))
    return LeaderPeptide
