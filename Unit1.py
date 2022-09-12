#Code is from Finding Hidden Messages in DNA (Bioinformatics I)

#PatternCount
#A method used to count the number of times Pattern occurs in Text
#Inputs: Text- a string and Pattern a string
#Output; Count: The number of times Pattern occurs in Text
def PatternCount(Text,Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        kmer=Text[i:i+len(Pattern)]
        if kmer==Pattern:
            count+=1
    print(count)
    return count

#A method used to determine the number of
#inputs: Text is a string, k is an int
#output: freq: a dictionary that contains each unique kmer within the 
#text and the number of times the kmer appears in the Text
def FrequencyTable(Text,k):
    freq={}
    for i in range(len(Text)-k+1):
        kmer=Text[i:i+k]
        if kmer in freq:
            freq[kmer]+=1
        else: freq[kmer]=1
    #print(freq)
    return freq

#This method returns the highest value in a dictionary
#Input: dict- dictionary
#Output: max_value - the maximum value in the dictionary
def MaxMap(dict):
    max_value=-1
    for keys in dict:
        if dict[keys]>max_value:
            max_value=dict[keys]
    return max_value

#This method 
#Input:Text: string 
#Output: returns a list freq_patterns that contains 
def FrequentWords (Text,k):
    freq_patterns=[]
    map=FrequencyTable(Text,k)
    max_of_map=MaxMap(map)
    for keys in map:
        if map[keys]==max_of_map:
            freq_patterns.append(keys)
    return freq_patterns

def PrintList(output):
    for i in range(len(output)):
        print(output[i])

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

def PatternMatch(pattern,genome):
    matches=[]
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)]==pattern or RevComp(genome[i:i+len(pattern)])==pattern:
            matches.append(i)
    PrintList(matches)
    return matches

def FindClumps(text,k, L,t):
    pattern=[]
    n=len(text)
    for i in range(n-L):
        window=text[i:i+L]
        map = FrequencyTable(window,k)
        for keys in map:
            if map[keys]>=t:
                pattern.append(keys)
    pattern=list(dict.fromkeys(pattern)) #This removes duplicates by converting to dict and back to list
    PrintList(pattern)
    return pattern

#Week 2

def Skew(genome):
    genome=genome.upper()
    skew=[0]
    for i in range(len(genome)):
        nucleotide=genome[i]
        curr_val=skew[i]
        if nucleotide=="C":
            curr_val-=1
        elif nucleotide=='G':
            curr_val+=1
        skew.append(curr_val)
    #PrintListFlat(skew)
    return skew

def PrintListFlat(list):
    out=''
    for  i in range(len(list)):
        out+= str(list[i])+" "
    print(out)

def MinSkew(genome):
    skew=Skew(genome)
    out=[]
    min_value = min(skew)
    for i in range(len(skew)):
        if skew[i]==min_value:
            out.append(i)
    PrintListFlat(out)
    return out

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

def ApproxPatternMatch(pattern,text,d):
    locations=[]
    n=len(text)
    k=len(pattern)
    for i in range(n-k+1):
        if Hamming(pattern,text[i:i+k]) <=d :
            locations.append(i)
    #PrintListFlat(locations)
    return locations

def ApproxPatternCount(text,pattern,d):
    count=0
    n=len(text)
    k=len(pattern)
    for i in range(n-k+1):
        if Hamming(text[i:i+k],pattern) <=d:
            count+=1
    #print(count)
    return count

def Neightbors(pattern,d):
    #Base cases of the recursive alg
    if d == 0:
        return pattern
    if len(pattern)==1:
        return ['A','C','G','T']

    Neightborhood=[]
    SuffixNeighbors=Neightbors(pattern[1::],d)
    for i in range(len(SuffixNeighbors)):
        text=SuffixNeighbors[i]

        if Hamming(text,pattern[1::]) < d:
            for nucleotide in ['A','C','G','T']:
                Neightborhood.append(nucleotide+text)
        else:
            Neightborhood.append(pattern[0]+text)
    #PrintList(Neightborhood)
    #print()
    return Neightborhood

def FreqWordsWithMismatches(text,k,d):
    pattern = []
    freqmap={}
    n=len(text)

    for i in range(n-k+1):
        kmer=text[i:i+k]
        neighborhood=Neightbors(kmer,d)
        for j in range(len(neighborhood)-1):
            neighbor=neighborhood[j]
            if neighbor not in freqmap:
                freqmap[neighbor]=1
            else:
                freqmap[neighbor]+=1

    max_freq=MaxMap(freqmap)
    for key in freqmap:
        if freqmap[key]== max_freq:
            pattern.append(key)
    PrintListFlat(pattern)
    return pattern

def FreqWordWithMismatchAndRevComp(text,k,d):
    #Initialize variables for this method
    pattern = []
    freqmap = {}
    n = len(text)

    #Begins the search
    for i in range(n - k + 1):
        kmer = text[i:i + k]
        neighborhood = Neightbors(kmer, d)
        neighborhood.append(Neightbors(RevComp(kmer),d))
        for j in range(len(neighborhood) - 1):
            neighbor = neighborhood[j]
            if neighbor not in freqmap:
                freqmap[neighbor] = 1
            else:
                freqmap[neighbor] += 1

    for i in range(n - k + 1):
        kmer = RevComp(text[i:i + k])
        neighborhood = Neightbors(kmer, d)
        neighborhood.append(Neightbors(RevComp(kmer),d))
        for j in range(len(neighborhood) - 1):
            neighbor = neighborhood[j]
            if neighbor not in freqmap:
                freqmap[neighbor] = 1
            else:
                freqmap[neighbor] += 1

    #Finds the most frequent and saves to final output
    max_freq = MaxMap(freqmap)
    for key in freqmap:
        if freqmap[key] == max_freq:
            pattern.append(key)
    PrintListFlat(pattern)
    return pattern

# Week 3
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

def MotifEnumeration(DNA,k,d):
    patterns=[]
    kmers=[]
    #Get all kmers
    for i in range(len(DNA)):
        for j in range(len(DNA[i])-k+1):
            kmers.append(DNA[i][j:j+k])

    #Generate mismatches
    neighborhood=[]
    for kmer in kmers:
        neighborhood+=Neightbors(kmer,d)

    #Check if in all seq
    for neighbor in neighborhood:

        isInAll = []
        for i in range(len(DNA)):
            isInAll.append(False)

        for i in range(len(DNA)):
            for j in range(len(DNA[i])-k+1):
                if Hamming(neighbor,DNA[i][j:j+k]) <= d:
                    isInAll[i]=True
        add= True
        for entry in isInAll:
            if entry == False:
                add=False
        if add:
            patterns.append(neighbor)

    #Removes the duplicates
    patterns=list(dict.fromkeys(patterns))

    #Prints the results
    #PrintListFlat(patterns)
    return patterns

def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(Matrix[i][j]) + " "
        print(out)

def DistBetweenPatternAndString(Pattern,DNA):
    k=len(Pattern)
    distance=0
    for line in DNA:
        ham=1000000000000000
        for i in range(len(line)-k+1):
            text=line[i:i+k]
            if ham > Hamming(text,Pattern):
                ham=Hamming(text,Pattern)
        distance+=ham
    #print(distance)
    return distance

def AllStrings(k):
    pattern=""
    for i in range(k):
        pattern+="A"
    return Neightbors(pattern,len(pattern))

import math
def MedianString(DNA,k):
    distance=math.inf
    patterns=AllStrings(k)
    for i in range(len(patterns)):
        pattern=patterns[i]
        if distance > DistBetweenPatternAndString(pattern,DNA):
            distance=DistBetweenPatternAndString(pattern,DNA)
            median=pattern
    print(median)
    return median

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j].upper()
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    profile = Count(Motifs)
    k = len(Motifs)
    t = len(profile["A"])
    for i in range(t):
        profile["A"][i] =profile["A"][i]/k
        profile["C"][i] =profile["C"][i]/k
        profile["G"][i] =profile["G"][i]/k
        profile["T"][i] =profile["T"][i]/k
    return profile

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    for j in range(len(Motifs)):
        for i in range(len(consensus)):
            if consensus[i] != Motifs[j][i]:
                score += 1
    return score

def Pr(Text, Profile):
    prob = 1.0
    for i in range(len(Text)):
        prob = prob * Profile[Text[i].upper()][i]
    return prob

# Then copy your ProfileMostProbablePattern(Text, k, Profile) and Pr(Text, Profile) functions here.
def ProfileMostProbablePattern(Text, k, Profile):
    #Generate Dictionary of all the kmers in the text
    kmerDict = {}
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        if kmer not in kmerDict:
            kmerDict[kmer] = 1
    #Calculate the probability
    for keys in kmerDict:
        kmerDict[keys] = Pr(keys, Profile)

    m = max(kmerDict.values())
    for keys in kmerDict:
        if kmerDict[keys] == m:
            return keys

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs with pseudocounts, as a dictionary of lists.
def PseudoProfile(Motifs):
    profile = Count(Motifs)
    k = len(Motifs)+4
    t = len(profile["A"])
    for i in range(t):
        profile["A"][i] =(profile["A"][i]+1)/k
        profile["C"][i] =(profile["C"][i]+1)/k
        profile["G"][i] =(profile["G"][i]+1)/k
        profile["T"][i] =(profile["T"][i]+1)/k
    return profile

def ImprovedGreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = PseudoProfile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs

###### WEEK 4
def Motifs(Profile,DNA):
    motif=[]
    k=len(Profile['A'])
    for i in range(len(DNA)):
        motif.append(ProfileMostProbablePattern(DNA[i],k,Profile))
    return motif

# pro= {'A': [4/5,0,0,1/5],'C':[0,3/5,1/5,0], 'G':[1/5,1/5,4/5,0],"T":[0,1/5,0,4/5]}
# print(pro)
# DNA=[
#     "TTACCTTAAC",
#     "GATGTCTGTC",
#     'ACGGCGTTAG',
#     'CCCTAACGAG',
#     'CGTCAGAGGT'
# ]
# print(Motifs(pro,DNA))

import random
def RandomizedMotifSearch(DNA, k ,t):
    temp_Motif=[]
    #Randomly Chooses a set of k-mers from all strings
    for i in range(t):
        var= random.randint(0,len(DNA[0])-k)
        temp_Motif.append(DNA[i][var:var+k])
    bestMotif=temp_Motif[:]

    #Begin Iteration
    while True:
        p=PseudoProfile(temp_Motif)
        temp_Motif=Motifs(p,DNA)[:]
        if Score(temp_Motif)<Score(bestMotif):
            bestMotif=temp_Motif[:]
        else:
            return bestMotif

    return bestMotif

def MultiRandomSearch(DNA,k,t,n):
    bestScore=math.inf
    for i in range(n):
        print("i is: "+str(i))
        motif=RandomizedMotifSearch(DNA,k,t)
        if Score(motif)<bestScore:
            print(Score(motif))
            bestScore=Score(motif)
            out=motif
    return out

import math
def Entropy(Motif):
    profile=Profile(Motif)
    entropy=0
    for j in range(len(profile['A'])):
        tempEnt=0
        for base in ['A','C','G','T']:
            prob=profile[base][j]
            if prob != 0:
                tempEnt+=prob*math.log2(prob)
        print(-tempEnt)
        entropy+=-tempEnt
    print()
    print(entropy)
    return entropy

p1 = float ( (600.-15) / (600.-15+1)   )
p2 = 1 - p1
from itertools import *
counter = 0
for seq in combinations(range(10),2):
    counter +=1
#print(pow(p2,2) * pow(p1,8) * counter)

def Random(prob):
    sum=0
    for i in range(len(prob)):
        sum+=prob[i]

    if sum>1:
        for i in range((len(prob))):
            prob[i]/=sum

    sum=round(sum)*100
    num = random.randint(0,sum)/100.
    distribution=0
    for i in range(len(prob)):
        if distribution<= num and num < distribution+prob[i]:
            return i
        else:
            distribution+=prob[i]

def ProfileRandomPattern(Text, k, Profile):
    #Generate Dictionary of all the kmers in the text
    kmerDict = {}
    for i in range(len(Text) - k + 1):
        kmer = Text[i:i + k]
        if kmer not in kmerDict:
            kmerDict[kmer] = 1
    kmerList=list(kmerDict)

    #Calculate the probability
    prob=[]
    for kmer in kmerList:
        prob.append(Pr(kmer, Profile))
    return random.choices(kmerList,prob)[0]


def GibbsSampler(DNA,k,t,n):
    temp_Motif = []
    # Randomly Chooses a set of k-mers from all strings
    for i in range(t):
        var = random.randint(0, len(DNA[0]) - k)
        temp_Motif.append(DNA[i][var:var + k])
    bestMotif = temp_Motif[:]

    for j in range(1,n):
        i = random.randint(0,t-1)
        temp_Motif.pop(i)
        profile=PseudoProfile(temp_Motif)
        temp_Motif.insert(i,ProfileRandomPattern(DNA[i],k,profile))
        if Score(temp_Motif)< Score(bestMotif):
            #print("Score of new Motif")
            #print(Score(temp_Motif))
            #print()
            bestMotif=temp_Motif[:]

    return bestMotif

def MultiGibbs(DNA,k,t,n,m):
    bestScore=math.inf

    for i in range(m):
        print("i is: "+str(i))
        motif=GibbsSampler(DNA,k,t,n)
        if Score(motif)<bestScore:
            print('Score of best motif so far: ' +str(Score(motif)))
            bestScore=Score(motif)
            out=motif
    return out


