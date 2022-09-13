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

#This method finds the frequent k lengthed substring patterns in a string Text
#Input:Text: string and k an int that is the length ot the 
#pattern sought for
#Output: returns a list freq_patterns that contains 
def FrequentWords (Text,k):
    freq_patterns=[]
    map=FrequencyTable(Text,k)
    max_of_map=MaxMap(map)
    for keys in map:
        if map[keys]==max_of_map:
            freq_patterns.append(keys)
    return freq_patterns

#This method prints a list, each element in the list is printed in a new line
def PrintList(output):
    for i in range(len(output)):
        print(output[i])

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

#This method finds instances of a pattern in a genome. This seaches for both the
#pattern and its reverse complement in the genome
#Input: Pattern - a string and Genome another string
#Output: A list containing the index of where the pattern appears in the genome
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

# This method calculates the skew of a DNA string.
# Skew is a measure of the G and C content of a DNA string
# and is defined by +1 for G and -1 for C. The skew can be used
# to find the origin of replication in a bacterial genome
# Input: A string genome
# Output: A list that shows the skew at a given index of the genome

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

#This method combines all elements in a list into a single string
#and prints the string on a single line
def PrintListFlat(list):
    out=''
    for  i in range(len(list)):
        out+= str(list[i])+" "
    print(out)

# This method takes a genome and returns the location where the skew 
# is the smallest
# Input: A string genome
# Output: a list containing the indices in genome where skew is the smallest
def MinSkew(genome):
    skew=Skew(genome)
    out=[]
    min_value = min(skew)
    for i in range(len(skew)):
        if skew[i]==min_value:
            out.append(i)
    PrintListFlat(out)
    return out

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

# This method finds when a pattern appears in a text with at most d
# mismatches
# Input: pattern is a string, text is a string, and d is an int
# Output: A list locations that contains the index of where the pattern
# can be found in the text
def ApproxPatternMatch(pattern,text,d):
    locations=[]
    n=len(text)
    k=len(pattern)
    for i in range(n-k+1):
        if Hamming(pattern,text[i:i+k]) <=d :
            locations.append(i)
    #PrintListFlat(locations)
    return locations

# This method counts the instances a pattern appears in a text with at most d
# mismatches
# Input: pattern is a string, text is a string, and d is an int
# Output: An int count of the number of times a pattern appears in a text
def ApproxPatternCount(text,pattern,d):
    count=0
    n=len(text)
    k=len(pattern)
    for i in range(n-k+1):
        if Hamming(text[i:i+k],pattern) <=d:
            count+=1
    #print(count)
    return count

# This is a recursive method that returns all possible patterns
# with at most d differences from pattern
# Input: A string pattern and a positive integer d
# Output: A list containing the patterns and strings that have at
# most d mistmatches from the pattern
def Neightbors(pattern,d):
    #Base cases of the recursive alg
    if d == 0: #If distance is 0, return original pattern
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

# This method finds the most frequent patterns of length k in a
# string text with at most d mismatches
# Input: A string text, a int k, and an int d
# Output: A list with all patterns that appear the most frequently
# in a string with at most d mismatches
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

# This method finds the most frequent patterns of length k in a
# string text with at most d mismatches. This method also checks
# if the pattern appears as a reverse complement in the text
# Input: A string text, a int k, and an int d
# Output: A list with all patterns that appear the most frequently
# in a string with at most d mismatches
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

# This method seaches for patterns that have a length k that
# appear with at most d mismatches in DNA sequences in a list 
# of DNA sequences
# Input: DNA is a list of dna sequences, k is a positive int, 
# and d is a positive int
# Output: A list of all patterns that appear in all DNA strings with
# at most d mimatches in a list of DNA sequences
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

# This method helps print a list of lists
def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(Matrix[i][j]) + " "
        print(out)


# This method calculates the distance between a pattern
# and a list of strings. Distance is defined as the sum of
# all minimum hamming distances of the pattern and a substring
# of DNA
# Input: A string pattern and a list of strings DNA
# Output: An int distance
import math
def DistBetweenPatternAndString(Pattern,DNA):
    k=len(Pattern)
    distance=0
    for line in DNA:
        ham= math.inf
        for i in range(len(line)-k+1):
            text=line[i:i+k]
            if ham > Hamming(text,Pattern):
                ham=Hamming(text,Pattern)
        distance+=ham
    #print(distance)
    return distance

# This method generates all possible DNA strings of length k
# Input: A positive int k
# Output: A list of every possible DNA sequence of length k
def AllStrings(k):
    pattern=""
    for i in range(k):
        pattern+="A"
    return Neightbors(pattern,len(pattern))


# This method find a k-mer Pattern that minimizes the hamming distance between a Pattern and a list of 
# DNA strings over all k-mers Patterns. This k-mer pattern is known as a median string for DNA
# Input: DNA is a list of strings and k is an int that represents the length of the sought pattern
# Output: A string median that is length k that has the minimum distance between the pattern and all
# DNA strings
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

# This method generates a count matrix that counts the number of
# instances a nucleotide appears at a given index in a series of strings
# Input: Motifs is a list of DNA strings of equal length
# Output: A dictionary count 
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

# This method generates a consensus string from a list of DNA
# Sequences
# Input: Motif is a list of strings of DNA
# Output: A consensus string

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

# This is a helper method for a Greedy Search algorithm
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

# This is a helper method for a Greedy Search algorithm
# Input: text is a string that represents a patterns and Profile is a
# list of lists that represents the probability of a given nucleotide
# at a certain index
# Output: The probability a text would appear in a string based on a
# a probability profile
def Pr(Text, Profile):
    prob = 1.0
    for i in range(len(Text)):
        prob = prob * Profile[Text[i].upper()][i]
    return prob

# This method generates the most likely k length pattern that appears in 
# a list of strings Text with a given probability matrix Profile.
# This is a helper method for a Greedy Search algorithm
# Input: Text is a list of strings, k is a positive int, and Profile
# a list of lists. 
# Output: returns a string of length k that is the most probable pattern
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

# This is a Greedy Search algorithm to search for hidden motifs found in a collection of strings
# That represent DNA.
# Input:  A list of kmers Dna, and integers k and t (where k is the length of the pattern
# and t is the number of kmers in Dna)
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

# This method generates pseudoprofiles which is similar to a normal profile, but each
# entry has 1 added to it to account for unlikely events
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

# This is a Greedy Search algorithm to search for hidden motifs found in a collection of strings
# That represent DNA. It is improved by using pseudocounts.
# Input:  A list of kmers Dna, and integers k and t (where k is the length of the pattern
# and t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
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

# This method generates motifs based on a profile matrix and a collection of
# strings
# This method is a helper method for a Monte Carlo random algorithm
# Input: A list profile and a list of strings DNA
# Output: A list motif that contains a 
def Motifs(Profile,DNA):
    motif=[]
    k=len(Profile['A'])
    for i in range(len(DNA)):
        motif.append(ProfileMostProbablePattern(DNA[i],k,Profile))
    return motif


# This method returns the best motif in a collection of strings DNA
# Input: A list of strings DNA, positive int k, and t is the number of DNA
# strings in the list
# Output: A best motif that approximates a solution

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

# This method runs RandomizedMotifSeach n times, and returns the
# best scoring motif
# Input: A list of strings DNA, positive int k, t is the number of DNA
# strings in the list, and n is the number of iterations
# Output: A best motif that is the best
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

# This method returns a random pattern
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


