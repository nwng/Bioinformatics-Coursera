# The following code was implemented to solve problems from 
# Comparing Genes, Proteins, and Genomes (Bioinformatics III)


###### The following are helper methods from previous chapters
# This method reads a matrix from a text file and
# saves it as a list of lists.
def MatrixFromText(filename, dest=0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    matrix=[]
    for line in f:
        rawLine=line.rsplit(' ')
        toadd=[]
        for entry in rawLine:
            toadd.append(int(entry.strip()))
        matrix.append(toadd)
    return matrix

# This method helps print a list of lists
def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(Matrix[i][j]) + " "
        print(out)

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
####

# This method reads a directed acyclic graph from a text file
# and saves it as a dictionary
def DAGfromFile(filename,top=0,dest=0):
    if top==1:
        topList=[]

    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    map={}
    topList=[]
    for line in f:
        rawLine=line.rsplit(' ')
        startNode=int(rawLine[0])
        endNode=int(rawLine[1])
        weight=int(rawLine[2].strip())
        if startNode not in map:
            map[startNode]={}
            map[startNode][endNode]=weight
            topList.append(startNode)
            #print(map)
        else:
            map[startNode][endNode] = weight
            #print(map)
    if top==1:
        #print(topList)
        return map,topList
    return map

# This method reads a directed acyclic graph from a text file
# and saves it as a list of lists
def ModDAGfromFile(filename,dest=0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    map=[]
    for line in f:
        rawLine=line.rsplit(' ')
        startNode=int(rawLine[0])
        endNode=int(rawLine[1])
        weight=int(rawLine[2].strip())
        map.append((startNode,endNode,weight))
    return map

# This method is a dynamic programming approach to calculate the
# minimum number of coins needed for an int amount of money
# Input: A positive int money and a list of int coins
# Ouput: A int of the minimum number of coins
import math
def DPChange(money,coins):
    minNum=[0]
    for i in range(1,money+1):
        #print('i is: '+str(i))
        curr_num= math.inf
        for j in range(len(coins)):
            #print("curr coin: "+ str(coins[j]))
            if i >= coins[j]:
                if minNum[i-coins[j]] + 1 < curr_num:
                    curr_num=minNum[i-coins[j]] + 1
        minNum.append(curr_num)
        #print(minNum)
    return minNum[money]

# This method finds a maximum-weight path connecting the source to the sink
# in an n x m grid
# Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. 
# Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid 
# whose edges are defined by the matrices Down and Right.
def ManhattanTourist(n, m, Down, Right):
    s=[[0]*(m+1) for i in range(n+1)]
    # print(s)
    for i in range(1,n+1):
        # print('i is: '+str(i))
        s[i][0]=s[i-1][0]+Down[i-1][0]
        # print(s)
    for j in range(1,m+1):
        # print('j is: ' + str(j))
        s[0][j]=s[0][j-1]+Right[0][j-1]
        # print(s)
    for i in range(1,n+1):
        for j in range(1,m+1):
            s[i][j]=max((s[i-1][j]+Down[i-1][j]),s[i][j-1]+Right[i][j-1])
    #MatrixPrint(s)
    return s[-1][-1]

# LCSBacktrack produces a matrix that store pointers that reconstruct
# the longest common string between strings v and w
# Input: Two strings v and w
# Output: A list of lists backTrack matrix
def LCSBacktrack(v,w):
    s = [[0] * (len(w)+1) for i in range(len(v)+1)]
    backTrack = [[0] * (len(w)+1) for i in range(len(v)+1)]
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            match=0
            if v[i-1]==w[j-1]:
                match=1
            s[i][j]=max(
                s[i-1][j],
                s[i][j-1],
                s[i-1][j-1] + match
            )
            if s[i][j]==s[i-1][j]:
                backTrack[i][j]="V"
            elif s[i][j]==s[i][j-1]:
                backTrack[i][j]='>'
            elif s[i][j]==s[i-1][j-1]+match:
                backTrack[i][j]="\\"
    #MatrixPrint(s)
    #print()
    #MatrixPrint(backTrack)
    return backTrack

import sys
sys.setrecursionlimit(1500)

# This method reconstructs the longest common string
# Input: A list backtrack that has pointers, string v,
# and two ints i and j
# Output: A string that is the longest common string
def OutputLCS(backtrack,v,i,j):
    if i==0 or j==0:
        return ""

    if backtrack[i][j]=="V":
        return OutputLCS(backtrack,v,i-1,j)
    elif backtrack[i][j]==">":
        return OutputLCS(backtrack,v,i,j-1)
    else:
        return OutputLCS(backtrack,v,i-1,j-1) + v[i-1]

# This method returns the longest path in a directed
# acyclic graph
# Input: Two ints startNode and endNode and a string
# filename
# Output: An int that is the length of the longest path
# and a list that shows the order of nodes from the start
# node to the end node
import copy
def LongestPathInDAG(startNode,endNode,filename):
    travelList = ModDAGfromFile(filename)
    travelDict= DAGfromFile(filename)
    allNodes={}
    map = DAGfromFile(filename)
    b={}
    for keys in map:

        b[keys]=[]
        for subkeys in map[keys]:
            b[subkeys]=[]

    #Finds all previous nodes
    for nodes in map:
        for outDegree in map[nodes]:
            b[outDegree].append(nodes)

    b[startNode]=[startNode]

    for entry in travelList:
        allNodes[entry[0]]=0
        allNodes[entry[1]] = 0
    allNodes=list(allNodes.keys())
    allNodes.sort()
    s={}
    for node in allNodes:
        s[node]=-math.inf

    s[startNode]=0
    b[startNode]=[startNode]
    copy_original_b=copy.deepcopy(b)
    # print(s)
    # print(b)

    for entry in travelList:
        previousNode=entry[0]
        currNode=entry[1]
        weight= entry[2]
        if s[currNode]<s[previousNode]+weight:
            s[currNode]=s[previousNode]+weight
            b[currNode] = previousNode
        #checks in case updated
        # print(currNode)
        if len(copy_original_b[currNode])>1:
            for previousNode in copy_original_b[currNode]:
                weight=travelDict[previousNode][currNode]
                if s[currNode] < s[previousNode] + weight:
                    s[currNode] = s[previousNode] + weight
                    b[currNode] = previousNode

    print(s[endNode])
    path=[]
    currNode=endNode
    while currNode!=startNode:
        path.append(currNode)
        currNode=b[currNode]
    path.append(startNode)
    PrintListFlat(path[::-1])
    return s[endNode], path[::-1]


BLOSUM62={'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2,
                'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2},
          'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1,
                'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},
          'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3,
                'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
          'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1,
                'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3},
          'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0,
                'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},
          'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3,
                'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},
          'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3,
                'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},
          'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1,
                'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},
          'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0,
                'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
          'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2,
                'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1},
          'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3,
                'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1},
          'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6,
                'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},
          'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0,
                'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1},
          'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3,
                'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3},
          'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1,
                'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},
          'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0,
                'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2},
          'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1,
                'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2},
          'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2,
                'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},
          'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3,
                'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},
          'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1,
                'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}



