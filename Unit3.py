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

def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(Matrix[i][j]) + " "
        print(out)

def PrintList(output):
    for i in range(len(output)):
        print(output[i])

def PrintListFlat(list):
    out=''
    for  i in range(len(list)):
        out+= str(list[i])+" "
    print(out)

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
#print(ModDAGfromFile('map.txt'))


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

# down=MatrixFromText('down.txt')
# right=MatrixFromText('right.txt')
# print(ManhattanTourist(17,11,down,right))

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

def OutputLCS(backtrack,v,i,j):
    if i==0 or j==0:
        return ""

    if backtrack[i][j]=="V":
        return OutputLCS(backtrack,v,i-1,j)
    elif backtrack[i][j]==">":
        return OutputLCS(backtrack,v,i,j-1)
    else:
        return OutputLCS(backtrack,v,i-1,j-1) + v[i-1]
# v=''
# w=''
#
# backtrack=LCSBacktrack(v,w)
# print()
# print(OutputLCS(backtrack,v,len(v),len(w)))


# def LongestPathinDAGv1(startNode,endNode,filename,top=0):
#     if top==1:
#         map,topList=DAGfromFile(filename,1)
#     else:
#         map=DAGfromFile(filename)
#
#     s={}
#     b={}
#     #print(topList)
#     topList.append(endNode)
#     for i in range(len(topList)):
#         s[topList[i]]=-1*math.inf
#         b[topList[i]]=-1
#
#     s[startNode]=0
#     b[startNode]=startNode
#
#     skipList=[-1]
#     for i in range(1,len(topList)):
#         currnode=topList[i]
#         max=-1*math.inf
#         for j in range(i):
#             prevNode=topList[j]
#             if prevNode not in skipList:
#                 if currnode in map[prevNode]:
#                     weight=map[prevNode][currnode]
#                     if max < s[prevNode]+weight:
#                         max=s[prevNode]+weight
#                         backNode=prevNode
#         if max==-1*math.inf:
#             skipList.append(currnode)
#         else:
#             s[currnode]=max
#             b[currnode]=backNode
#
#     #Get Path
#     currnode=endNode
#     path=[]
#
#     while currnode!=startNode:
#         path.append(currnode)
#         currnode=b[currnode]
#     path.append(startNode)
#
#     PrintListFlat(path[::-1])
#     return s[endNode],path[::-1]
# #LongestPathinDAGv1(0,4,'map.txt',1)
#
#
# def ZeroInDegree(map):
#     allNodes={}
#     for entry in map:
#         sourceNode=entry[0]
#         sinkNode=entry[1]
#         allNodes[sourceNode]=0
#         allNodes[sinkNode]=0
#     allNodes=list(allNodes.keys())
#     allNodes.sort()
#     inDegree={}
#     outDegree={}
#     for node in allNodes:
#         inDegree[node]=0
#         outDegree[node]=0
#         for entry in map:
#             if entry[1]==node:
#                 inDegree[node]+=1
#             if entry[0]==node:
#                 outDegree[node]+=1
#
#     possibleSource=[]
#     possibleSink=[]
#     for key in inDegree:
#         if inDegree[key]==0:
#             possibleSource.append(key)
#         if outDegree[key]==0:
#             possibleSink.append((key))
#     return possibleSource, possibleSink
#
#
# def LongestPathinDAGv2(startNode, endNode, filename):
#
#     map = ModDAGfromFile(filename)
#     map.sort()
#     PrintList(map)
#     print()
#     inD,outD=ZeroInDegree(map)
#     inD.remove(startNode)
#     outD.remove(endNode)
#     #print(NoInDegree)
#
#     s = {}
#     b = {}
#
#     for i in range(len(map)):
#         sourceNode=map[i][0]
#         sinkNode = map[i][1]
#         s[sourceNode] = 0
#         b[sourceNode] = -1
#         s[sinkNode] = 0
#         b[sinkNode] = -1
#         if sourceNode in inD:
#             s[sourceNode]=-1*math.inf
#             b[sourceNode]=-1
#
#
#     s[startNode] = 0
#     b[startNode] = startNode
#
#     if endNode not in s:
#         s[endNode]=0
#     if endNode not in b:
#         b[endNode]=-1
#
#
#     print("These have no in degree")
#     print(inD)
#     print("These have no out degree")
#     print(outD)
#     print()
#
#     for i in range(len(map)):
#         entry=map[i]
#         currnode = entry[1]
#         prevNode = entry[0]
#         weight=entry[2]
#
#
#         if currnode not in outD:
#             if s[currnode]<s[prevNode]+weight:
#                     print("node " +str(currnode)+'\'s new value = '+str(s[prevNode]+weight))
#                     s[currnode]=s[prevNode]+weight
#                     b[currnode]=prevNode
#
#     # Get Path
#     currnode = endNode
#     path = []
#     print(s)
#     print(b)
#
#     while currnode != startNode:
#         path.append(currnode)
#         currnode = b[currnode]
#     path.append(startNode)
#
#     print(s[endNode])
#     PrintListFlat(path[::-1])
#     return s[endNode], path[::-1]
#
# #LongestPathinDAGv2(0,44,'map.txt')
#
# def LongestPathinDAGv3(startNode, endNode, filename):
#     travelList= ModDAGfromFile(filename)
#     travel=[]
#     for entry in travelList:
#         if entry[1] not in travel:
#             travel.append(entry[1])
#
#     map = DAGfromFile(filename)
#     b={}
#     for keys in map:
#
#         b[keys]=[]
#         for subkeys in map[keys]:
#             b[subkeys]=[]
#
#     #Finds all previous nodes
#     for nodes in map:
#         for outDegree in map[nodes]:
#             b[outDegree].append(nodes)
#
#     b[startNode]=[startNode]
#
#     #Begin calculating nodes
#     s={}
#     for node in b:
#         s[node]=-1*math.inf
#     s[startNode]=0
#     print(b)
#     for entry in travelList:
#         currNode=entry[1]
#         values=[]
#         node=[]
#         predocessors=b[currNode]
#         for predocessor in predocessors:
#             values.append(s[predocessor]+map[predocessor][currNode])
#             node.append(predocessor)
#         print("Current Node: "+str(currNode))
#         print("we needed these nodes to get here:")
#         print(predocessors)
#         print("Current Scores: ")
#         PrintListFlat(values)
#         best=max(values)
#         print('best value: '+str(best))
#
#         for i in range(len(values)):
#             if values[i]==best:
#                 s[currNode]=best
#                 b[currNode]=[node[i]]
#                 print('best node: '+str(node[i]))
#                 break
#         print()
#
#     path=[]
#     currNode=endNode
#     while currNode!=startNode:
#         path.append(currNode)
#         currNode=b[currNode][0]
#     path.append(startNode)
#
#     print(s)
#     print(b)
#     print()
#     print(s[endNode])
#     PrintListFlat(path[::-1])
# LongestPathinDAGv3(0,44,'map.txt')
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

    # #Troubleshooting
    # keyss=list(copy_original_b.keys())
    # keyss.sort()
    # for key in keyss:
    #     print(str(key)+": "+str(copy_original_b[key]))
    #
    # print()
    # print(travelDict)
    # print()
    # print(s)
    # print(b)
    # print()
    print(s[endNode])
    path=[]
    currNode=endNode
    while currNode!=startNode:
        path.append(currNode)
        currNode=b[currNode]
    path.append(startNode)
    PrintListFlat(path[::-1])


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


def OutputLCS_Global(backtrack,v,i,j):
    if i == 0 and j == 0:
        return ""
    if i == 0 or j == 0:
        return "-"

    if len(v)==len(backtrack)-1:
        if backtrack[i][j]=="V":
            return OutputLCS_Global(backtrack,v,i-1,j) + v[i-1]
        elif backtrack[i][j]==">":
            return OutputLCS_Global(backtrack,v,i,j-1) + '-'
        else:
            return OutputLCS_Global(backtrack,v,i-1,j-1) + v[i-1]
    else:
        if backtrack[i][j]=="V":
            return OutputLCS_Global(backtrack,v,i-1,j) + '-'
        elif backtrack[i][j]==">":
            return OutputLCS_Global(backtrack,v,i,j-1) + v[i-1]
        else:
            return OutputLCS_Global(backtrack,v,i-1,j-1) + v[i-1]


def GlobalAlignmentProblem(match,mismatchPenalty,inDelPenalty,v,w):
    s = [[0] * (len(w)+1) for i in range(len(v)+1)]
    backTrack = [[0] * (len(w)+1) for i in range(len(v)+1)]

    for i in range(1,len(s[0])):
        s[0][i]=s[0][i-1]-inDelPenalty
    for i in range(1,len(s)):
        s[i][0]=s[i-1][0]-inDelPenalty

    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            score=0
            if v[i-1]==w[j-1]:
                score=match
            else:
                score=-mismatchPenalty
            s[i][j]=max(
                s[i-1][j]-inDelPenalty,
                s[i][j-1]-inDelPenalty,
                s[i-1][j-1] + score
            )
            #Begin Constructing the backtrack
            if s[i][j]==s[i-1][j]-inDelPenalty:
                backTrack[i][j]="V"
            elif s[i][j]==s[i][j-1]-inDelPenalty:
                backTrack[i][j]='>'
            elif s[i][j]==s[i-1][j-1]+score:
                backTrack[i][j]="\\"

    MatrixPrint(s)
    print()
    MatrixPrint(backTrack)
    print()
    print(s[-1][-1])
    print(OutputLCS_Global(backTrack,v,len(v),len(w)))
    print(OutputLCS_Global(backTrack, w, len(v), len(w)))

v='TCA'
w='CA'
GlobalAlignmentProblem(2,5,1,v,w)

