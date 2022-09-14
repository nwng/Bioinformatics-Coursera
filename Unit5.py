# This code was implemented to solve problems from 
# Genomic Data Science and Clustering (Bioinformatics V)


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
            toadd.append(float(entry.strip()))
        matrix.append(toadd)
    return matrix

# This method helps print a list of lists
def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(round(Matrix[i][j],1)) + " "
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
##################

# This is a helper method to read points from a textfile
def PointsFromFile(filename, dest=0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    k=0
    m=0
    first=True
    points=[]
    for line in f:
        rawLine=line.rsplit(' ')
        toAdd=[]
        if first:
            first=False
            k=int(rawLine[0].strip())
            m=int(rawLine[1].strip())
        else:
            for entry in rawLine:
                toAdd.append(float(entry.strip()))
            points.append(toAdd)
    return k,m,points

# This method computes the Euclidean Distance between a DataPoint
# and centers and returns the smallest distance of that point from all centers
# Input: A list of ints DataPoint and a list of lists of ints Centers
# Output: A int representing the smallest distance between the datapoints and all
# centers
def Distance(DataPoint,Centers):
    k = len(Centers)
    m=len(Centers[0])
    d = []
    d_temp = 0
    if k!=1:
        for i in range(k):
            center=Centers[i]
            for j in range(m):
                d_temp+=(center[j]-DataPoint[j])**2
            d.append(d_temp**.5)
            d_temp=0
    else:
        for j in range(m):
            d_temp += (Centers[0][j] - DataPoint[j]) ** 2
        d.append(d_temp ** .5)

    return min(d)

# This method reads points from a text file and selects
# centers from the list. The first center is the first point
# in the file and subsequent centers are the furthest euclidian
# distance from all existing centers
# Input: A string filename
# Output: A list of centers and their coordinates
def FarthestFirstTraversal(filename,dest=0):
    #Gets Data from the file
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    k=0
    m=0
    first=True
    points=[]
    for line in f:
        rawLine=line.rsplit(' ')
        toAdd=[]
        if first:
            first=False
            k=int(rawLine[0].strip())
            m=int(rawLine[1].strip())
        else:
            for entry in rawLine:
                toAdd.append(float(entry.strip()))
            points.append(toAdd)

    #Begins to select centers
    centers=[]
    centers.append(points[0])

    while len(centers)<k:
        dist=[]
        for i in range(len(points)):
            dist.append(Distance(points[i],centers))
        maxSize=max(dist)
        for i in range(len(dist)):
            if dist[i]==maxSize:
                centers.append(points[i])

    #Prints the Output of the Centers
    for i in range(len(centers)):
        out=""
        for j in range(len(centers[i])):
            out+=str(centers[i][j])+" "
        print(out)
    return centers

# This method reads points from a text file and calculates
# the the mean squared distance from each data point to its nearest center
# Input: A string filename
# Output: A float distort
def SqErrorDistortionFromFile(filename,dest=0):
    if dest==0:
        location = "/Users/nathanng/PycharmProjects/bioinfo/"
    else:
        location = "/Users/nathanng/Downloads/"
    file = location + filename
    f = open(file, "r")
    k=0
    m=0
    first=True
    arePoints=False
    areCenters=True
    centers=[]
    points=[]
    for line in f:
        if line.strip() == "--------":
            areCenters = False
        else:
            rawLine=line.rsplit(' ')
            toAdd=[]
            if first:
                first=False
                k=int(rawLine[0].strip())
                m=int(rawLine[1].strip())
            elif areCenters:
                center=[]
                for entry in rawLine:
                    center.append(float(entry.strip()))
                centers.append(center)
            else:
                for entry in rawLine:
                    toAdd.append(float(entry.strip()))
                points.append(toAdd)

    distort=0
    for point in points:
        distort+=(Distance(point,centers))**2
    distort/=len(points)
    print(round(distort,3))
    return distort

# This method calculates the the mean squared distance 
# from each data point to its nearest center
# Input: Lists centers and points containing coordinates
# Output: A float distort
def SqErrorDistortion(centers,points):
    distort=0
    for point in points:
        distort+=(Distance(point,centers))**2
    distort/=len(points)
    return distort


# This method calculates the the point whose i-th coordinate 
# is the average of the i-th coordinates of all points
# Input: A list points containing the coordinates for all points
# Output: A list grav listing the coordinates for the center of gravity
import copy
def CenterOfGravity(points):
    n=len(points)
    m=len(points[0])
    grav=[]
    for i in range(m):
        sum=0
        for j in range(n):
            sum+=points[j][i]
        sum/=n
        grav.append(sum)
    return grav

# An implementation of the Lloyd algorithm, a clustering heuristics for the k-Means Clustering Problem. 
#It first chooses k arbitrary distinct points Centers from Data as centers and then iteratively performs 
# the following two steps: Centers to Clusters and Clusters to Centers
# Centers to Clusters: After centers have been selected, assign each data point to the cluster corresponding 
# to its nearest center; ties are broken arbitrarily.
# Clusters to Centers: After data points have been assigned to clusters, assign each cluster’s center of 
# gravity to be the cluster’s new center.

# Input: Integers k and m followed by a set of points Data in m-dimensional space.
# Output: A set Centers consisting of k points (centers) resulting from applying 
# the Lloyd algorithm to Data and Centers, where the first k points from Data are 
# selected as the first k centers.

def LLoyd(k,m,points):
    #Choose initial centers
    centers=copy.deepcopy(points[0:k])

    count=0
    while count<k*2:
        clusters = [[] for i in range(k)]
        #Assign points to clusters
        for point in points:
            dist=[]
            for center in centers:
                dist.append(Distance(point,[center]))
            minDist=min(dist)
            for i in range(len(dist)):
                if dist[i]==minDist:
                    clusters[i].append(point)
        #Choose new Centers
        for i in range(len(clusters)):
            centers[i]=CenterOfGravity(clusters[i])
        count+=1
        #print("Count is: " +str(count))

    for i in range(len(centers)):
        out = ""
        for j in range(len(centers[0])):
            out+=str(round(centers[i][j],3))+" "
        print(out)
    return centers

# This is a method to help initialize a k means clustering by picking k centers one 
# at a time, but instead of choosing the point farthest from those picked so far, it 
# chooses each point at random in such a way that distant points are more likely to be chosen than nearby points
# Input: A list data with coordinates and a positive int k
# Output: A list with coordinates for k centers
import random
def kMeansInit(Data,k):
    centers=[]
    centers.append(Data[random.randint(0,len(Data))])
    print(centers)
    while len(centers) < k:
        dist=[]
        sum=0
        for datapoint in Data:
            value=Distance(datapoint,centers)**2
            dist.append(value)
            sum+=value
        for i in range(len(dist)):
            dist[i]/=sum
        centers.append(random.choices(Data,dist))
    return centers

# This method calculates the farthest distance between a 
# datapoint and all avaliable centers
# Input: A list data and centers
# Output: A float maxDistance

def MaxDistance(data,centers):
    maxDistance=0
    k=len(centers)
    clusters = [[] for i in range(k)]
    # Assign points to clusters
    for point in data:
        dist = []
        for center in centers:
            dist.append(Distance(point, [center]))
        minDist = min(dist)
        for i in range(len(dist)):
            if dist[i] == minDist:
                clusters[i].append(point)

    for i in range(len(clusters)):
        cluster=clusters[i]
        center=centers[i]
        dist=[]
        for point in cluster:
            dist.append(Distance(point,[center]))
        tempMax= max(dist)
        if tempMax>maxDistance:
            maxDistance=tempMax
    return maxDistance

# This returns the conditional probability of generating the outcome Data given a coin with bias θ
# Input: An int n, float Data, and a float Theta
# Output: A float conditional probability
def Pr(n,Data,Theta):
    return Theta**(n*Data)*(1-Theta)**(n*(1-Data))

# This methodprogressively generates n different partitions of the underlying data into clusters, 
# all represented by a tree in which each node is labeled by a cluster of genes. 
# In general, the i-th partition merges the two closest clusters from the (i - 1)-th partition and 
# has n - i + 1 clusters.
# Input: A list of lists D representing a distance matrix and a positive int n
# Output: The result of applying HierarchicalClustering to this distance matrix (using Davg), 
# with each newly created cluster listed on each line.
import math
def HierarchicalClustering(D,n):
    # MatrixPrint(D)
    # print()
    Clusters=[]
    for i in range(n):
        Clusters.append(str(i+1))
    T=[]
    while len(Clusters)>1:
        smallestValue=math.inf
        node=[]
        for i in range(1,len(D)):
            for j in range(0,i):
                # print("i is: " + str(i+1))
                # print("j is: " + str(j+1))
                # print("Dij is: "+str(D[i][j]))
                # print()

                if D[j][i]<smallestValue:
                    smallestValue=D[j][i]
                    node=[j,i,Clusters[j],Clusters[i]]
                    #print(node)
        # MatrixPrint(D)
        # print()
        # T.append([node[3],node[2]])
        Clusters.pop(node[1])
        # print(Clusters)
        Clusters.pop(node[0])
        # print(Clusters)

        tempMerge=[D[node[0]],D[node[1]]]
        toMerge=[]

        numElementsNode0 = len(node[2].split(" "))
        numElementsNode1 = len(node[3].split(" "))

        for i in range(len(tempMerge[0])):
            if i == node[0]:
                toMerge.append(0.00)
            else:
                top = tempMerge[0]
                bot = tempMerge[1]
                numerator=(top[i]*numElementsNode0)+(bot[i]*numElementsNode1)
                totElem=numElementsNode1+numElementsNode0
                toMerge.append(numerator/totElem)

        if len(node[3])<len(node[2]):
            T.append(str(node[3]) + " " + str(node[2]))
            print(str(node[3]) + " " + str(node[2]))
            Clusters.insert(node[0], str(node[3]) + " " + str(node[2]))
        else:
            T.append(str(node[2]) + " " + str(node[3]))
            print(str(node[2]) + " " + str(node[3]))
            Clusters.insert(node[0],str(node[2])+" "+str(node[3]))
        print(Clusters)
        print()

        #Update D
        D.pop(node[1])
        D.pop(node[0])
        D.insert(node[0],toMerge)
        for i in range(len(D)):
            D[i].pop(node[1])
            D[i][node[0]]=toMerge[i]
        MatrixPrint(D)
        print()
    #     print()
    # print()
    # print(Clusters)
    PrintList(T)
