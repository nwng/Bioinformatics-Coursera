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

def MatrixPrint(Matrix):
    for i in range(len(Matrix)):
        out=''
        for j in range(len(Matrix[i])):
            out+=str(round(Matrix[i][j],1)) + " "
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


#FarthestFirstTraversal('trav.txt')

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

def SqErrorDistortion(centers,points):
    distort=0
    for point in points:
        distort+=(Distance(point,centers))**2
    distort/=len(points)
    return distort

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

# k,m,points = PointsFromFile('trav.txt')
# LLoyd(k,m,points)

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

#print(kMeansInit(points,k))

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


# data=[[2, 6], [4, 9], [5, 7], [6, 5], [8, 3] ]
# centers=[ [4, 5], [7, 4]]
# print("Max Distance")
# print(MaxDistance(data,centers))
# print()
#
# data1=[[2, 6], [4, 9], [5, 7], [6, 5], [8, 3] ]
# centers1=[[4, 5], [7, 4]]
# print("Distortion")
# print(SqErrorDistortion(centers1,data1))
# print()
#
# print("Center of Gravity")
# points=[[17, 0, -4], [3, 14, 23], [9, 7, 16], [7, 3, 5]]
# print(CenterOfGravity(points))

def Pr(n,Data,Theta):
    return Theta**(n*Data)*(1-Theta)**(n*(1-Data))


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


#matrix= MatrixFromText('trav.txt')
#HierarchicalClustering(matrix,20)
# quest1=(0.7**5)*(0.3**1)
# print(round(quest1,3))
# print()
#
# d1=(0**2)+(2**2)
# d2=(2**2)+(3**2)
# prob=(1/d1)/((1/d1)+(1/d2))
# print(round(prob,3))
# print()
#
# num = (2*.6) + (4*.1) + (5*.8) + (6*.5) + (8*.7)
# denom=.6+.1+.8+.5+.7
# print(round(num/denom,3))
# num = (6*.6) + (9*.1) + (7*.8) + (5*.5) + (3*.7)
# print(round(num/denom,3))
# print()

print(.75**4)