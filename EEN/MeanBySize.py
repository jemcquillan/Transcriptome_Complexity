import sys
from scipy.stats import sem

InFile=open(sys.argv[1], 'r')
OutFile=open(sys.argv[1]+".fmt", 'w')

Exons={}
for line in InFile:
    A=line.split()
#    print A
    gene=A[0]
    trans=A[1]
    Ne=float(A[2])
    exons=int(A[3])
    if exons>1 and exons<20:
        if Exons.has_key(exons):
            Exons[exons].append(Ne)
        else:
            Exons[exons]=[Ne]



for num in Exons.keys():
    print >>OutFile, num, sum(Exons[num])/len(Exons[num]), sem(Exons[num])
