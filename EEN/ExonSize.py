from __future__ import division

import gzip, sys

#"GCF_000001405.39_GRCh38.p13_genomic.gtf.gz"
InFile=gzip.open(sys.argv[1], 'r')
out=sys.argv[1][0:-6]
OutFile=open(out+"effectiveIntrons", 'w')

TransDict={}
i=0
for line in InFile:
    if line[0]!="#":
        i=i+1
#        if i>1000:
#            break
        A=line.split()
        if A[2]=="exon":
            chrom=A[0]
            start=int(A[3])
            stop=int(A[4])
            strand=A[6]
            gene=A[9][1:-2]
            trans=A[11][1:-2]
            if TransDict.has_key(trans):
                TransDict[trans].append([chrom, start, stop, strand, gene, trans])
            else:
                TransDict[trans]=[[chrom, start, stop, strand, gene, trans]]


for item in TransDict:
    ExonList=TransDict[item]
    Gene=ExonList[0][4]
    if len(Gene)==0:
        Gene="unknown"
    length=0
    LengthList=[]
    if ExonList[0][3]=="-":
        ExonList.reverse()
    for thing in ExonList:
        stop=thing[2]
        start=thing[1]
        exonlen=stop-start+1
        LengthList.append(exonlen)
        length=length+exonlen
    ne=0
    for ex in LengthList:
        ne=ne+(ex/length)**2
    if len(ExonList)>1:
        print >>OutFile, Gene, item, "{:.2f}".format(1/ne), len(ExonList)#, LengthList
        
    
