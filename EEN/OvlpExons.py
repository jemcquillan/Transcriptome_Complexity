from __future__ import division

import gzip, sys

#"GCF_000001405.39_GRCh38.p13_genomic.gtf.gz"
InFile=gzip.open(sys.argv[1], 'r')
out=sys.argv[1][0:-6]
OutFile=open(out+"ovlpExons", 'w')

TransDict={}
GeneDict={}
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
            trans=A[13][1:-2]
            if TransDict.has_key(trans):
                TransDict[trans].append([chrom, start, stop, strand, gene, trans])
            else:
                TransDict[trans]=[[chrom, start, stop, strand, gene, trans]]

GeneDict={}
for item in TransDict:
    Exons=TransDict[item]
    chrom=Exons[0][0]
    start=int(Exons[0][1])
    stop=int(Exons[-1][2])
    gene=Exons[0][4]
    stramd=Exons[0][3]
    if not gene in GeneDict:
        GeneDict[gene]=[chrom, start, stop, strand]

print "CheckGeneOvlp", len(GeneDict.keys()), "genes"
Overlp=[]
Pairs={}
for gene1 in GeneDict.keys():
    for gene2 in GeneDict.keys():
        if gene1!=gene2:
            exon1=GeneDict[gene1]
            exon2=GeneDict[gene2]
            chrom1=exon1[0]
            chrom2=exon2[0]
            if chrom1==chrom2:
                strand1=exon1[3]
                strand2=exon2[3]
                if strand1==strand2:
                    dummy=0
                    start1=int(exon1[1])
	            stop1=int(exon1[2])
	            start2=int(exon2[1])
                    stop2=int(exon2[2])
                    if start1 <= start2 and start2<=stop1:
                        dummy=1
                    elif start1 <= stop2 and stop2 <= stop1:
                        dummy=1
                    elif start2 <= stop1 and stop1 <=stop2:
                        dummy=1
                    if dummy==1:
                        Pairs[gene1]=gene2
                        Pairs[gene2]=gene1
                        if not gene1 in Overlp:
                            Overlp.append(gene1)
                        if not gene2 in Overlp:
                            Overlp.append(gene2)
                        break

print len(Overlp), "overlapping genes out of ", len(GeneDict.keys())                

TransDictOvlp={}
for item in TransDict.keys():
    ExonList1=TransDict[item]
    Gene1=ExonList1[0][4]
    if Gene1 in Overlp:
        TransDictOvlp[item]=TransDict[item]
print "DictDownsized to ", len(TransDictOvlp.keys()), "out of " , len(TransDict.keys())
Final=[]
for item in TransDictOvlp.keys():
    dummy=0
    for thingy in TransDictOvlp.keys():
        if item!=thingy:
            ExonList1=TransDict[item]
            Gene1=ExonList1[0][4]
            if len(Gene1)==0:
                Gene1="unknown"
            ExonList2=TransDict[thingy]
	    Gene2=ExonList2[0][4]
            if len(Gene2)==0:
                Gene2="unknown"
            if Gene2==Pairs[Gene1] and Gene1 in Overlp and Gene2 in Overlp and ExonList1[0][0]==ExonList2[0][0] and ExonList1[0][3]==ExonList2[0][3] and Gene1!=Gene2:
                for exon1 in ExonList1:
                    for exon2 in ExonList2:
                        start1=int(exon1[1]) 
                        stop1=int(exon1[2])
                        start2=int(exon2[1])
                        stop2=int(exon2[2])
                        if start1 <=start2 and start2<=stop1:
                            dummy=1
                        elif start1<= stop2	and stop2 <= stop1:
                            dummy=1
                        elif start2	<= stop1 and stop1 <=stop2:
                            dummy=1
                            if dummy==1:
                                break
                    if dummy==1:
                        break
                if dummy==1:
                    print >>OutFile, item, thingy
                    if not item in Final:
                        Final.append(item)
                    if not thingy in Final:
                        Final.append(thingy)
                    break

print len(Final), "overlapping exons out of " , len(TransDict.keys())
