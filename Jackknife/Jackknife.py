import matplotlib.pyplot as plt
import sys, random, gzip
InFile=open(sys.argv[1],'r')
GTF=gzip.open(sys.argv[3], 'rt')
TransFile=open(sys.argv[2], 'r')

sample=int(sys.argv[4])
thresh=float(sys.argv[5])
thresh2=float(sys.argv[6])
thresh3=float(sys.argv[7])

TransDict={}
for line in GTF:
  if line[0]!="#":
    A=line.split('\t')
    if A[2]=="transcript":
        info=A[8].split()
        gene=info[1][1:-2]
        trans=info[3][1:-2]
#        print gene, trans
        if gene in TransDict:
            TransDict[gene].append(trans)
        else:
            TransDict[gene]=[trans]

EpTDict={}
TransFile.readline()
for line in TransFile:
    A=line.split(',')
    EpTDict[A[0]]=float(A[1])
        
InFile.readline()
TransFile.readline()
EpG=[]
EpT=[]
Genes=[]
for line in InFile:
    A=line.split(',')
    EpG.append(float(A[3]))
    EpT.append(float(A[1]))
    Genes.append(A[0])

plt.hist(EpG, bins=100)
plt.savefig('EpG_hist.png')
plt.close()

plt.hist(EpT, bins=100)
plt.savefig('EpT_hist.png')
plt.close()

plt.hist(Genes, bins=100)
plt.savefig('Genes_hist.png')
plt.close()

j=float(0.0000)
k=float(0)
l=float(0)
reps=100
means_EpG = []
means_EpT= []
means_TpG = []
for i in range(0,reps):
    foo =random.sample(EpG, sample)
    mean1 = sum(foo)/len(foo)
    means1.append(mean1)
    if sum(foo)/len(foo)> thresh2:
        j=j+1
    foo2=random.sample(EpT, sample)
    mean2 = sum(foo2)/len(foo2)
    means2.append(mean2)
    if sum(foo2)/len(foo2) > thresh:
        k=k+1
    foo3=random.sample(Genes, sample)
    findem=[]
    for gene in Genes:
        for item in TransDict[gene]:
#            print item, gene, EpTDict[item]
            findem.append(EpTDict[item])
    mean3 = sum(findem)/len(findem)
    means3.append(mean3)
    if sum(findem)/len(findem)> thresh3:
        l=l+1
        
print("EpG", j/reps)
print("TpG", k/reps)
print("EpT", l/reps)
plt.hist(means_EpG)
plt.savefig('EpG_jackknife.png')
plt.close()

plt.hist(means_EpT)
plt.savefig('TpG_jackknife.png')
plt.close()

plt.hist(means_TpG)
plt.savefig('EpT_jackknife.png')
plt.close()

