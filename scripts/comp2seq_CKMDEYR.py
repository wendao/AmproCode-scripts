import sys
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from scipy import spatial
from tqdm import tqdm

#greatest common divisor
def gcd(*num):
    num = [x for x in num if x>0]
    if len(num)==0: return 0
    gcdl = []
    for i in range(1, sorted(num)[0] + 1):
        for index, j in enumerate(num):
            if j % i == 0:
                if (index + 1) == len(num):
                    gcdl.append(i)
                    break
                continue
            else:
                break
    if not gcdl:
        return 1
    else:
        return sorted(gcdl)[-1]

recordlist = SeqIO.parse('db.fasta','fasta')
recorddic = SeqIO.to_dict(recordlist)

prot = []
code = []

for gene in recorddic.keys():
  record = recorddic[gene]
  seq = record.seq
  nC = seq.count('C')
  nM = seq.count('M')
  nK = seq.count('K')
  nD = seq.count('D')
  nE = seq.count('E')
  nY = seq.count('Y')
  nR = seq.count('R')
  n = gcd( nC, nM, nK, nD+nE+1, nY, nR)
  if n==0: continue
  prot.append(record.id)
  code.append( [ float(nC)/n,
                 float(nM)/n,
                 float(nK)/n,
                 float(nD+nE+1)/n,
                 float(nY)/n,
                 float(nR)/n,
             ] )

N = len(prot)
code = np.asarray(code)
Nc = 6

sigma = float(sys.argv[1])
top1 = 0
top3 = 0
for p in tqdm(range(N),ncols=100):
    #norm, add nosie
    anchor = -1
    min_ct = 9999
    for i in range(Nc):
        if code[p,i]>0 and code[p,i]<min_ct:
            anchor = i
            min_ct = code[p,i]
    target = code[p,:] / min_ct
    #print("norm", prot[p], target)
    target = target * np.random.normal(1.0, sigma, Nc)
    #print("noise", target)

    dist = []
    for i in range(N):
        dist.append(spatial.distance.cosine(code[i,:], target))
    dist_ndx = [(d,i) for i,d in enumerate(dist)]
    dist_ndx = sorted(dist_ndx, key=lambda x : x[0])
    
    #check rank
    if dist_ndx[0][1] == p:
        top1 += 1
        top3 += 1
    elif dist_ndx[1][1] == p or dist_ndx[2][1] == p:
        top3 += 1

print(N, "top1=", float(top1)/N, "top3=", float(top3)/N)
