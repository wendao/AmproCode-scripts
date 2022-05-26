import sys
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from scipy import spatial
from tqdm import tqdm

#greatest common divisor
def gcd(*num):
    num = [x for x in num if x>0]
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

#recordlist = SeqIO.parse('UP000005640_9606.fasta','fasta')
recordlist = SeqIO.parse('secreted_seq.fasta','fasta')
recorddic = SeqIO.to_dict(recordlist)

prot = []
code = []

for gene in recorddic.keys():
  record = recorddic[gene]
  prot.append(record.id)
  seq = record.seq
  nC = seq.count('C')
  nM = seq.count('M')
  nK = seq.count('K')
  nD = seq.count('D')
  nE = seq.count('E')
  nY = seq.count('Y')
  n = gcd( nC, nM, nK+1, nD+nE+1 )
  code.append( [ float(nC)/n,
                 float(nM)/n,
                 float(nK)/n,
                 float(nD+nE+1)/n,
                 float(nY)/n])

N = len(prot)
code = np.asarray(code)

#target = np.array([1.89, 0, 1.0, 0.95, 1.02]) #ACFWKYCV
target = np.array([2.06, 0.01, 1.0, 0.95, 1.00])

dist = []
for i in range(N):
    dist.append(spatial.distance.cosine(code[i,:], target))
dist_ndx = [(d,i) for i,d in enumerate(dist)]
dist_ndx = sorted(dist_ndx, key=lambda x : x[0])

r = 0
for d, i in dist_ndx[:10]:
    r += 1
    print( "rank=", r, prot[i], d)

