import sys
from math import sqrt
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

recordlist = SeqIO.parse(sys.argv[1], 'fasta')
recorddic = SeqIO.to_dict(recordlist)

prot = []
code = []

for gene in recorddic.keys():
  record = recorddic[gene]
  prot.append("|".join(record.id.split('|')[:-1]))
  seq = record.seq
  nC = seq.count('C')
  nM = seq.count('M')
  nK = seq.count('K')
  nD = seq.count('D')
  nE = seq.count('E')
  n = gcd( nC, nM, nK, nD+nE+1 )
  code.append( [ float(nC)/n,
                 float(nM)/n,
                 float(nK)/n,
                 float(nD+nE+1)/n ])

N = len(prot)
code = np.asarray(code)

target = [ float(i) for i in sys.argv[2:] ]

dist = []
for i in range(N):
    dist.append(spatial.distance.cosine(code[i,:], target))
dist_ndx = [(d,i) for i,d in enumerate(dist)]
dist_ndx = sorted(dist_ndx, key=lambda x : x[0])

r = 0
check_id = []
save_top = []
for d, i in dist_ndx:
    if prot[i] in check_id:
        continue
    else:
        r += 1
        check_id.append(prot[i])
        save_top.append( (r, prot[i], d) )
        if r==10: break
save_top.append( (11, "FAKE", 999) ) #close

map_redundancy = {}
last_dist = -1.0
last_rank = 0
redundancy = 1
for r, p, d in save_top:
    #print("#", r, p, d )
    if d-last_dist < 1e-6:
        #same level
        #print("same level")
        redundancy += 1
    else:
        #new level
        #print("new level")
        if last_rank>0:
            #print("check:", last_rank, r)
            for i in range(last_rank, r):
                map_redundancy[i] = redundancy
        last_rank = r
        last_dist = d
        redundancy = 1

for rank, redundancy in map_redundancy.items():
    print("rank=", rank, save_top[rank-1][1], save_top[rank-1][2], redundancy)
