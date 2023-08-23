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

recordlist = SeqIO.parse(sys.argv[1],'fasta')
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
  nY = seq.count('Y')
  nR = seq.count('R')
  nW = seq.count('W')
  nH = seq.count('H')
  nS = seq.count('S')
  n = gcd( nC, nM, nK, nD+nE+1, nY, nR, nW, nH, nS)
  if n==0: continue
  code.append( [ float(nC)/n,
                 float(nM)/n,
                 float(nK)/n,
                 float(nD+nE+1)/n,
                 float(nY)/n,
                 float(nR)/n,
                 float(nW)/n,
                 float(nH)/n,
                 float(nS)/n,
             ] )

N = len(prot)
code = np.asarray(code)
Nc = 9

sigma = float(sys.argv[2])
top1 = 0
top3 = 0
for p in tqdm(range(N), ncols=100):
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
    for r, _, d in save_top:
        if d-last_dist < 1e-6:
            redundancy += 1
        else:
            if last_rank>0:
                for i in range(last_rank, r):
                    map_redundancy[i] = redundancy
            last_rank = r
            last_dist = d
            redundancy = 1

    #print(p, prot[p])
    #print(save_top[0])
    #print(map_redundancy)
    #break
    
    if save_top[0][1] == prot[p] and map_redundancy[1]==1:
        top1 += 1
        top3 += 1
    elif save_top[1][1] == prot[p] and map_redundancy[2]==1:
        top3 += 1
    elif save_top[2][1] == prot[p] and map_redundancy[3]==1:
        top3 += 1

print(N, "top1=", float(top1)/N, "top3=", float(top3)/N)
