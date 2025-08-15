from collections import defaultdict
import sys
import os

n = int(sys.argv[1])
q = 2**n
block = int(sys.argv[2].split('_')[-1])
F.<a> = GF(q)
R.<x> = PolynomialRing(F)
final_hyp = defaultdict(list)
with open(os.path.join("data"+str(2**n), "block_"+str(block)), 'r') as file:
    for line in file:
        polys = line.strip().split(';')
        final_hyp[R(sage_eval(polys[0], locals={'a':a}))].append(R(sage_eval(polys[1], locals={'a':a})))
all_hyp = []
for v in final_hyp:
    for u in final_hyp[v]:
        all_hyp.append((v, u))
out = defaultdict(int)

for i, (v, u) in enumerate(all_hyp):
    C = HyperellipticCurve(u, v)
    frob = C.frobenius_polynomial()
    out[frob] += 1

print(len(out))

out = [str(frob.list())+'\n' for frob in out.keys()]

with open(os.path.join("data"+str(2**n), "out_"+str(block)), 'w') as outfile:
    outfile.writelines(out)

print("done with block " + str(block))