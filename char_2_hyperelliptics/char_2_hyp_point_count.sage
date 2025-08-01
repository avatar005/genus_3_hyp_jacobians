from collections import defaultdict
import sys
import os

q = 8
n = 3
block = int(sys.argv[1].split('_')[-1])
F.<a> = GF(q)
R.<x> = PolynomialRing(F)
final_hyp = defaultdict(list)
with open(os.path.join("data"+str(2**n), "block_"+str(block)), 'r') as file:
    for line in file:
        polys = line.strip().split(';')
        # print(polys[1])
        # print(sage_eval(polys[1], locals={'a':a}))
        final_hyp[R(sage_eval(polys[0], locals={'a':a}))].append(R(sage_eval(polys[1], locals={'a':a})))
# print(final_hyp[1])
all_hyp = []
for q in final_hyp:
    for p in final_hyp[q]:
        all_hyp.append((q, p))
# all_hyp = [(q, p) for p in final_hyp[q] for q in final_hyp]
# print(len(all_hyp))
out = defaultdict(int)

for i, (q, p) in enumerate(all_hyp):
    # print(q, p, type(p), type(q))
    # if i % 200 == 0:
    #     print(i)
    C = HyperellipticCurve(p, q)
    frob = C.frobenius_polynomial()
    out[frob] += 1

print(len(out))
# print(out)

out = [str(frob.list())+'\n' for frob in out.keys()]

with open(os.path.join("data_"+str(2**n), "out_"+str(block)), 'w') as outfile:
    outfile.writelines(out)

print("done with block " + str(block))