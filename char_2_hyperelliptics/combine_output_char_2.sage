import sys
import os
import ast

n = int(sys.argv[1])
q = 2**n
blocks = int(sys.argv[2])

out = dict()
R.<x> = PolynomialRing(ZZ)

for i in range(blocks):
    with open(os.path.join("data"+str(q), "out_"+str(i)), 'r') as weil_file:
            weil_polys = weil_file.readlines()
            for weil in weil_polys:
                weil_poly = R(ast.literal_eval(weil[:-1]))
                out[weil_poly] = True

with open(os.path.join("data"+str(q), "all_weil_polys_"+str(q)), 'w') as outfile:
    for weil in out:
        outfile.write(str(weil.list()) + '\n')
