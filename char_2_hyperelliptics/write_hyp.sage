from sage.misc.cython import cython_import_all

gen_hyp = cython_import_all("gen_hyp.spyx", globals())
import os
n = 5
final_hyp = g3_hyp_even_5(n)
print(sum(len(x) for x in final_hyp.values()))
num_blocks = 4096
files = []
for i in range(num_blocks):
    files.append(open(os.path.join("data"+str(2**n), "block_"+str(i)), 'a+'))
i = 0
for q in final_hyp:
    print(q)
    for p in final_hyp[q]:
        files[i%num_blocks].write(str(q.list()) + ';' + str(p.list()) + '\n')
        i += 1

for file in files:
    file.close()
