from sage.misc.cython import cython_import_all
import os

gen_hyp = cython_import_all("gen_hyp.spyx", globals())
import os
n = 2
final_hyp = generate_hyperelliptics_char2(n)
print(sum(len(x) for x in final_hyp.values()))
num_blocks = 4096
files = []
try:
    os.rmdir("data"+str(2**n))
except:
    pass
os.mkdir("data"+str(2**n))
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
