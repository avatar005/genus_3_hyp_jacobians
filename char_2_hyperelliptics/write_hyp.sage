import os, shutil, sys
gen_hyp = load("gen_hyp_using_gap.sage") 

# set parameters
n = int(sys.argv[1])

# prepare output location
num_blocks = int(sys.argv[2])
files = []
try:
    shutil.rmtree("data"+str(2**n))
except:
    pass
os.mkdir("data"+str(2**n))

# create file handles
for i in range(num_blocks):
    files.append(open(os.path.join("data"+str(2**n), "block_"+str(i)), 'a+'))

try:
    # write data
    i = 0
    for v, u in generate_hyperelliptics_char2(n):
        files[i%num_blocks].write(str(v.list()) + ';' + str(u.list()) + '\n')
        i += 1
    # close all files
    for file in files:
        file.close()

# if user interrupts, write output anyways
except KeyboardInterrupt:
    # close all files
    for file in files:
        file.close()


