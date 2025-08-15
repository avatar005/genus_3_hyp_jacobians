import os, shutil, sys
gen_hyp = load("gen_hyp_using_gap.sage") 

# set parameters
n = int(sys.argv[1])
v_file = sys.argv[2]
num_blocks = int(sys.argv[3])
process_num = int(sys.argv[4]) - 1
files = []

with open(v_file) as file:
    v_str = file.readline()

# create file handles
for i in range(num_blocks):
    files.append(open(os.path.join("data"+str(2**n), "block_"+str(i + num_blocks*process_num)), 'a+'))

try:
    # write data
    i = 0
    for v, u in generate_hyperelliptics_char2(n, v_str):
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


