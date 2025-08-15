#!/bin/bash

n=3
num_blocks=150
num_cores=10

q=$((2**$n))

# prepare files
rm -rf data$q
mkdir data$q

# generate hyperellitpics
num_vs=$(sage -c "load('gen_hyp_using_gap.sage'); print(write_vs($n))")
num_blocks=$(($num_blocks + ($num_vs - $num_blocks % $num_vs)))
num_blocks_per_process=$(($num_blocks / $num_vs))
echo "$num_vs v's generated"
sage --preparse write_hyp.sage
ls data$q/poly* | parallel -j $num_cores sage -python write_hyp.sage.py $n {} $num_blocks_per_process {#}


# compute point counts (Frobenius polynomials)
sage --preparse char_2_hyp_point_count.sage
ls data$q/block_* | parallel -j $num_cores sage -python char_2_hyp_point_count.sage.py $n {}
sage combine_output_char_2.sage $n

rm *.sage.py

echo "Done"