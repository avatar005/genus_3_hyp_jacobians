#!/bin/bash

n=3
num_blocks=409
num_cores=10
q=$((2**$n))

sage write_hyp.sage $n $num_blocks
sage --preparse char_2_hyp_point_count.sage
ls data$q/* | parallel -j $num_cores python3 char_2_hyp_point_count.sage.py $n {}
sage combine_output_char_2.sage $n $num_blocks

echo "Done"