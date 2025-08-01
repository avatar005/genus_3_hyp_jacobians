`gen_hyp.spyx` is the main computational file. It is written in Cython.

`write_hyp.sage` compiles and makes calls to `gen_hyp.spyx`. Based on the $n$ value and number of blocks (for multithreading) set within this `write_hyp.sage`, hyperelliptics are generated and stored in a directory called `dataq` where `q` is replaced by the size of the field ($2^n$). Note this directory MUST be made before sage can write to it.

`char_2_hyp_point_count.sage` then needs to be run using GNU parallel with the number of blocks from the above.
```bash
ls data32/* | parallel sage char_2_hyp_point_count.sage {}
```
This takes the longest. It computes the point counts for each hyperelliptic and produces all of the characteristic polynomials. 

Finally, `combine_output.sage` combines all of the files produced by the many threads run above into a single file called `all_weil_polys_32`.