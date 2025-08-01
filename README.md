<h1> Weil Polynomials of Jacobians of Genus 3 Hyperelliptic Curves</h1>

This repository is a collection of functions to check which isogeny classes of genus 3 abelian varieties contain the Jacobians of hyperelliptic curves. It uses the data found on LMFDB.org. The accompanying paper is [Borodin, May, 2025] (Classifying Weil Polynomials of Jacobians of Genus 3 Hyperelliptic Curves).

This works best in an interactive sage session. In particular, to get started, open the home directory in a live sage terminal and running 

```python
%attach abvar_data.sage
```

The first thing one must do is get the data. This can be done by running

```python
data = get_data()
```

to store the data as a dictionary in the variable `data`. The data is stored as a dictionar keyed by `[presence of hyperelliptic][(factors), slope type][field]`. For more details, see source code. 

To see data as text files, run

```python
write_data(data)
```

This separates data by the dictionary keys described above and stores everything in the directory `data`. 

To sort data based on the rules found in our paper, run

```python
sort_data(data)
```

This takes a while to run and stores the results in `data`. All polynomials corresponding to isogeny classes that do not contain Jacobians of hyperelliptic curves are stored under `data['undetected']`. The data is also automatically stored in `data.dict` so that the sorting computation doesn't have to be run every time one opens a new terminal (so future runs of `get_data()` will load the results of the sorting in addition to the raw data found on the LMFDB). 

The `sort_data()` function uses the rules found in `our_jacobi_rules()`. These rules can be modified to test new conjectured rules. The source and destination of `sort_data()` can also be modified. That is, if instead of sorting all available data, one wants to, say, check how many of the unclassified polynomials are classified by some new rule, one can give the optional arguments `source` and `destination` (to be, say, `'undetected'` and `'new'`, respectively). The rules are still checked against all false positives, so this does not yield a notable speedup, but can be helpful for organization. 

<h2>Statistics</h2>

One of the most useful functions of this codebase is producing statistics regarding which polynomials have and have not been classified.

To get coarse statistics on how many polynomials get classified by the rules, use 

```python
print_statistics(data)
```

This prints out statsitics on number of isogeny classes without hyperelliptic curves classified by `our_jacobi_rules()` divided by the factoring type and $p$-rank.

For statistics divided by rule instead, use

```python
print_statistics_by_rule(data)
```

This produces a list of the rules and the number of polynomials conforming to this rule (`'hits'`). It also produces a `'unique hits'` column, which is the number of polynomials conforming to this rule and no other rules. 

Finally, if one wants to divide the output by factoring type, $p$-rank and rule (a combination of the above functions), use

```python
print_statistics_by_rule_and_type(data)
```

All three functions have an optional argument `table=False`, which, when set to `True`, produces a $\LaTeX$ table with the output in addition to the normal printed output. 

<h2>Plotting</h2>

The final important function of this code is plotting. In particular, one can plot the three parameters which define a genus 3 algebraic variety. These are the $a, b, c$ parameters from our paper, which are defined therein. The

```python
plot(data, factoring, slope_type, q, factored=True):
```

function allows 3-dimensional plots of the data for a given factoring type (given as `(2,2,2)`, `(2,4)` or `(6,)`), $p$-rank (given in terms of slope type &mdash; see the beginnign of the source code for the definitions of `slope_type` in terms of the Newton polygon), and $q$ (size of finite field). The `factored` flag allows one to control whether the Weil polynomial is factored when the three parameters are extracted. There are also special values for $q$ allowing multiple plots to be produced simultaneously &mdash; see the docstring of the function for details. This function is highly customizable from within the source code. 

