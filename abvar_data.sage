from ast import expr
from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict
import pickle
import os
from functools import partial
from termcolor import cprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import os, shutil
attach("isogeny_classes.sage")

RZ = PolynomialRing(ZZ, 'x')

# allows us to store slope-type of Newton Polygon as a single integer
slope_types = {('0A', '0B', '0C', '1A', '1B', '1C'): 0,
               ('1/2A', '1/2B', '1/2C', '1/2D', '1/2E', '1/2F'): 4,
               ('0A', '1/2A', '1/2B', '1/2C', '1/2D', '1A'): 2,
               ('0A', '0B', '1/2A', '1/2B', '1A', '1B'): 3,
               ('1/3A', '1/3B', '1/3C', '2/3A', '2/3B', '2/3C'): 1}

# reference dict between internal slope-type numerical identifier and p-rank
p_rank_dict = {0:3,
               4:0,
               2:1,
               3:2,
               1:0}

# list of primes and dictionary of their respective characteristics
primes = [2,3,4,5,7,8,9,11,13,16,17,19,23,25]
char = {2:2, 3:3, 4:2, 5:5, 7:7, 8:2, 9:3, 11:11, 13:13, 16:2, 17:17, 19:19, 23:23, 25:5}

# types of factoring that can exist for genus 3 Weil polynomials
exp_patterns = ((6,), (2, 4), (2, 2, 2))

def N(g, q):
    '''predicted counts of isogeny classes based on [DiPippo, Howe 2000]'''
    count = q**(g*(g+1)/4)
    count *= 2**g / factorial(g)
    count *= prod([((2*i)/(2*i - 1))**(g + 1 - i) for i in range(1, g + 1)])
    return float(count)

def _prefetch_polys():
    '''Prepare a dictionary of label:polynomial pairs to speed up downloading data from the database'''

    data = defaultdict(list)
    primes = [2,3,4,5,7,8,9,11,13,16,17,19,23,25]
    for field in primes:
        for rec in db.av_fq_isog.search({"q":field}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        print(field)
    return data

def get_data_initial(prefetch=None):
    '''Downloads data from database. Pickles resulting fetched data as a dictionary keyed by hyperelliptic classification; (factoring, slope type); field. 
    Each entry contains polynomial, coefficients in factored form, and LMFDB label'''
    
    def get_abc(factors, classification):
        '''Get values of factored coefficients from factors and factoring type'''
        factors.sort(key=lambda x: len(x))
        
        if classification == (6, ):
            a = factors[0][1]
            b = factors[0][2]
            c = factors[0][3]
            
        elif classification == (2, 4):
            a = factors[0][1]
            b = factors[1][1]
            c = factors[1][2]

        elif classification == (2, 2, 2):
            a = factors[0][1]
            b = factors[1][1]
            c = factors[2][1]
        
        return (a, b, c)

    if os.path.isfile("data.dict"):
        raise Exception("The data already exists. Are you sure you want to re-import. Delete the current 'data.dict' file to re-import.")
    if not prefetch:
        print("prefetch")
        prefetch = _prefetch_polys()
    # querry by hyperelliptic classification; (factoring, slope type); field q: returns the polynomial, coefficients in factored form, and label
    data = defaultdict(partial(defaultdict, partial(defaultdict, list)))
    
    # set up keys for database querry
    fetch_keys = {"g":3}

    # fetch data
    for rec in db.av_fq_isog.search(fetch_keys, ["label", "q", "hyp_count", "poly", "slopes", "simple_factors"]):
        q = rec["q"]
        slope = slope_types[tuple(rec["slopes"])]
        have_hyp = (rec["hyp_count"] > 0)
        factors = rec["simple_factors"]
        for i in range(len(factors)):
            factors[i] = prefetch[factors[i][:-1]]
        factoring_type = tuple([len(factors[i]) - 1 for i in range(len(factors))])
        coeffs = get_abc(factors, factoring_type)
        data[have_hyp][factoring_type, slope][q].append((rec["poly"], coeffs, rec["label"]))

    # store data in "data.dict"
    pickle_data(data)
    
def get_data():
    """Loads saved data file from 'data.dict' if it exists. Otherwise, imports data from LMFDB."""
    try:
        with open("data.dict", 'rb') as filename:
            return(pickle.load(filename))
    except FileNotFoundError:
        get_data_initial()
        return get_data()

def pickle_data(data):
    '''Save data to drive as a pickled dictionary'''
    with open("data.dict", 'wb') as filename:
        pickle.dump(data, filename)
    print("Data saved to 'data.dict' as a pickled dictionary.")

def write_data(data, excluded={}):
    '''Write data as txt files in 'data' folder in working directory. Overwrites all existing files.'''
    # deal with existing files
    try:
        for dir in os.listdir('data'):
            try:
                shutil.rmtree(os.path.join('data', dir))
            except NotADirectoryError:
                pass
    except FileNotFoundError:
        os.mkdir('data')

    # write data as text files
    for have_hyp, D in data.items():
        if have_hyp in excluded:
            continue
        for (factoring, slopes), L2 in D.items():
            dir = os.path.join('data', f'{factoring}_{slopes}')
            if not os.path.isdir(dir):
                os.makedirs(dir)
            for q, polys in L2.items():
                filename = os.path.join(dir, f"data_{have_hyp}_{q}_{factoring}_{slopes}.txt")
                with open(filename, "w") as F:
                    for poly in polys:
                        F.write(f"{poly[1]}, {poly[0]}, {poly[2]} \n")

def our_jacobi_rules(poly, slope, q, factors, coeffs):
    '''Based on input parameters, returns a dictionary containing evaluations of all relevant rules.'''
    # initialize parameters necessary to apply rules.
    p = char[q]
    p_rank = p_rank_dict[slope]
    r_poly = _real_weil(poly, q).coefficients(sparse=False)
    if len(factors) == 1:
        c, b, a = r_poly[:-1]
    elif len(factors) == 3:
        alpha, beta, gamma = coeffs
    else:
        alpha = coeffs[0]
        delta = coeffs[1]
        # easier to convert epsilon from Weil poly -> Real Weil poly than derive from real weil poly
        epsilon = coeffs[2] - 2*q
    
    s,t,u = poly[1:4]

    vq = log(q, p)
        
    # allows use of pre-implemented functions from https://github.com/roed314/abvar-fq/blob/master/isogeny_classes.sage#L3115
    iso_class = IsogenyClass(poly=poly)

    # initialize rule dictionary
    rules = {}

    # compute rules. Formatted descriptions of rules can be found in 
    # "Classifying Weil Polynomials of Jacobians of Genus 3 Hyperelliptic Curves over Finite Fields" by Borodin and May

    if p == 2:
        rules['0.N.N.0'] = (s%2, t%2, u%2) in ((0, 1, 1), (1, 0, 1), (1, 1, 0))
        rules['0.N.0.0'] = slope == 4

        if len(factors) == 2:
            rules['0.2.2.0'] = alpha == 0 and abs(epsilon) == 3

        elif len(factors) == 3:
            if p_rank == 1:
                rules['0.3.1.0'] = (gamma - alpha in (1, sqrt(p*q) or (beta-alpha, gamma-beta) in ((sqrt(p*q), 1), (1, sqrt(p*q)))))
            elif p_rank == 2:
                rules['0.3.2.0'] = alpha + 5 > gamma
                rules['0.3.2.1'] = alpha == -p*vq and gamma > beta and gamma < p*vq - 1
                rules['0.3.2.2'] = gamma == p*vq and beta > alpha and alpha > -p*vq + 1

    elif p != 2:
        rules['1.N.N.0'] = (t%2, u%2) == (0, 1)
        if len(factors) == 1:
            if p_rank == 0:
                rules['1.1.0.0'] = is_prime(q) and q % 3 == 1 and slope == 1 and b > -q and b%q == 0 and (b/q) % 2 == 0
                # rules['1.1.0.0'] = is_prime(q) and q % 3 == 1 and p_rank == 0 and len(exps) == 1 and slopes == 1 and b > 2*q and b%q == 0 and (b/q) % 2 == 1
                rules['1.1.0.1'] = q == p**2 and ((b == -3*q and c%q == 0 and (c/q)%2 == 1) or (b in (-q, -2*q) and abs(a) == 2*p))#or (b > 2*q and c % q == 0 and (c/q) % 2 == 1))

        elif len(factors) == 2:
            if p_rank == 1:
                rules['1.2.1.0'] = is_prime(q) and alpha % 2 == 1 and delta % 4 == 0 and epsilon % 4 == (2 - 2*q) % 4
            elif p_rank == 2:
                rules['1.2.2.0'] = is_prime(q) and delta % 2 == 1 and abs(epsilon) in (2, 3)

        elif len(factors) == 3:
            if p_rank == 0:
                rules['1.3.0.0'] = False
                for k in range(1, vq):
                    rules['1.3.0.0'] |= (beta, gamma) == (3*k, 3*vq)
                    rules['1.3.0.0'] |= (alpha, beta) == (-3*k, -3*vq)
                    rules['1.3.0.0'] |= beta == 0 and (alpha, gamma) in ((-3, 3*vq), (-3*vq, 3))
            elif p_rank == 1:
                rules['1.3.1.0'] = is_square(q) and abs(s) <= q**0.5 and s%2 == 1
                rules['1.3.1.1'] = q % 4 == 3 and q > 3 and is_prime(q) and abs(alpha) <= 3 and abs(gamma) <= 3
            elif p_rank == 2:
                rules['1.3.2.0'] = p != q and ((gamma == vq*p and (beta - alpha in (1, 3))) or (alpha == -vq*p and (gamma - beta in (1, 3))))
                rules['1.3.2.1'] = q % 4 == 1 and (alpha, beta, gamma) in ((-2, -2, 0), (0, 2, 2))
            rules['1.3.N.0'] = alpha**2 + beta**2 + gamma**2 == 9

            forbidden = [(-4, -1, 0), (0, 1, 4), (-4, -3, 0), (0, 3, 4), (-3, -2, 0), (-2, 0, 1), (-1, 0, 2), (0, 2, 3), ]
            rules['1.3.N.1'] = (alpha, beta, gamma) in forbidden
    # resultant 1 method
    rules['N.N.N.1'] = bool(iso_class._nojac_serre())
    # if len(factors) == 2:
    #     rules['N.N.N.1'] = abs(alpha**2 + alpha*delta + epsilon) == 1
    # elif len(factors) == 3:
    #     rules['N.N.N.1'] = (abs(alpha*beta - gamma*(alpha+beta) + gamma**2) == 1 or
    #                         abs(alpha*gamma - beta*(alpha+gamma) + beta**2) == 1 or
    #                         abs(beta*gamma - alpha*(beta + gamma) + alpha**2) == 1)

    rules['N.N.N.0'] = iso_class._nojac_pointcounts()

    rules['N.N.N.2'] = iso_class._nojac_howe_lauter()

    # rules['N.1.N.0'] = len(exps) == 1 and c < 0 and is_prime(-c) and b**3 + c == -q**2

    obstructed_discriminants = [-99, -91, -83, -75, -67, -59, -51, -43, -39, -35, -27, -23, -20, -19, -15, -12, -11, -8, -7, -4, -3, 0]
    rules['N.3.N.0'] = len(factors) == 3 and alpha==beta and beta==gamma and (alpha**2 - 4*q in obstructed_discriminants)

    rules['N.3.0.0'] = p in (2, 3, 5) and len(factors) == 3 and p_rank == 0 and abs(alpha) != p*vq and abs(gamma) != p*vq
    rules['N.3.0.1'] = p_rank == 0 and len(factors) == 3 and ((alpha, beta) == (-p*vq, -p*vq) or (beta, gamma) == (p*vq, p*vq))
    rules['N.2.0.0'] = is_prime(q) and p_rank == 0 and q % 8 != 7 and len(factors) == 2 and (alpha,delta,epsilon) == (0, 0, -4*q)

    return rules


def sort_data(data=None, source=False, dest='undetected'):
    '''Sort data according to rules found in our_jacobi_rules. Original and destination data sets can be set with optional arguments.'''
    if not data:
        data = get_data()
    
    data[dest].clear()
    data["false positives"].clear()
    counts = {dest: defaultdict(int),
              source: defaultdict(int),
              "false positives": defaultdict(int),
              True: defaultdict(int)}
    
    # check all isogeny classes that do not contain a Jacobian and see if they 
    # are classified by our rules. 
    
    for (exps, slope), D in data[source].items():
        print(exps, slope)
        for q, polys in D.items():
            for poly in polys:
                counts[source][(exps, slope)] += 1
                guess = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                
                if not (True in guess.values()):
                    counts[dest][(exps, slope)] += 1
                    data[dest][(exps, slope)][q].append(poly)
    
    print('False done')

    # check for false positives (these occur if a rule claims that an isogeny class
    # with a Jacobian contains a Jacobian)
    for (exps, slope), D in data[True].items():
        print(exps, slope)
        for q, polys in D.items():
            for poly in polys:
                counts[True][(exps, slope)] += 1
                guess = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                if True in guess.values():
                    counts['false positives'][(exps, slope)] += 1
                    data['false positives'][(exps, slope)][q].append(poly)
    
    # print statistics about the sort
    fp = sum(counts['false positives'].values())
    true = sum(counts[True].values())
    print(f"Percent false positive: {fp}/{true} = {float(fp/true)}")
    undetected = sum(counts[dest].values())
    false = sum(counts[source].values())
    print(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}")

    # store data for future use
    pickle_data(data)

def _real_weil(poly, q):
    '''Helper function to compute the real Weil polynomial based on the Weil polynomial'''
    a, b, c = poly[1:4]
    c1, c2, c3 = a, b - 3*q, c - 2*q*a
    return RZ([c3, c2, c1, 1])

def print_statistics(data=None, table=False):
    '''Print coarse statsitics on how many polynomials we have classified in each category. 
    Optional table option to print a LaTeX table with computed data.'''

    # color-coding function
    color = lambda ratio: (255*(ratio)**0.5, 255 - 255*(ratio)**0.5, 0)

    # if no data provided, load data from drive
    if not data:
        data = get_data()

    # initialize variables
    counts = defaultdict(lambda: defaultdict(int))
    headers = ['(Factoring, $p$-rank)', 'Undetected by our rules', 'Undetected by all rules']
    array = defaultdict(list)

    # compute counts to compute statistics
    for have_hyp, D1 in data.items():
        for (exps, slope), D2 in D1.items():
            for q, polys in D2.items():
                counts[have_hyp][(exps, slope)] += len(polys)

    # formatting and printing of counts
    for (exps, slope) in sorted(counts[False], key=lambda x:(x[0], p_rank_dict[x[1]])):
        # code to combine the two p_rank = 0 cases into one output
        if slope == 1:
            continue
        if slope == 4:
            counts[False][(exps, slope)] += counts[False][(exps, 1)]
            counts[False][exps, 1] = 0
            counts['undetected'][(exps, slope)] += counts['undetected'][(exps, 1)]
            counts['undetected'][(exps, 1)] = 0

        # print counts with colored formatting
        print(f"factoring pattern: {exps}, p-rank: {p_rank_dict[slope]}")
        undetected = counts['undetected'][(exps, slope)]
        false = counts[False][(exps, slope)]
        cprint(f"\t Percent negatives undetected by our rules: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
        # store in array for LaTeX table output
        array[exps, slope] = [(len(exps), p_rank_dict[slope]), f'{undetected}/{false} = {float(undetected/false):.4f}']#, f'{unknown}/{false} = {float(unknown/false):.4f}']

    # print cumulative statistics
    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    print()
    cprint(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
    # store cumulative statistics for LaTeX table
    array['last'] = [r'\textbf{Total}', f"{undetected}/{false} = {float(undetected/false):.4f}"]#, f"{unknown}/{false} = {float(unknown/false):.4f}"]

    # print LaTeX table
    if table:
        table = Table([array[l] for l in array], headers)
        print(table)

def print_statistics_by_rule(data=None, table=False):
    '''Print number of hits and unique hits per rule. Optional table option to print a latex table containing the results.'''
    # if data not provided, get data
    if not data:
        data = get_data()
    
    # initialize variables
    counts = defaultdict(set)

    # initialize LaTeX output
    headers = ['Rule', 'Number of hits', 'Number of unique hits']
    array = dict()

    # generate counts
    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule].add(poly[2])
    
    # print results
    for rule in sorted(counts.keys(), key=lambda x:str(x)):
        other_hits = set()
        for r in counts:
            if r != rule:
                other_hits = other_hits.union(counts[r])
        print(f"Rule: {rule}, Number of hits: {len(counts[rule])}, Number of unique hits: {len(counts[rule] - other_hits.intersection(counts[rule]))}")
        array[rule] = [r'\textbf{' + rule + '}', len(counts[rule]), len(counts[rule] - other_hits.intersection(counts[rule]))]

    # print LaTeX table
    if table:
        table = Table([array[l] for l in array], headers)
        print(table)

def print_statistics_by_rule_and_type(data=None, table=False):
    '''Split by polynominal type, print the relevant rules and how many hits (total + unique) they get. 
    Optional table option to print a latex table with the information.'''

    # if no data provided, get data
    if not data:
        data = get_data()

    # intitialize LaTeX output
    headers = ['Type', 'Rule', 'Hits', 'Unique hits']
    table = Table([], headers)

    # comput counts for each output type
    for (exps, slope), D in data[False].items():
        counts = defaultdict(set)
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule].add(poly[2])
        # print results 
        for rule in sorted(counts.keys(), key=lambda x:str(x)):
            other_hits = set()
            for r in counts:
                if r != rule:
                    other_hits = other_hits.union(counts[r])
            unique_hits = len(counts[rule] - other_hits.intersection(counts[rule]))
            print(f"Type: {exps}, {p_rank_dict[slope]}, Rule: {rule}, Number of hits: {len(counts[rule])}, Number of unique hits: {unique_hits}")
            table.add([f"({len(exps)}, {p_rank_dict[slope]})",  rule, len(counts[rule]), unique_hits])
    # print LaTeX tablej
    print(table)

def plot(data, factoring, slope_type, q, factored=True):
    '''Plots coefficients in 3d. Plots two datasets set by positive_key and negative_key. 
    q can either specify the cardinality of the field or special values according to
    q == 0: even
    q == 1: odd
    q == -1: squares
    Plots everything in a single figure with multiple 3d plots. 
    Allows toggles for the two datasets/ Parameter "factored" controls whether the 
    factored or unfactored coefficients are used.'''
    def restriction(x, y, z, q):
        '''Returns boolean for a restriction as a function of x, y, z, q. '''
        xn = x_mod(x, y, z, q)
        yn = y_mod(x, y, z, q)
        zn = z_mod(x, y, z, q)
        return True
    def x_mod(x, y, z, q):
        '''Allows one to apply a modification to the x variable. '''
        return x
    def y_mod(x, y, z, q):
        '''Allows one to apply a modification to the y variable. '''
        return y
    def z_mod(x, y, z, q):
        '''Allows one to apply a modification to the z variable. '''
        return z
    positive_key = True
    negative_key = 'undetected'
    sc1 = []
    sc2 = []
    if q == 0:
        pr, fig, rows, cols = [2, 4, 8, 16], plt.figure(num=1), 2, 2
    elif q == 1:
        pr, fig, rows, cols = [3, 5, 7, 9, 11, 13, 17, 19, 23, 25], plt.figure(num=1), 3, 4
    elif q == -1:
        pr, fig, rows, cols = [4, 9, 16, 25], plt.figure(num=1), 2, 2
    else:
        pr, rows, cols, fig = [q], 1, 1, plt.figure(num=1)
    for i in range(len(pr)):
        q = pr[i]
        
        x1 = []; y1 = []; z1 = []
        x2 = []; y2 = []; z2 = []

        if factored:
            # “negative” group
            for _, (x, y, z), _ in data[negative_key][factoring, slope_type][q]:
                if restriction(x, y, z, q):
                    x1.append(x_mod(x, y, z, q)); y1.append(y_mod(x, y, z, q)); z1.append(z_mod(x, y, z, q))
            # “positive” group
            for _, (x, y, z), _ in data[positive_key][factoring, slope_type][q]:
                if restriction(x, y, z, q):
                    x2.append(x_mod(x, y, z, q)); y2.append(y_mod(x, y, z, q)); z2.append(z_mod(x, y, z, q))
        
        else:
            # “negative” group
            for (_, x, y, z, _, _, _), _, _ in data[negative_key][factoring, slope_type][q]:
                if restriction(x, y, z, q):
                    x1.append(x_mod(x, y, z, q)); y1.append(y_mod(x, y, z, q)); z1.append(z_mod(x, y, z, q))
            # “positive” group
            for (_, x, y, z, _, _, _), _, _ in data[positive_key][factoring, slope_type][q]:
                if restriction(x, y, z, q):
                    x2.append(x_mod(x, y, z, q)); y2.append(y_mod(x, y, z, q)); z2.append(z_mod(x, y, z, q))
                    print(y2[-1]/sqrt(q), q)

        # create the figure & 3D axes
        ax = fig.add_subplot(rows, cols, i + 1, projection='3d')
        ax.set_proj_type('ortho')
        ax.set_xlabel('a')
        ax.set_ylabel('b')
        ax.set_zlabel('c')
        ax.set_title(str(q))

        # plot both groups
        sc1.append(ax.scatter(x1, y1, z1, marker='x', c='red', label=str(negative_key), s=100))
        sc2.append(ax.scatter(x2, y2, z2, c='green', label=str(positive_key), s=100))

    # place checkbuttons in a small inset axes
    # [left, bottom, width, height] in fraction of figure
    ax_checkbox = fig.add_axes([0.02, 0.4, 0.12, 0.15])
    labels     = [str(negative_key), str(positive_key)]
    visibility = [True, True]
    check      = CheckButtons(ax_checkbox, labels, visibility)

    # toggle visibility callback
    def _toggle(label):
        if label == str(negative_key):
            for sc in sc1:
                sc.set_visible(not sc.get_visible())
        else:
            for sc in sc2:
                sc.set_visible(not sc.get_visible())
        plt.draw()

    check.on_clicked(_toggle)

    plt.show()

class Table:
    '''Helper class to convert an array into a latex table.'''
    def __init__(self, array, headers):
        self.table = array
        self.headers = headers
    def __str__(self):
        '''defines printing behavior'''
        string = r'\begin{center}' + '\n'
        string += r'\begin{tabular}{|' + 'c|'*len(self.table[0]) + '}\n'
        string += r'\hline' + '\n'
        for i in range(len(self.headers)):
            string += r'\textbf{' + str(self.headers[i]) + '}'
            if i < len(self.headers) - 1:
                string += r' & '
            else:
                string += r' \\' + '\n\\hline \n'
        for line in self.table:
            for i in range(len(line)):
                string += str(line[i])
                if i < len(line) - 1:
                    string += r' & '
                else:
                    string += r' \\ \hline' + '\n'
        string += r'\end{tabular}' + '\n' + r'\end{center}'
        return string
    def add(self, line):
        '''append a line to an already created Table object.'''
        self.table.append(line)
