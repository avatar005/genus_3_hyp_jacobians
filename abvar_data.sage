from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict
import itertools
import pickle
import os
from functools import partial
from termcolor import cprint
from pprint import pprint
from itertools import product, chain, combinations
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons


slope_types = {('0A', '0B', '0C', '1A', '1B', '1C'): 0,
               ('1/2A', '1/2B', '1/2C', '1/2D', '1/2E', '1/2F'): 4,
               ('0A', '1/2A', '1/2B', '1/2C', '1/2D', '1A'): 2,
               ('0A', '0B', '1/2A', '1/2B', '1A', '1B'): 3,
               ('1/3A', '1/3B', '1/3C', '2/3A', '2/3B', '2/3C'): 1}

p_rank_dict = {0:3,
               4:0,
               2:1,
               3:2,
               1:0}

modulo = 4

primes = [2,3,4,5,7,8,9,11,13,16,17,19,23,25]
char = {2:2, 3:3, 4:2, 5:5, 7:7, 8:2, 9:3, 11:11, 13:13, 16:2, 17:17, 19:19, 23:23, 25:5}

def prefetch_polys():
    data = defaultdict(list)
    primes = [2,3,4,5,7,8,9,11,13,16,17,19,23,25]
    for field in primes:
        for rec in db.av_fq_isog.search({"q":field}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        print(field)
    return data

def pickle_data(data):
    with open("data.dict", 'wb') as filename:
        pickle.dump(data, filename)

def get_data_initial(prefetch=None):
    if os.path.isfile("data.dict"):
        raise Exception("The data already exists. Are you sure you want to re-import")
    if not prefetch:
        print("prefetch")
        prefetch = prefetch_polys()
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

    pickle_data(data)
    
def get_data():
    with open("data.dict", 'rb') as filename:
        return(pickle.load(filename))

def write_data(data):
    for have_hyp, D in data.items():
        for (factoring, slopes), L2 in D.items():
            for q, polys in L2.items():
                filename = f"new/data_{have_hyp}_{q}_{factoring}_{slopes}.txt"
                with open(filename, "w") as F:
                    for poly in polys:
                        F.write(f"{poly[1]}, {poly[0]}, {poly[2]}, {poly[1][0]**2 - poly[1][1]**2  + poly[1][2]} \n")

def all_coeff_patterns(mod):
    for i in range(mod):
        for j in range(mod):
            for k in range(mod):
                yield (i, j, k)

def restriction(x, y, z, q):
    return True
def x_modification(x_value):
    return x_value
def y_modification(y_value):
    return y_value
def z_modification(z_value):
    return z_value


def plot(data, factoring, slope_type, q):
    positive_key = True
    negative_key = 'unknown'
    sc1 = []
    sc2 = []
    if q == 0:
        pr = [2, 4, 8, 16]
        fig = plt.figure(num=1)
        rows = 2
        cols = 2

    elif q == 1:
        pr = [3, 5, 7, 9, 11, 13, 17, 19, 23, 25]
        fig = plt.figure(num=1)
        rows = 3
        cols = 4
    elif q == -1:
        pr = [9, 25]
        fig = plt.figure(num=1)
        rows = 1
        cols = 2
    else:
        pr = [q]
        rows = 1
        cols = 1
        fig = plt.figure(num=1)
    for i in range(len(pr)):
        q = pr[i]
        
        x1 = []; y1 = []; z1 = []
        x2 = []; y2 = []; z2 = []

        # “negative” group
        for _, (x, y, z), _ in data[negative_key][factoring, slope_type][q]:
            if restriction(x, y, z, q):
                x1.append(x); y1.append(y); z1.append(z)
        # “positive” group
        for _, (x, y, z), _ in data[positive_key][factoring, slope_type][q]:
            if restriction(x, y, z, q):
                x2.append(x); y2.append(y); z2.append(z)

        # create the figure & 3D axes
        # ax = fig.add_subplot(111, projection='3d')
        ax = fig.add_subplot(rows, cols, i + 1, projection='3d')
        ax.set_proj_type('ortho')
        ax.set_xlabel('a')
        ax.set_ylabel('b')
        ax.set_zlabel('c')
        ax.set_title(str(q))

        # plot both groups
        sc1.append(ax.scatter(x1, y1, z1, marker='x', c='red', label=str(negative_key), s=50))
        sc2.append(ax.scatter(x2, y2, z2, c='green', label=str(positive_key), s=50))

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

    

exp_patterns = ((6,), (2, 4), (2, 2, 2))
def pattern_keys():
    for i in exp_patterns:
        for j in range(5):
            yield (i, j, True)

def n_tuple(length, lower, upper):
    tup = [range(lower, upper + 1) for i in range(length)]
    return itertools.product(*tup)

def evaluate(coeffs, input):
    dic = {0:0, 1:1, 2:-1, 3:4, 4:-4}
    a, b, c, d, e, f= coeffs
    a = dic[a]
    b = dic[b]
    c = dic[c]
    d = dic[d]
    e = dic[e]
    f = dic[f]
    x1, x2, x3 = input
    return x1**2 * a + x2**2 * b + x3**2 * c + d*x1 + e*x2 + f*x3


def poly_rule_check(data):
    # for mod in range(2, 20):
    #     print(mod)
    p = 2
    for pattern in n_tuple(6, 0, 2):
        if pattern[2] != p:
            print(pattern)
            p = pattern[1]
        for (exps, slopes), D1 in data['norule'].items():
            flag_no = defaultdict(lambda: defaultdict(int))
            for q, polys in D1.items():
                for poly in polys:
                    coeffs = poly[1]
                    if flag_no[q][evaluate(pattern, coeffs)] == 0:
                        flag_no[q]['remaining'] += 1
                    flag_no[q][evaluate(pattern, coeffs)] += 1
                for key in flag_no[q]:
                    if flag_no[q][key] < 4 and key != 'remaining':
                        flag_no[q][key] = 0
                        flag_no[q]['remaining'] -= 1
                for poly in data[True][(exps, slopes)][q]:
                    coeffs = poly[1]
                    if flag_no[q][evaluate(pattern, coeffs)] != 0:
                        flag_no[q][evaluate(pattern, coeffs)] = 0
                        flag_no[q]['remaining'] -= 1
                    if flag_no[q]['remaining'] == 0:
                        break
            count = [0,0]
            for q in flag_no.keys():
                if flag_no[q]['remaining'] > 0:
                    count[q % 2] += 1
            if count[0] >= 3:
                print(pattern, (exps, slopes), end=': ')
                for q in [2, 4, 8, 16]:
                    for c in flag_no[q]:
                        if flag_no[q][c] > 0 and c != 'remaining':
                            print(q, c, flag_no[q][c], end='; ')
                print()
            if count[1] > 4:
                print(pattern, (exps, slopes), end=': ')
                for q in [3, 5, 7, 9, 11, 13, 17, 19, 23, 25]:
                    for c in flag_no[q]:
                        if flag_no[q][c] > 0 and c != 'remaining':
                            print(q, c, flag_no[q][c], end='; ')
                print()

EXP_PATTERNS = {(6,)}
SLOPES      = {0,1,2,3,4}
ALL_RULES   = [
    (eps1,eps2,eps3,eps4, del1,del2,del3)
    for eps1,eps2,eps3,eps4 in product((1,2,3), repeat=4)
    for del1,del2,del3 in product((-1,0,1), repeat=3)
]

def rule_check(data):
    # print(f"R = {r}")
    alive = set(ALL_RULES)
    obs   = {r: [] for r in ALL_RULES}
    unobs = {r: [] for r in ALL_RULES}

    for q in set(primes)-{2,4,8,16}:
        print(q)
        p = char[q]
        for slope in SLOPES:
            for exp_pat in EXP_PATTERNS:
                for have_hyp in (True, False):
                    for poly in data[have_hyp][(exp_pat, slope)][q]:
                        a,b,c = poly[1]
                        if c < 0:
                            # div = Integer(c).prime_divisors()
                            # div_ = [p for p in div if p < abs(c)]
                            # if all([p%12==1 for p in div_]) and b**3+c==-q**2:
                            if is_prime(abs(c)):
                                a_p = {e: a**e for e in (1,2,3)}
                                b_p = {e: b**e for e in (1,2,3)}
                                c_p = {e: c**e for e in (1,2,3)}
                                q_p = {e: q**e for e in (1,2,3)}
                                for rule in list(alive):
                                    e1,e2,e3,e4, d1,d2,d3 = rule
                                    if d1*a_p[e1] + d2*b_p[e2] + d3*c_p[e3] == q_p[e4]:
                                        if have_hyp:
                                            unobs[rule].append((a,b,c))
                                            alive.discard(rule)
                                        else:
                                            obs[rule].append((a,b,c))
    results = [
        (rule, len(obs[rule]))
        for rule in alive
        if obs[rule] and not unobs[rule]
    ]
    return results


def power_rules(results, data):
    for q in {2,4,8,16}:
        print(q)
        for slope in SLOPES:
            for exp_pat in EXP_PATTERNS:
                for have_hyp in (True, False):
                    for poly in data[have_hyp][(exp_pat, slope)][q]:
                        a,b,c = poly[1]
                        if c<0 and is_prime(abs(c)):
                            rules_good = []
                            for e1,e2,e3,e4,d1,d2,d3 in all_rules_odd:
                                if d1*a**e1+d2*b**e2+d3*c**e3 == q**e4:
                                    # rules_good.append(rules[(e1,e2,e3,e4,d1,d2,d3)])
                                    rules_good.append((e1,e2,e3,e4,d1,d2,d3))
                            if rules_good:
                                print(f'poly: {(a,b,c)}, rules: {rules_good}')

def point_count(data):
    obstructed, unobstructed = {}, {}
    for q, D in data.items():
        obstructed[q], unobstructed[q] = {"neg": [], "div": [], "over": []}, []
        if q<=25:
            for _, polys in D.items():
                for i in range(len(polys)):
                    iso_class = IsogenyClass(label=polys[i][2])
                    obs, type_, params = _nojac_pointcounts(iso_class)
                    if obs:
                        obstructed[q][type_].append((iso_class, params))
                        # print(iso_class, "obs")
                    else:
                        unobstructed[q].append(iso_class)
                        # print(iso_class, "unobs")
            print(q, len(obstructed[q]["neg"]), len(obstructed[q]["div"]), len(unobstructed[q]))
    return obstructed, unobstructed

def powerset(iterable):
    """
    Return the powerset of the given iterable as an iterator of tuples.
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def boolean_powerset_funcs(funcs_dict):
    """
    Given a dict mapping label->boolean function,
    return a dict mapping each subset-label to a new function that ANDs
    all funcs in that subset.

    - Empty subset gets label 'TRUE' and always returns True.
    """
    result = {}
    items = list(funcs_dict.items())
    for subset in powerset(items):
        if subset:
            labels, funcs = zip(*subset)
            combo_label = ' and '.join(labels)
        else:
            combo_label = 'TRUE'
            funcs = ()
        def combo_fn(*args, _funcs=funcs, **kwargs):
            return all(f(*args, **kwargs) for f in _funcs)

        result[combo_label] = combo_fn

    return result

def power_rule_check(rule, coeffs, q):
    e1, e2, e3, e4, d1, d2, d3 = rule
    a, b, c = coeffs
    return d1*a**e1 + d2*b**e2 + d3*c**e3 == q**e4

HYPOTHESES = {'a prime': lambda a,b,c: is_prime(abs(a)),
              'b prime': lambda a,b,c: is_prime(abs(b)),
              'c prime': lambda a,b,c: is_prime(abs(c)),
              'a<0': lambda a,b,c: a<0,
              'b<0': lambda a,b,c: b<0,
              'c<0': lambda a,b,c: c<0}

ALL_HYPOTHESES = boolean_powerset_funcs(HYPOTHESES)

def our_jacobi_rules(poly, slopes, q, exps, coeffs):
    p = trial_division(q)

    coeffs4 = tuple([c % 4 for c in coeffs])
    poly2 = tuple([c % 2 for c in poly[1:4]])
    sorted_coeffs = tuple(sorted(coeffs))
    if len(exps) == 3:
        a,b,c = sorted_coeffs
    else:
        a, b, c = coeffs

    p_rank = p_rank_dict[slopes]
    
    rules = {}

    vq = log(q, p)

    if p == 2:
        # Case 1 of Theorem 1
        rules['0.N.N.0'] = poly2 in ((0, 1, 1), (1, 0, 1), (1, 1, 0))

        # Prop 1
        rules['0.N.0.0'] = slopes == 4

        rules['0.2.2.0'] = p_rank == 2 and len(exps) == 2 and a == 0 and c in (2*q - 3, 2*q + 3)
        # rules['0.3.3.0'] = p_rank == 3 and len(exps) == 3 and a == b and b == c and a%2 == 1 and abs(a) > ceil(q/4)

        rules['0.3.1.0'] = p_rank == 1 and len(exps) == 3 and (c - a in (1, sqrt(p*q) or (b-a, c-b) in ((sqrt(p*q), 1), (1, sqrt(p*q)))))

    elif p != 2:
        # Costa Theorem 2.8
        rules['1.N.N.0'] = poly2[1:] == (0, 1)

        # 2x4 splitting, slope type 2, prime field != 2, (1, 0, 2) mod 4 or (3, 0, 2) mod 4 pattern
        rules['1.2.1.0'] = is_prime(q) and p_rank == 1 and len(exps) == 2 and coeffs4 in ((1, 0, 2), (3, 0, 2))

        rules['1.2.2.0'] = is_prime(q) and p_rank == 2 and len(exps) == 2 and b % 2 == 1 and c in (2*q - 2, 2*q - 3, 2*q + 2, 2*q + 3)

        rules['1.3.0.0'] = False
        if p == 3 and p_rank == 0 and len(exps) == 3:
            for k in range(1, vq):
                rules['1.3.0.0'] |= (b, c) == (3*k, 3*vq)
                rules['1.3.0.0'] |= (a, b) == (-3*k, -3*vq)
                rules['1.3.0.0'] |= b == 0 and (a, c) in ((-3, 3*vq), (-3*vq, 3))

        # rules['1.3.1.0'] = p_rank == 1 and len(exps) == 3 and (sorted_coeffs == (-3, 0, 0) or sorted_coeffs == (0, 0, 3))

        rules['1.3.1.1'] = is_square(q) and p_rank == 1 and len(exps) == 3 and abs(poly[1]) <= q**0.5 and poly[1]%2 == 1

        # rules['1.3.1.2'] = q % 4 == 3 and is_prime(q) and p_rank == 1 and len(exps) == 3 and sorted_coeffs in ((-2, 0, 0), (0, 0, 2))
        rules['1.3.1.2'] = q % 4 == 3 and q > 3 and is_prime(q) and p_rank == 1 and len(exps) == 3 and abs(a) <= 3 and abs(c) <= 3

        rules['1.3.1.3'] = q % 4 == 1 and is_prime(q) and p_rank == 1 and len(exps) == 3 and (abs(a) in (1, 3) or abs(c) in (1, 3))

        forbidden = [(-4, -1, 0), (0, 1, 4), (-4, -3, 0), (0, 3, 4), (-3, -2, 0), (-2, 0, 1), (-1, 0, 2), (0, 2, 3), ]
        rules['1.3.2.0'] = p_rank == 2 and len(exps) == 3 and (a, b, c) in forbidden

        rules['1.3.2.1'] = q % 4 == 1 and p_rank == 2 and len(exps) == 3 and (a, b, c) in ((-2, -2, 0), (0, 2, 2))

        # rules['1.3.3.0'] = p_rank == 3 and len(exps) == 3 and a**2 == b**2 and b**2 == c**2 and a%2 == 1

        rules['1.3.3.1'] = p_rank == 3 and len(exps) == 3 and a == b and b == c and abs(a) > ceil(q/3)

        rules['1.3.3.2'] = p_rank == 3 and len(exps) == 3 and a**2 + b**2 + c**2 == 9
        
    rules['N.1.N.0'] = len(exps) == 1 and c < 0 and is_prime(-c) and b**3 + c == -q**2

    # 2x2x2 splitting of the form (x**2 - tx + q)**3, conditions on discriminant t**2 - 4q
    obstructed_discriminants = [-99, -91, -83, -75, -67, -59, -51, -43, -39, -35, -27, -23, -20, -19, -15, -12, -11, -8, -7, -4, -3, 0]
    rules['N.3.N.0'] = len(exps) == 3 and a==b and b==c and (a**2 - 4*q in obstructed_discriminants)

    # 2x2x2 splitting, if two are equal and third is not (so abs(a - b) + abs(b - c) + abs(a - c) == 2)
    rules['N.3.N.1'] = len(exps) == 3 and abs(a - b) + abs(b - c) + abs(a - c) == 2

    # 2x2x2 splitting with parameters a >= b >= c, checks that a - b = b - c = 1
    rules['N.3.N.2'] = len(exps) == 3 and c - b == 1 and b - a == 1

    rules['N.3.0.1'] = p_rank == 0 and len(exps) == 3 and ((a, b) == (-p*vq, -p*vq) or (b, c) == (p*vq, p*vq))

    rules['S.2.0.0'] = is_prime(q) and q % 8 != 7 and p_rank == 0 and len(exps) == 2 and (a,b,c) == (0, 0, -2*q)

    rules['S.3.0.0'] = p in (2, 3, 5) and len(exps) == 3 and p_rank == 0 and abs(a) != p*vq and abs(c) != p*vq

    # 2x2x2 splitting, p-rank=0, for q=2,3,4,9 if a**2 = b**2 = c**2 and p | a**2
    # rules['S.3.0.1'] = p in {2, 3, 5} and p_rank == 0 and len(exps) == 3 and a**2 == b**2 and b**2 == c**2
    
    return rules

def our_jacobi_rules_old(poly, slopes, q, exps, coeffs):
    p = trial_division(q)

    coeffs4 = tuple([c % 4 for c in coeffs])
    poly2 = tuple([c % 2 for c in poly[1:4]])
    sorted_coeffs = tuple(sorted(coeffs))
    if len(exps) == 3:
        a,b,c = sorted_coeffs
    else:
        a, b, c = coeffs

    prime_nums = Primes()
    
    rules = {}

    # Costa Theorem 2.8
    rules['Costa'] = p%2 != 0 and poly2[1:] == (0, 1)

    # Case 1 of Theorem 1
    rules['Thm 1.1'] = p == 2 and poly2 == (0, 1, 1) 

    # Case 2 in Theorem 1
    rules['Thm 1.2'] = p == 2 and slopes != 3 and poly2 == (1, 0, 1)

    # Case 3 in Theorem 1
    rules['Thm 1.3'] = p == 2 and slopes == 3 and poly2 == (1, 1, 0)

    # Prop 1
    rules['Proposition 1'] = p == 2 and slopes == 4

    # 2x4 splitting, slope type 2, prime field != 2, (1, 0, 2) mod 4 or (3, 0, 2) mod 4 pattern
    rules[2] = (p != 2 and is_prime(q)) and slopes == 2 and len(exps) == 2 and (coeffs4 == (1, 0, 2) or coeffs4 == (3, 0, 2))

    # 2x2x2 splitting of the form (x**2 - tx + q)**3, conditions on discriminant t**2 - 4q
    obstructed_discriminants = [-99, -91, -83, -75, -67, -59, -51, -43, -39, -35, -27, -23, -20, -19, -15, -12, -11, -8, -7, -4, -3, 0]
    rules[3] = len(exps) == 3 and a==b==c and (a**2 - 4*q in obstructed_discriminants)
    
    # 2x2x2 splitting with parameters a >= b >= c, checks that a - b = b - c = 1
    rules[4] = len(exps) == 3 and c - b == 1 and b - a == 1
    
    # 2x2x2 splitting, if two are equal and third is not (so abs(a - b) + abs(b - c) + abs(a - c) == 2)
    rules[5] = len(exps) == 3 and abs(a - b) + abs(b - c) + abs(a - c) == 2

    # 2x2x2 splitting, p-rank=0, for q=2,3,4,9 if a**2 = b**2 = c**2 and p | a**2
    rules[6] = p in {2, 3, 5} and slopes in {1, 4} and len(exps) == 3 and a**2 == b**2 and b**2 == c**2

    forbidden = [(-4, -1, 0), (0, 1, 4), (-4, -3, 0), (0, 3, 4), (-3, -2, 0), (-2, 0, 1), (-1, 0, 2), (0, 2, 3), ]
    
    rules[7] = p != 2 and is_prime(q) and slopes == 2 and len(exps) == 3 and (sorted_coeffs == (-3, 0, 0) or sorted_coeffs == (0, 0, 3))
    rules[8] = q % 4 == 3 and is_prime(q) and slopes == 2 and len(exps) == 3 and (sorted_coeffs == (-2, 0, 0) or sorted_coeffs == (0, 0, 2))
    rules[9] = p != 2 and is_square(q) and slopes == 2 and len(exps) == 3 and (abs(poly[1]) <= q**0.5 and poly[1]%2 == 1)
    rules[10] = p == q and slopes == 3 and len(exps) == 2 and b % 2 == 1 and c in (2*q - 2, 2*q - 3, 2*q + 2, 2*q + 3)
    rules[11] = p == 2 and slopes == 3 and len(exps) == 2 and a == 0 and c in (2*q - 3, 2*q + 3)
    rules[12] = p == 2 and slopes == 3 and len(exps) == 3 and a**2 + c**2 + a - c == 6
    rules[13] = p != 2 and slopes == 0 and len(exps) == 3 and a**2 + b**2 + c**2 == 9
    rules[14] = p != 2 and slopes == 3 and len(exps) == 3 and (a, b, c) in forbidden
    rules[15] = q % 4 == 1 and slopes == 3 and len(exps) == 3 and (a, b, c) in [(-2, -2, 0), (0, 2, 2)]

    rules[16] = len(exps) == 1 and c < 0 and -c in prime_nums and b + c == -q**2
    rules[17] = len(exps) == 1 and c < 0 and -c in prime_nums and b**3 + c == -q**2
    
    # if iso_class.has_jacobian == -1 or rule1 or rule2 or rule3:
    return rules

def get_abc(factors, classification):
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

def sort_data(data=None):
    if not data:
        data = get_data()
    
    data["undetected"].clear()
    data["false positives"].clear()
    counts = {"undetected": defaultdict(int),
              False: defaultdict(int),
              "false positives": defaultdict(int),
              True: defaultdict(int)}
    
    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                counts[False][(exps, slope)] += 1
                guess = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                if not (True in guess.values()):
                    counts['undetected'][(exps, slope)] += 1
                    data['undetected'][(exps, slope)][q].append(poly)
    
    for (exps, slope), D in data[True].items():
        for q, polys in D.items():
            for poly in polys:
                counts[True][(exps, slope)] += 1
                guess = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                if True in guess.values():
                    counts['false positives'][(exps, slope)] += 1
                    data['false positives'][(exps, slope)][q].append(poly)
    
    fp = sum(counts['false positives'].values())
    true = sum(counts[True].values())
    print(f"Percent false positive: {fp}/{true} = {float(fp/true)}")
    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    print(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}")

    pickle_data(data)

    # return data


def test_new_rule(data=None):
    if not data:
        data = get_data()
    
    data["norule"].clear()
    data["false positives"].clear()
    counts = {"norule": defaultdict(int),
              False: defaultdict(int),
              "false positives": defaultdict(int),
              True: defaultdict(int)}
    
    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                counts[False][(exps, slope)] += 1
                # test rule !!!!!!!!!!!!!!!!!!!!
                # evaluates to True if we don't expect a jacobian
                p = trial_division(q)
                a, b, c = poly[1]
                polyn = tuple([i%2 for i in poly[0][1:4]])
                newrule = [
                    polyn == (1,0, 1) and p == 2 and slope == 0
                ]
                if not (True in newrule):
                    counts['norule'][(exps, slope)] += 1
                    data['norule'][(exps, slope)][q].append(poly)
    
    for (exps, slope), D in data[True].items():
        for q, polys in D.items():
            for poly in polys:
                counts[True][(exps, slope)] += 1
                # test rule !!!!!!!!!!!!!!!!!!!!
                # evaluates to True if we don't expect a jacobian
                p = trial_division(q)
                a, b, c = poly[1]
                polyn = tuple([i%2 for i in poly[0][1:4]])
                newrule = [
                    polyn == (1,0,1) and p == 2 and slope == 0
                ]
                if True in newrule:
                    counts['false positives'][(exps, slope)] += 1
                    data['false positives'][(exps, slope)][q].append(poly)
    
    fp = sum(counts['false positives'].values())
    true = sum(counts[True].values())
    print(f"Percent false positive: {fp}/{true} = {float(fp/true)}")
    undetected = sum(counts['norule'].values())
    false = sum(counts[False].values())
    print(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}")

    pickle_data(data)

    # return data


class Table:
    def __init__(self, array, headers):
        self.table = array
        self.headers = headers
    def __str__(self):
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
color = lambda ratio: (255*(ratio)**0.5, 255 - 255*(ratio)**0.5, 0)

def print_statistics(data=None, table=False):
    if not data:
        data = get_data()

    counts = defaultdict(lambda: defaultdict(int))
    headers = ['(Factoring, $p$-rank)', 'Undetected by our rules', 'Undetected by all rules']
    array = defaultdict(list)

    for have_hyp, D1 in data.items():
        for (exps, slope), D2 in D1.items():
            for q, polys in D2.items():
                counts[have_hyp][(exps, slope)] += len(polys)

    for (exps, slope) in sorted(counts[False], key=lambda x:(x[0], p_rank_dict[x[1]])):
        if slope == 1:
            continue
        if slope == 4:
            counts[False][(exps, slope)] += counts[False][(exps, 1)]
            counts[False][exps, 1] = 0
            counts['undetected'][(exps, slope)] += counts['undetected'][(exps, 1)]
            counts['undetected'][(exps, 1)] = 0
            counts['unknown'][(exps, slope)] += counts['unknown'][(exps, 1)]
            counts['unknown'][(exps, 1)] = 0
            counts['norule'][(exps, slope)] += counts['norule'][(exps, 1)]
            counts['norule'][(exps, 1)] = 0
        print(f"factoring pattern: {exps}, p-rank: {p_rank_dict[slope]}")
        undetected = counts['undetected'][(exps, slope)]
        false = counts[False][(exps, slope)]
        cprint(f"\t Percent negatives undetected by our rules: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
        unknown = counts['unknown'][(exps, slope)]
        cprint(f"\t Percent negatives undetected by our rules + known obstructions: {unknown}/{false} = {float(unknown/false)}", color(unknown/false))
        norule = counts['norule'][(exps, slope)]
        # cprint(f"\t Percent negatives undetected by our rules + known obstructions: {norule}/{false} = {float(norule/false)}", color(norule/false))
        array[exps, slope] = [(len(exps), p_rank_dict[slope]), f'{undetected}/{false} = {float(undetected/false):.4f}', f'{unknown}/{false} = {float(unknown/false):.4f}']

    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    unknown = sum(counts['unknown'].values())
    norule = sum(counts['norule'].values())
    print()
    cprint(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
    cprint(f"Percent negatives undetected: {unknown}/{false} = {float(unknown/false)}", color(unknown/false))
    array['last'] = [r'\textbf{Total}', f"{undetected}/{false} = {float(undetected/false):.4f}", f"{unknown}/{false} = {float(unknown/false):.4f}"]
    # cprint(f"Percent negatives undetected: {norule}/{false} = {float(norule/false)}", color(norule/false))
    if table:
        table = Table([array[l] for l in array], headers)
        print(table)

def print_statistics_by_rule(data=None, table=False):
    if not data:
        data = get_data()
    
    counts = defaultdict(set)

    headers = ['Rule', 'Number of hits', 'Number of unique hits']
    array = dict()

    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule].add(poly[2])
    
    for rule in sorted(counts.keys(), key=lambda x:str(x)):
        other_hits = set()
        for r in counts:
            if r != rule:
                other_hits = other_hits.union(counts[r])
        print(f"Rule: {rule}, Number of hits: {len(counts[rule])}, Number of unique hits: {len(counts[rule] - other_hits.intersection(counts[rule]))}")
        array[rule] = [r'\textbf{' + rule + '}', len(counts[rule]), len(counts[rule] - other_hits.intersection(counts[rule]))]

    if table:
        table = Table([array[l] for l in array], headers)
        print(table)

def print_statistics_by_rule_and_type(data=None, table=False):
    if not data:
        data = get_data()

    for (exps, slope), D in data[False].items():
        counts = defaultdict(set)
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule].add(poly[2])
    
        for rule in sorted(counts.keys(), key=lambda x:str(x)):
            other_hits = set()
            for r in counts:
                if r != rule:
                    other_hits = other_hits.union(counts[r])
            print(f"Type: {exps}, {p_rank_dict[slope]}, Rule: {rule}, Number of hits: {len(counts[rule])}, Number of unique hits: {len(counts[rule] - other_hits.intersection(counts[rule]))}")

def known_jacobi_sort(data=None):
    if not data:
        data = get_data()

    data["unknown"].clear()
    count = 0
    for (exps, slope), D in data['undetected'].items():
        print(exps, slope)
        for q, polys in D.items():
            for poly in polys:
                iso_class = IsogenyClass(label = poly[2])
                point_counts_restriction = False
                q = iso_class.q
                for n in range(len(iso_class.curve_counts)):
                    if iso_class.curve_counts[n] > 2*(q**(n+1) + 1):
                        point_counts_restriction = True
                if iso_class.has_jacobian == 0 and not point_counts_restriction:
                    data["unknown"][(exps, slope)][q].append(poly)
                    count += 1
                
    print(count)
    pickle_data(data)
    # return data

def mod_2_checks(data):
    for pattern in all_coeff_patterns(2):
        for (exps, slope), D in data[False].items():
            flag = 0
            for q in [2, 4, 8, 16]:
                for poly in data[False][(exps, slope)][q]:
                    if tuple([c % 2 for c in poly[0][1:4]]) == pattern:
                        flag += 1
                    # if tuple([c % 2 for c in poly[1]]) == pattern:
                    #     flag += 1
                        
            if flag:
                flag2 = True
                for q in [2,4,8,16]:
                    for poly in data[True][exps, slope][q]:
                        if tuple([c%2 for c in poly[0][1:4]]) == pattern:
                        # if tuple([c%2 for c in poly[1]]) == pattern:
                            flag2 = False
                            # flag2 = False
            if flag2 and flag:
                print(exps, slope, pattern, flag)

def rule_set_comp(data):
    for (exps, slope), D in data[False].items():
        for q in [2, 4, 8, 16]:
            for poly in data[False][exps, slope][q]:
                pol = tuple([c % 2 for c in poly[0][1:4]])
                rules1 = [
                    exps == (6,) and slope == 0 and pol in ((0, 1, 1), (1, 0, 1)),
                    exps == (6,) and slope == 3 and pol == (1, 1, 0),
                    exps == (2, 4) and slope == 3 and pol == (1, 1, 0),
                ]
                rules1 = True in rules1
                coeff = tuple([c % 2 for c in poly[1]])
                rules2 = [
                    exps == (2, 4) and coeff == (0, 1, 1),
                    exps == (2, 4) and slope != 0 and coeff == (1, 0, 1),
                    exps == (6,) and slope != 0 and coeff in ((0, 0, 1), (1, 0, 1)),
                    exps == (6,) and slope == 0 and coeff in ((1, 0, 1), (0, 1, 1)),
                    exps == (6,) and slope == 3 and coeff[0] == 1,
                ]
                rules2 = True in rules2
                if rules1 != rules2:
                    print(exps, slope, rules1, rules2, poly)


def write_rule_statistics_table(stats, name):
    headers = ['Hypothesis','# Hits', '# Unique Rules', '# Rules Hit All Fields', '# Rules w/ Unique Hits', 'Top Rules']
    max_label_width = max(len(headers[0]), max(len(str(label)) for label in stats))
    col_widths = []
    for i in range(1, 6):
        max_cell = max(len(headers[i]), max((len(str(v[i-1])) for v in stats.values() if len(v) >= i), default=0))
        col_widths.append(max_cell)
    filename = f"{name}.txt"
    with open(filename, 'w') as f:
        f.write(' '.join([headers[0].ljust(max_label_width)] +
                         [headers[i].ljust(col_widths[i-1]) for i in range(1, 6)]) + '\n')
        for label, values in sorted(stats.items(), key=lambda item: -item[1][0]):
            cells = [str(values[i]) if i < len(values) else '' for i in range(5)]
            f.write(' '.join([label.ljust(max_label_width)] +
                             [cells[i].ljust(col_widths[i]) for i in range(5)]) + '\n')


def rule_overlap(data):
    slope_types, fac_patterns = [0,1,2,3,4], [(2,2,2), (2,4), (6,)]
    rule_names = set(range(2, 18)).union({'Costa', 'Thm 1.2', 'Thm 1.3', 'Thm 1.5', 'Proposition 1'})
    classified = {rule: set() for rule in rule_names}
    all_hits = set()
    for q in primes:
        for slope in slope_types:
            for exp_pat in fac_patterns:
                for poly in data[False][(exp_pat, slope)][q]:
                    rules = our_jacobi_rules(poly[0], slope, q, exp_pat, poly[1])
                    for k, v in rules.items():
                        if v:
                            classified[k].add(tuple(poly[0]))
                            all_hits.add(tuple(poly[0]))
    hits = {rule: set() for rule in rule_names}
    for k, v in classified.items():
        others = set().union(*(classified[r] for r in rule_names if r != k))
        hits[k] = (v - others, len(v - others), len(v))

    ## pairwise intersections
    filename = f"Intersections.txt"
    pairs = sorted([(A, B) for A in rule_names for B in rule_names if A!=B], key=lambda pair: -len(classified[pair[0]] & classified[pair[1]]))
    with open(filename, 'w') as f:
        for A, B in pairs:
            f.write(f"Rules {A} and {B}: size of {A}: {len(classified[A])}, size of {B}: {len(classified[B])}, size of intersection: {len(classified[A] & classified[B])}. \n")





    return hits

def rule_intersections(data, rules, hyp, primes_, exp_patterns, slopes, name):
    k = len(rules)
    num_to_rules, rule_to_num = {n: rules[n] for n in range(k)}, {rules[n]: n for n in range(k)}
    rule_hits = {n: set() for n in range(k)}
    for q in primes_:
        for slope in slopes:
            for exp_pat in exp_patterns:
                for poly in data[False][(exp_pat, slope)][q]:
                    if hyp(*poly[1]):
                        for rule in rules:
                            if power_rule_check(rule, poly[1], q):
                                rule_hits[rule_to_num[rule]].add(poly[1])
    pairs = sorted([(A, B) for A in range(k) for B in range(k) if A < B], key=lambda pair: -len(rule_hits[pair[0]] & rule_hits[pair[1]]))
    filename = f"Intersection of {name}.txt"
    with open(filename, 'w') as f:
        for A, B in pairs:
            f.write(f"Rules {num_to_rules[A]} and {num_to_rules[B]} of sizes {len(rule_hits[A])} and {len(rule_hits[B])}, size of intersection: {len(rule_hits[A] & rule_hits[B])}. \n")
