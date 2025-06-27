from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict
import itertools
import pickle
import os
from functools import partial
from termcolor import cprint
from pprint import pprint
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

primes = [1, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101, 103, 107, 109, 113, 121, 125, 127, 128, 131, 137, 139, 149, 151, 157, 163, 167, 169, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227]


def prefetch_polys():
    data = defaultdict(list)
    for field in primes:
        for rec in db.av_fq_isog.search({"q":field, "g":2}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        for rec in db.av_fq_isog.search({"q":field, "g":1}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        print(field)
    return data

def pickle_data(data):
    with open("data_g2.dict", 'wb') as filename:
        pickle.dump(data, filename)

def get_data_initial(prefetch=None):
    if os.path.isfile("data_g2.dict"):
        raise Exception("The data already exists. Are you sure you want to re-import")
    if not prefetch:
        print("prefetch")
        prefetch = prefetch_polys()
    # querry by (factoring, slope type); field q; hyperelliptic classification: returns the polynomial, coefficients in factored form, and label
    data = defaultdict(partial(defaultdict, partial(defaultdict, list)))
    
    # set up keys for database querry
    for q in primes:
        fetch_keys = {"g":2, "q":q}

        # fetch data
        for rec in db.av_fq_isog.search(fetch_keys, ["label", "hyp_count", "poly", "p_rank", "simple_factors"]):
            p_rank = rec["p_rank"]
            have_hyp = (rec["hyp_count"] > 0)
            factors = rec["simple_factors"]
            for i in range(len(factors)):
                factors[i] = prefetch[factors[i][:-1]]
            factoring_type = tuple([len(factors[i]) - 1 for i in range(len(factors))])
            coeffs = get_ab(factors, factoring_type)
            data[have_hyp][factoring_type, p_rank][q].append((rec["poly"], coeffs, rec["label"]))

    pickle_data(data)
    
def get_data():
    with open("data_g2.dict", 'rb') as filename:
        return(pickle.load(filename))

def write_data(data):
    for have_hyp, D in data.items():
        for (factoring, slopes), L2 in D.items():
            for q, polys in L2.items():
                filename = f"genus2/data_{have_hyp}_{q}_{factoring}_{slopes}.txt"
                with open(filename, "w") as F:
                    for poly in polys:
                        F.write(f"{poly[1]}, {poly[0]}, {poly[2]}, {poly[1][0]**2 - poly[1][1]**2  + poly[1][2]} \n")

def get_ab(factors, classification):
    if classification == (4, ):
        a = factors[0][1]
        b = factors[0][2]
        
    elif classification == (2, 2):
        a = factors[0][1]
        b = factors[1][1]

    return (a, b)

def our_jacobi_rules(poly, p_rank, q, exps, coeffs):
    p = trial_division(q)

    s,t = coeffs
    if abs(s) < abs(t):
        s, t = t, s
    s *= -1
    t *= -1
    a, b = coeffs
    
    rules = {}
    
    if len(exps) == 2:
        rules[0] = abs(s - t) == 1
        rules[1] = p_rank == 2 and s == t and t**2 - 4*q in {-3, -4, -7}
        rules[2] = p_rank == 2 and q == 2 and abs(s) == 1 and abs(t) == 1 and s != t
        rules[3] = p_rank == 1 and is_square(q) and s**2 == 4*q and is_squarefree(s-t)
        if p_rank == 0:
            rules[5] = p > 3 and s**2 != t**2
            rules[6] = p == 3 and not is_square(q) and s**2 == t**2 and s**2 == 3*q
            rules[7] = p == 3 and is_square(q) and (s - t) % (3*(q**0.5)) != 0
            rules[8] = p == 2 and (s**2 - t**2) % (2*q) != 0
            rules[9] = q in (2,3) and s == t
            rules[10] = q in (4, 9) and s**2 == t**2 and s**2 == 4*q

    else:
        rules[11] = a**2 - b == q and b < 0 and not (False in [div % 3 == 1 for div in prime_divisors(b)])
        rules[12] = p_rank == 2 and a == 0 and b == 1 - 2*q
        rules[13] = p_rank == 2 and p > 2 and a == 0 and b == 2 - 2*q
        if p_rank == 0:
            rules[14] = p % 12 == 11 and is_square(q) and a == 0 and b == -q
            rules[15] = p == 3 and is_square(q) and a == 0 and b == -q
            rules[16] = p == 2 and not is_square(q) and a == 0 and b == -q
            rules[17] = q in (2, 3) and a == 0 and b == -2*q

    
    return rules

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

def print_statistics_by_rule(data=None):
    if not data:
        data = get_data()
    
    counts = defaultdict(int)

    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule] += 1
    
    pprint(counts)