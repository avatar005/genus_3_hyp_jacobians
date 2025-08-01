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

primes = [1, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25,
          27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64,
          67, 71, 73, 79, 81, 83, 89, 97, 101, 103, 107, 109, 
          113, 121, 125, 127, 128, 131, 137, 139, 149, 151, 157, 
          163, 167, 169, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227]

def prefetch_polys():
    '''Prepare a dictionary of label:polynomial pairs to speed up downloading data from the database'''

    data = defaultdict(list)
    primes = [2,3,4,5,7,8,9,11,13,16,17,19,23,25]
    for field in primes:
        for rec in db.av_fq_isog.search({"q":field}, ["label", "poly", "hyp_count"]):
            data[rec["label"]] = (rec["poly"], rec['hyp_count'] > 0)
        print(field)
    return data


def test_jacobians(prefetch=None):
    '''Downloads data from database. Pickles resulting fetched data as a dictionary keyed by hyperelliptic classification; (factoring, slope type), field. 
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

    if not prefetch:
        print("prefetch")
        prefetch = prefetch_polys()
    # querry by hyperelliptic classification; (factoring, slope type); field q: returns the polynomial, coefficients in factored form, and label
    data = defaultdict(partial(defaultdict, partial(defaultdict, list)))
    
    # set up keys for database querry
    fetch_keys = {"g":3, "is_simple": False}

    # fetch data
    for rec in db.av_fq_isog.search(fetch_keys, ["label", "q", "hyp_count", "poly", "slopes", "simple_factors"]):
        # q = rec["q"]
        # have_hyp = (rec["hyp_count"] > 0)
        factors = rec["simple_factors"]
        for i in range(len(factors)):
            factors[i] = prefetch[factors[i][:-1]][1]
        if rec["hyp_count"] == 0 and not all(factors):
            print(rec)
        # factoring_type = tuple([len(factors[i]) - 1 for i in range(len(factors))])
        # coeffs = get_abc(factors, factoring_type)
        # data[have_hyp][factoring_type, slope][q].append((rec["poly"], coeffs, rec["label"]))



def prefetch_polys_g2():
    data = defaultdict(list)
    for field in primes:
        for rec in db.av_fq_isog.search({"q":field, "g":2}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        for rec in db.av_fq_isog.search({"q":field, "g":1}, ["label", "poly"]):
            data[rec["label"]] = rec["poly"]
        print(field)
    return data

def pickle_data_g2(data):
    with open("data_g2.dict", 'wb') as filename:
        pickle.dump(data, filename)

def get_data_initial_g2(prefetch=None):
    if os.path.isfile("data_g2.dict"):
        raise Exception("The data already exists. Are you sure you want to re-import")
    if not prefetch:
        print("prefetch")
        prefetch = prefetch_polys_g2()
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
            coeffs = get_ab_g2(factors, factoring_type)
            data[have_hyp][factoring_type, p_rank][q].append((rec["poly"], coeffs, rec["label"]))

    pickle_data_g2(data)
    
def get_data_g2():
    with open("data_g2.dict", 'rb') as filename:
        return(pickle.load(filename))

def write_data_g2(data):
    for have_hyp, D in data.items():
        for (factoring, slopes), L2 in D.items():
            for q, polys in L2.items():
                filename = f"genus2/data_{have_hyp}_{q}_{factoring}_{slopes}.txt"
                with open(filename, "w") as F:
                    for poly in polys:
                        F.write(f"{poly[1]}, {poly[0]}, {poly[2]}, {poly[1][0]**2 - poly[1][1]**2  + poly[1][2]} \n")

def get_ab_g2(factors, classification):
    if classification == (4, ):
        a = factors[0][1]
        b = factors[0][2]
        
    elif classification == (2, 2):
        a = factors[0][1]
        b = factors[1][1]

    return (a, b)

def p_rank_from_coeffs(coeffs, p):
    if (coeffs[0] % p != 0 or coeffs[0] == 0) and (coeffs[1] % p != 0 or coeffs[1] == 0):
        return 2
    if (coeffs[0] % p != 0 or coeffs[0] == 0) and (coeffs[1] % p == 0 and coeffs[1] != 0):
        return 1
    return 0

def test_g3_breakdowns():
    def get_data():
        with open("data.dict", 'rb') as filename:
            return(pickle.load(filename))
    
    data = get_data()

    for s in range(5):
        count = 0
        for q, polys in data[False][(2, 4), s].items():
            p = trial_division(q)
            for poly in polys:
                p_rank = p_rank_from_coeffs(poly[1][1:], p)
                if True in our_jacobi_rules_g2(p_rank, q, (4,), poly[1][1:]).values():
                    count += 1
                    # print(poly)
        print(count)
        count = 0
        for q, polys in data[True][(2, 4), s].items():
            for poly in polys:
                p_rank = p_rank_from_coeffs(poly[1][1:], p)
                if True in our_jacobi_rules_g2(p_rank, q, (4,), poly[1][1:]).values():
                    count += 1
                    print(poly, db.av_fq_isog.lookup(poly[2], "hyp_count"))
        print(count)


def our_jacobi_rules_g2(p_rank, q, exps, coeffs):
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

def sort_data_g2(data=None):
    if not data:
        data = get_data_g2()
    
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
                guess = our_jacobi_rules_g2(slope, q, exps, poly[1])
                if not (True in guess.values()):
                    counts['undetected'][(exps, slope)] += 1
                    data['undetected'][(exps, slope)][q].append(poly)
    
    for (exps, slope), D in data[True].items():
        for q, polys in D.items():
            for poly in polys:
                counts[True][(exps, slope)] += 1
                guess = our_jacobi_rules_g2(slope, q, exps, poly[1])
                if True in guess.values():
                    counts['false positives'][(exps, slope)] += 1
                    data['false positives'][(exps, slope)][q].append(poly)
    
    fp = sum(counts['false positives'].values())
    true = sum(counts[True].values())
    print(f"Percent false positive: {fp}/{true} = {float(fp/true)}")
    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    print(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}")

    pickle_data_g2(data)

def print_statistics_by_rule_g2(data=None):
    if not data:
        data = get_data_g2()
    
    counts = defaultdict(int)

    for (exps, slope), D in data[False].items():
        for q, polys in D.items():
            for poly in polys:
                r = our_jacobi_rules_g2(poly[0], slope, q, exps, poly[1])
                for rule in r:
                    if r[rule] == True:
                        counts[rule] += 1
    
    pprint(counts)