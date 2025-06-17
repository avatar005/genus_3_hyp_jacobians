from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict
import itertools
import pickle
import os
from functools import partial
from termcolor import cprint

slope_types = {('0A', '0B', '0C', '1A', '1B', '1C'): 0,
               ('1/2A', '1/2B', '1/2C', '1/2D', '1/2E', '1/2F'): 4,
               ('0A', '1/2A', '1/2B', '1/2C', '1/2D', '1A'): 2,
               ('0A', '0B', '1/2A', '1/2B', '1A', '1B'): 3,
               ('1/3A', '1/3B', '1/3C', '2/3A', '2/3B', '2/3C'): 1}

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
    # querry by (factoring, slope type); field q; hyperelliptic classification: returns the polynomial, coefficients in factored form, and label
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
                        F.write(f"{poly[1]}, {poly[0]}, {poly[2]} \n")

def all_coeff_patterns(mod):
    for i in range(mod):
        for j in range(mod):
            for k in range(mod):
                yield (i, j, k)

exp_patterns = ((6,), (2, 4), (2, 2, 2))
def pattern_keys():
    for i in exp_patterns:
        for j in range(5):
            yield (i, j, True)

def n_tuple(length, lower, upper):
    tup = [range(lower, upper + 1) for i in range(length - 1)]
    tup.append(range(lower*4, upper*4))
    return itertools.product(*tup)

def evaluate(coeffs, input):
    a, b, c, d, e, f, g = coeffs
    x1, x2, x3 = input
    return x1**2 * a + x2**2 * b + x3**2 * c + a*b*d + b*c*e + c*a*f + g

def poly_rule_check(data, data_sorted):
    # for mod in range(2, 20):
    #     print(mod)
    for pattern in n_tuple(4, -1, 4):
        print(pattern)
        for (exps, slopes, have_hyp) in pattern_keys():
            flag_no = {}
            for q in primes:
                flag_yes = True
                # if pattern doesn't appear in data True, flag_yes = True
                for poly in data[q][(exps, slopes, have_hyp)]:
                    # coeffs = tuple((c % mod for c in poly[1]))
                    coeffs = poly[1]
                    if 0 == evaluate(pattern, coeffs):
                        flag_yes = False
                        break
                # if pattern does appear in data False, record an entry in flag_no
                if flag_yes:
                    if not (exps, slopes, False) in data[q]:
                        continue
                    for poly in data_sorted[q][(exps, slopes, False)]:
                        # coeffs = tuple((c % mod for c in poly[0]))
                        coeffs = poly[0]
                        if 0 == evaluate(pattern, coeffs):
                            print(pattern, coeffs)
                            if q in flag_no:
                                flag_no[q] += 1
                            else:
                                flag_no[q] = 1
                            # flag_no[q] = len(data_sorted[q][(exps, slopes)])
                            
            if len(flag_no) >= 1 and max(flag_no.values()) >= 1:
                print(pattern, (exps, slopes, have_hyp), flag_no)

def poly_rule_check_optimized(data, data_sorted):
    # for mod in range(2, 20):
    #     print(mod)
    p = 0
    for pattern in n_tuple(7, -1, 4):
        if pattern[1] != p:
            p = pattern[1]
            print(pattern)
        for (exps, slopes, have_hyp) in pattern_keys():
            flag_no = {}
            for q in primes:
                for poly in data_sorted[(exps, slopes)][q][False]:
                    # coeffs = tuple((c % mod for c in poly[0]))
                    coeffs = poly[0]
                    if 0 == evaluate(pattern, coeffs):
                        # print(pattern, coeffs)
                        if q in flag_no:
                            flag_no[q] += 1
                        else:
                            flag_no[q] = 1
                if q in flag_no:
                    # if pattern doesn't appear in data True, flag_yes = True
                    for poly in data[q][(exps, slopes, have_hyp)]:
                        # coeffs = tuple((c % mod for c in poly[1]))
                        coeffs = poly[1]
                        if 0 == evaluate(pattern, coeffs):
                            del flag_no[q]
                            break
                            
            if len(flag_no) >= 1 and max(flag_no.values()) >= 3:
                print(pattern, (exps, slopes, have_hyp), flag_no)

def discriminant_check(data):
    """
    for curves of the form E^3, check for restrictions on the discriminant t^2-4q.
    (x^2 - tx + q)^3.
    """
    exp_pattern = (2,2,2)
    disc_results = {}
    obstructed_discriminants = []        
    for disc in range(-100, 0):
        disc_results[(disc, True)] = []
        disc_results[(disc, False)] = []
        for slopes in range(5):
            for have_hyp in {True, False}:
                for q in primes:
                    for polyn in data[q][(exp_pattern, slopes, have_hyp)]:
                        coeffs = polyn[1]
                        #check the form E^3
                        if coeffs[0] == coeffs[1] == coeffs[2]:
                            discriminant = coeffs[0]**2 - 4*q
                            if discriminant == disc:
                                disc_results[(disc, have_hyp)].append(polyn)
        
        if disc_results[(disc, True)] == [] and disc_results[(disc, False)] != []:
            obstructed_discriminants.append(disc)
            print(f"discriminant {disc} with # {len(disc_results[(disc, False)])}")

    print(obstructed_discriminants)


def patterns_mod():
    for p in [2,3,5]:
        for mod in range(p):
            for coeff in ['a', 'b', 'c']:
                yield(p, mod, coeff)


def prime_divisors_check(data):
    divisor_results = {}
    obstructed_div_mod_factoring = []
    for q in primes:
        for exp in exp_patterns: 
            for p, mod, coeff in patterns_mod():
                for have_hyp in {True, False}:
                    for slopes in slope_types:
                        divisor_results[(p, mod, exp, have_hyp, slopes)] = []
                        print((p, mod, exp, have_hyp, coeff))
                        for polyn in data[q][(exp, slopes, have_hyp)]:
                            a, b, c = polyn[1][0], polyn[1][1], polyn[1][2]
                            if coeff == 'a':
                                if a != 0:
                                    prime_div = prime_divisors(a)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                            if coeff == 'b':
                                if b != 0:
                                    prime_div = prime_divisors(b)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                            if coeff == 'c':
                                if c != 0:
                                    prime_div = prime_divisors(c)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, slopes)].append(polyn)
                if divisor_results[(p, mod, exp, True, slopes)] == [] and divisor_results[(p, mod, exp, False, slopes)] != []:
                    obstructed_div_mod_factoring.append((p, mod, exp, slopes))

    print('done')
    print(obstructed_div_mod_factoring)

def distances():
    for i in range(100):
        for j in range(100):
            yield (i, j)

def distance_check(data):
    exp_pattern = (2,2,2)
    dist_results = {}
    obstructed_dist = []
    for q in primes:
        for slopes in range(5):
            for dist in distances():
                dist_results[(dist, True)] = []
                dist_results[(dist, False)] = []
                for have_hyp in {True, False}:
                        for polyn in data[q][(exp_pattern, slopes, have_hyp)]:
                            coeffs = sorted(polyn[1])
                            a, b, c = coeffs[2], coeffs[1], coeffs[0]
                            if a - b == dist[0] and b - c == dist[1]:
                                dist_results[(dist, have_hyp)].append(polyn)
                
                if dist_results[(dist, True)] == [] and dist_results[(dist, False)] != []:
                    obstructed_dist.append(dist)
                    print(q, slopes, dist, len(dist_results[(dist, False)]))

    print(obstructed_dist)
    # print(dist_results[0, False], dist_results[0,True])

def parity_check_2(data):
    exp_pattern = (2,4)
    parity_results = {}
    obstructed_polys = []
    slopes = 2
    for q in {2,4,8,16}:
        for hyp in {True, False}:
            for div in range(2,4):
                for mod1 in range(div):
                    for mod2 in range(div):
                        for mod3 in range(div):
                            parity_results[(div, mod1, mod2, mod3, hyp)] = []
                            for polyn in data[q][(exp_pattern, slopes, hyp)]:
                                coeffs = polyn[1]
                                if coeffs[0]%div == mod1 and coeffs[1]%div == mod2 and coeffs[2]%div == mod3:
                                    parity_results[(div, mod1, mod2, mod3, hyp)].append(polyn)
    for div in range(2,4):
        for mod1 in range(div):
            for mod2 in range(div):
                for mod3 in range(div):
                    if parity_results[(div, mod1, mod2, mod3, True)]==[] and parity_results[(div, mod1, mod2, mod3, False)]!=[]:
                        print((div, mod1, mod2, mod3), len(parity_results[(div, mod1, mod2, mod3, False)]))
    # print(parity_results[True], "break", parity_results[False])
    print('done')

def equality_check(data):
    #this is for p-rank zero, irreducible
    exp_pattern = (2,2,2)
    exp_slopes = {1,4}
    obstructed_q_results = {}
    obs_q = []
    for q in primes:
        # p = char[q]
        for have_hyp in {True, False}:
            obstructed_q_results[(q, have_hyp)] = []
            for slopes in exp_slopes:
                for polyn in data[q][(exp_pattern, slopes, have_hyp)]:
                    a, b, c = polyn[1][0], polyn[1][1], polyn[1][2]
                    if a**2 == b**2 and b**2 == c**2:
                        obstructed_q_results[(q, have_hyp)].append([a,b,c])
        if obstructed_q_results[(q, True)] == [] and obstructed_q_results[(q, False)] != []:
            obs_q.append(q)
            print(q, obstructed_q_results[(q, False)])

    print(obs_q)


def our_jacobi_rules(poly, slopes, q, exps, coeffs):
    p = trial_division(q)

    coeffs4 = tuple([c % 4 for c in coeffs])
    poly2 = tuple([c % 2 for c in poly[1:4]])
    sorted_coeffs = sorted(coeffs)
    a,b,c = sorted_coeffs
    
    rules = {}

    # Costa Theorem 2.8
    rules['Costa'] = poly2[1:] == (0, 1) and p%2 != 0

    # Case 1 in Theorem 1
    rules['Thm 1.1'] = poly2 == (0, 0, 1) and slopes != 0 and p == 2

    # Case 2 of Theorem 1
    rules['Thm 1.2'] = poly2 == (0, 1, 1) and p == 2

    # Case 3 in Theorem 1
    rules['Thm 1.3'] = poly2 == (1, 0, 1) and slopes != 3 and p == 2

    # Case 4 in Theorem 1
    rules['Thm 1.4'] = poly2 == (1, 1, 1) and (not slopes in (0, 2)) and p == 2

    # Case 5 in Theorem 1
    rules['Thm 1.5'] = (slopes == 3) and p == 2 and poly2 == (1, 1, 0)


    # STRANGE AND WORTH EXPLORING FURTHER 
    # if slope type 4 and characteristic is two
    rules[1] = (slopes == 4) and (q % 2 == 0)

    # 2x4 splitting, slope type 2, prime field != 2, (1, 0, 2) mod 4 or (3, 0, 2) mod 4 pattern
    rules[2] = len(exps) == 2 and (coeffs4 == (1, 0, 2) or coeffs4 == (3, 0, 2)) and p == q and p != 2 and slopes == 2
    # STRANGE AND WORTH EXPLORING FURTHER ^


    # 2x2x2 splitting of the form (x^2 - tx + q)^3, conditions on discriminant t^2 - 4q
    rules[3] = False
    obstructed_discriminants = [-99, -91, -83, -75, -67, -59, -51, -43, -39, -35, -27, -23, -20, -19, -15, -12, -11, -8, -7, -4, -3, 0]
    if len(exps) == 3:
        if a==b==c and (a**2 - 4*q in obstructed_discriminants):
            rules[3] = True
    
    # 2x2x2 splitting with parameters a >= b >= c, checks that a - b = b - c = 1
    rules[4] = False
    if len(exps) == 3:
        rules[4] = c - b == 1 and b - a == 1

    
    # 2x2x2 splitting, if two are equal and third is not (so abs(a - b) + abs(b - c) + abs(a - c) == 2)
    rules[5] = len(exps) == 3 and abs(a - b) + abs(b - c) + abs(a - c) == 2

    # 2x2x2 splitting, p-rank=0, for q=2,3,4,9 if a^2 = b^2 = c^2 and p | a^2
    rules[6] = False
    if len(exps) == 3 and slopes in {1,4} and q in {2,3,4,9}:
        rules[6] = a**2 == b**2 and b**2 == c**2 and b % p == 0

    # rules[7] = abs(poly[1]) > q + 1

    # if iso_class.has_jacobian == -1 or rule1 or rule2 or rule3:
    if True in rules.values():
        return False
    return True


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
                if guess:
                    counts['undetected'][(exps, slope)] += 1
                    data['undetected'][(exps, slope)][q].append(poly)
    
    for (exps, slope), D in data[True].items():
        for q, polys in D.items():
            for poly in polys:
                counts[True][(exps, slope)] += 1
                guess = our_jacobi_rules(poly[0], slope, q, exps, poly[1])
                if not guess:
                    counts['false positives'][(exps, slope)] += 1
                    data['false positives'][(exps, slope)][q].append(poly)
    
    fp = sum(counts['false positives'].values())
    true = sum(counts[True].values())
    print(f"Percent false positive: {fp}/{true} = {float(fp/true)}")
    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    print(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}")

    pickle_data(data)

    return data

def print_statistics(data=None):
    color = lambda ratio: (255*(ratio)**0.5, 255 - 255*(ratio)**0.5, 0)
    if not data:
        data = get_data()

    counts = defaultdict(lambda: defaultdict(int))

    for have_hyp, D1 in data.items():
        for (exps, slope), D2 in D1.items():
            for q, polys in D2.items():
                counts[have_hyp][(exps, slope)] += len(polys)

    for (exps, slope) in sorted(counts[False]):
        print(exps, slope)
        undetected = counts['undetected'][(exps, slope)]
        false = counts[False][(exps, slope)]
        cprint(f"\t Percent negatives undetected by our rules: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
        unknown = counts['unknown'][(exps, slope)]
        cprint(f"\t Percent negatives undetected by our rules + known obstructions: {unknown}/{false} = {float(unknown/false)}", color(unknown/false))

    undetected = sum(counts['undetected'].values())
    false = sum(counts[False].values())
    unknown = sum(counts['unknown'].values())
    print()
    cprint(f"Percent negatives undetected: {undetected}/{false} = {float(undetected/false)}", color(undetected/false))
    cprint(f"Percent negatives undetected: {unknown}/{false} = {float(unknown/false)}", color(unknown/false))


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
                    if iso_class.curve_counts[n] > min(2*(q**(n+1) + 1), q**(n+1) + 1 + 2*3*q**((n+1)/2)):
                        point_counts_restriction = True
                if iso_class.has_jacobian == 0 and not point_counts_restriction:
                    data["unknown"][(exps, slope)][q].append(poly)
                    count += 1
                
    print(count)
    pickle_data(data)
    return data



# def write_new_obstructions_data(data):
#     for q, D in data[0].items():
#         for (exps, slopes), polys in D.items():
            
#             filename = f"temp/data_norule_{q}_{exps}_{slopes}.txt"
#             with open(filename, "w") as F:
#                 for poly in polys:
#                     _ = F.write(", ".join(str(c) for c in poly) + "\n")

# def write_false_positives(data):
#     for q, D in data[1].items():
#         for (exps, slopes), polys in D.items():
            
#             filename = f"temp/false_pos_{q}_{exps}_{slopes}.txt"
#             with open(filename, "w") as F:
#                 for poly in polys:
#                     _ = F.write(", ".join(str(c) for c in poly) + "\n")
    
                    
# def write_exp_data(data):
#     write_new_obstructions_data(data)
#     write_false_positives(data)


                

