from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict

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

def get_data(field = None, prefetch=None):
    if not prefetch:
        print("prefetch")
        prefetch = prefetch_polys()
    R = PolynomialRing(ZZ, 'x')
    data = defaultdict(lambda: defaultdict(list))
    
    # set up keys for database querry
    fetch_keys = {"g":3}
    if field:
        fetch_keys["q"] = field

    # fetch data
    for rec in db.av_fq_isog.search(fetch_keys, ["label", "q", "hyp_count", "poly", "slopes", "simple_factors"]):
        q = rec["q"]
        slopes = tuple(sorted(rec["slopes"]))
        have_hyp = (rec["hyp_count"] > 0)
        # poly = R(rec["poly"])
        factors = rec["simple_factors"]
        for i in range(len(factors)):
            # factors[i] = db.av_fq_isog.lookup(factors[i][:-1], "poly")
            factors[i] = prefetch[factors[i][:-1]]
        exps = tuple([len(factors[i]) - 1 for i in range(len(factors))])
        coeffs = get_abc(factors, exps)
        data[q][exps, slope_types[slopes], have_hyp].append((rec["poly"], coeffs, rec["label"]))
    return data

def write_data(data):
    for q, D in data.items():
        for (exps, slopes, have_hyp), polys in D.items():
            
            filename = f"temp/data_{have_hyp}_{q}_{exps}_{slopes}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    # _ = F.write(", ".join(str(c) for c in poly) + "\n")
                    remainders = [r%modulo for r in poly[1]]
                    F.write(f"{poly[1]}, {remainders}, {poly[2]} \n")

def _in(arr, search, key=lambda x:x):
    for entry in arr:
        if key(entry) == search:
            return True
    return False

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

def mod_rule_check(data, data_sorted):
    for mod in range(2, 20):
        print(mod)
        for pattern in all_coeff_patterns(mod):
            for (exps, slopes, have_hyp) in pattern_keys():
                flag_no = {}
                for q in primes:
                    flag_yes = True
                    # if pattern doesn't appear in data True, flag_yes = True
                    for poly in data[q][(exps, slopes, have_hyp)]:
                        # coeffs = tuple((c % mod for c in poly[1]))
                        coeffs = tuple((poly[0][i+1] % mod for i in range(3)))
                        if coeffs == pattern:
                            flag_yes = False
                            break
                    # if pattern does appear in data False, record an entry in flag_no
                    if flag_yes:
                        if not (exps, slopes, False) in data[q]:
                            continue
                        for poly in data_sorted[q][(exps, slopes)]:
                            # coeffs = tuple((c % mod for c in poly[0]))
                            coeffs = tuple((poly[1][i+1] % mod for i in range(3)))
                            if coeffs == pattern:
                                flag_no[q] = len(data_sorted[q][(exps, slopes)])
                                break
                        # for poly in data[q][(exps, slopes, False)]:
                        #     coeffs = tuple((c % mod for c in poly[1]))
                        #     if coeffs == pattern:
                        #         flag_no[q] = True
                        #         break
                if len(flag_no) >= 1 and max(flag_no.values()) >= 1:
                    print(mod, pattern, (exps, slopes, have_hyp), flag_no)

def discriminant_check(data):
    """
    for curves of the form E^3, check for restrictions on the discriminant t^2-4q.
    (x^2 - tx + q)^3.
    """
    exp_pattern = (2,2,2)
    disc_results = {}
    obstructed_discriminants = []        
    for disc in range(-100, 100):
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
            print("discriminant {disc} with", len(disc_results[(disc, False)]))

    print(obstructed_discriminants)


def patterns_mod():
    for p in [2,3]:
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
                        divisor_results[(p, mod, exp, have_hyp, coeff)] = []
                        print((p, mod, exp, have_hyp, coeff))
                        for polyn in data[q][(exp, slopes, have_hyp)]:
                            a, b, c = polyn[1][0], polyn[1][1], polyn[1][2]
                            if coeff == 'a':
                                if a != 0:
                                    prime_div = prime_divisors(a)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                            if coeff == 'b':
                                if b != 0:
                                    prime_div = prime_divisors(b)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                            if coeff == 'c':
                                if c != 0:
                                    prime_div = prime_divisors(c)
                                    check = all([div%p == mod for div in prime_div])
                                    if check:
                                        divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                                elif mod == 0:
                                    divisor_results[(p, mod, exp, have_hyp, coeff)].append(polyn)
                if divisor_results[(p, mod, exp, True, coeff)] == [] and divisor_results[(p, mod, exp, False, coeff)] != []:
                    obstructed_div_mod_factoring.append((p, mod, exp, coeff))

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


def our_jacobi_rules(iso_class, poly, slopes, q, exps, coeffs):
    p = trial_division(q)
    R = PolynomialRing(ZZ, 'x')

    coeffs2 = tuple([c % 2 for c in coeffs])
    coeffs4 = tuple([c % 4 for c in coeffs])
    a,b,c = coeffs
    sorted_coeffs = sorted(coeffs)
    
    # if slope type 4 and simple
    rule1 = len(exps) == 1 and (slopes == 4) and (coeffs[0] == 0) and (coeffs[1] == 0) and ((coeffs[2] == -p*q) or (coeffs[2] == p*q))
    
    # if slope type 4 and characteristic is two
    rule2 = (slopes == 4) and (q % 2 == 0)

    # if slope type 3 and characteristic 2 and simple
    rule3 = (slopes == 3) and len(exps) == 1 and q%2 == 0 and coeffs[0]%2 == 1

    # if slope not type 0 and simple and coeffs are even even odd
    rule4 = coeffs2 == (0, 0, 1) and slopes != 0 and len(exps) == 1

    # if slope type 0 and simple and slopes are even even odd and characteristic is odd
    rule4_1 = coeffs2 == (0, 0, 1) and slopes == 0 and len(exps) == 1 and p%2 != 0

    # simple and coeffs follow odd even odd pattern
    rule5 = coeffs2 == (1, 0, 1) and len(exps) == 1
    
    # 2-4 splitting pattern and coeffs follow odd-even-odd and not ordinary
    rule6 = coeffs2 == (1, 0, 1) and len(exps) == 2 and slopes != 0

    # odd even odd pattern for 2-4 splitting, ordinary, and odd characteristic
    rule7 = coeffs2 == (1, 0, 1) and len(exps) == 2 and slopes == 0 and p != 2

    # split into quadratics and coefficients follow odd odd odd pattern and characteristic not 2
    rule8 = len(exps) == 3 and coeffs2 == (1, 1, 1) and p != 2

    # 2x4 splitting, slope type 2, prime field != 2, (1, 0, 2) mod 4 or (3, 0, 2) mod 4 pattern
    rule9 = len(exps) == 2 and (coeffs4 == (1, 0, 2) or coeffs4 == (3, 0, 2)) and p == q and p != 2 and slopes == 2

    # simple and slope type 0 and coefficients are even, odd, odd and characteristic 2
    rule10 = len(exps) == 1 and slopes == 0 and coeffs2 == (0, 1, 1) and p == 2

    # 2x4 splitting and slope type 3 and coefficients are even, odd, odd and characteristic 2
    rule11 = len(exps) == 2 and slopes == 3 and coeffs2 == (0, 1, 1) and p == 2

    # 2x2x2 splitting of the form (x^2 - tx + q)^3, conditions on discriminant t^2 - 4q
    rule12 = False
    obstructed_discriminants = [-99, -91, -83, -75, -67, -59, -51, -43, -39, -35, -27, -23, -20, -19, -15, -12, -11, -8, -7, -4, -3, 0]
    if len(exps) == 3:
        if coeffs[0]==coeffs[1]==coeffs[2] and (coeffs[0]**2 - 4*q in obstructed_discriminants):
            rule12 = True
    
    # 2x2x2 splitting with parameters a >= b >= c, checks that a - b = b - c = 1
    rule13 = False
    if len(exps) == 3:
        rule13 = sorted_coeffs[2] - sorted_coeffs[1] == 1 and sorted_coeffs[1] - sorted_coeffs[0] == 1

    
    # 2x2x2 splitting, if two are equal and third is not (so abs(a - b) + abs(b - c) + abs(a - c) == 2)
    rule14 = len(exps) == 3 and abs(coeffs[0] - coeffs[1]) + abs(coeffs[1] - coeffs[2]) + abs(coeffs[0] - coeffs[2]) == 2

    # 2x2x2 splitting, p-rank=0, for q=2,3,4,9 if a^2 = b^2 = c^2 and p | a^2
    rule15 = False
    if len(exps) == 3 and slopes in {1,4} and q in {2,3,4,9}:
        rule15 = coeffs[0]**2 == coeffs[1]**2 and coeffs[1]**2 == coeffs[2]**2 and coeffs[1]%p == 0

    '''
    # 2x4 with slope type 4, if first two coefficients are 0 and last coefficient is -2q with q > 3 and squarefree
    rule16 = len(exps) == 2 and q > 3 and q == p and slopes == 4 and sorted_coeffs == [-2*q, 0, 0]
    '''

    # if iso_class.has_jacobian == -1 or rule1 or rule2 or rule3:
    if rule1 or rule2 or rule3 or rule4 or rule5 or rule6 or rule7 or rule8 or rule4_1 or rule9 or rule10 or rule11 or rule12 or rule13 or rule14 or rule15 or rule16:
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

def sort_data(data):
    R = PolynomialRing(ZZ, 'x')
    new_obstructions = []
    known_obstructions = []
    jacobians = []
    false_positives = []
    false_positives_data = defaultdict(lambda: defaultdict(list))
    new_obstr_data = defaultdict(lambda: defaultdict(list))
    for q, D in data.items():
        for (exps, slope, have_hyp), polys in D.items():
            for polynom_ in polys:
                polynom = polynom_[0]
                # iso_class = IsogenyClass(poly=polynom)
                iso_class = None
                coeffs = polynom_[1]
                coeffs_mod_2 = tuple([c % modulo for c in coeffs])
                guess = our_jacobi_rules(iso_class, polynom, slope, q, exps, coeffs)
                if have_hyp and guess == True:
                    jacobians.append(polynom)
                elif have_hyp and guess == False:
                    false_positives.append(polynom)
                    
                    false_positives_data[q][exps, slope].append((coeffs, polynom))
                elif not have_hyp and guess==True:
                    new_obstructions.append(polynom)
                    # iso_class = IsogenyClass(poly=polynom)
                    # if iso_class.has_jacobian == 0:
                    #     new_obstr_data[q][exps, slope].append((coeffs, polynom))
                    new_obstr_data[q][exps, slope].append((coeffs, polynom))
                
                
                else:
                    known_obstructions.append(polynom)
    
    # return new_obstructions, known_obstructions, jacobians
    print(f"Ratio false positive: {len(false_positives)/len(jacobians)}")
    print(f"Percent detected negatives: {(len(known_obstructions)/len(new_obstructions))/(len(known_obstructions)/len(new_obstructions) + 1)}")
    # print(false_positives_data)
    print([[(key, len(false_positives_data[f][key]), f) for key in false_positives_data[f]] for f in false_positives_data])
    return new_obstr_data, false_positives_data


def write_new_obstructions_data(data):
    for q, D in data[0].items():
        for (exps, slopes), polys in D.items():
            
            filename = f"temp/data_norule_{q}_{exps}_{slopes}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    _ = F.write(", ".join(str(c) for c in poly) + "\n")

def write_false_positives(data):
    for q, D in data[1].items():
        for (exps, slopes), polys in D.items():
            
            filename = f"temp/false_pos_{q}_{exps}_{slopes}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    _ = F.write(", ".join(str(c) for c in poly) + "\n")
    
                    
def write_exp_data(data):
    write_new_obstructions_data(data)
    write_false_positives(data)


                

