from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict

slope_types = {('0A', '0B', '0C', '1A', '1B', '1C'): 0,
               ('1/2A', '1/2B', '1/2C', '1/2D', '1/2E', '1/2F'): 4,
               ('0A', '1/2A', '1/2B', '1/2C', '1/2D', '1A'): 2,
               ('0A', '0B', '1/2A', '1/2B', '1A', '1B'): 3,
               ('1/3A', '1/3B', '1/3C', '2/3A', '2/3B', '2/3C'): 1}

modulo = 4

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


def our_jacobi_rules(iso_class, poly, slopes, q, exps, coeffs):
    p = trial_division(q)
    R = PolynomialRing(ZZ, 'x')
    
    # if slope type 4 and simple
    rule1 = len(exps) == 1 and (slopes == 4) and (coeffs[0] == 0) and (coeffs[1] == 0) and ((coeffs[2] == -p*q) or (coeffs[2] == p*q))
    
    # if slope type 4 and characteristic is two
    rule2 = (slopes == 4) and (q % 2 == 0)

    # if slope type 3 and characteristic 2 and simple
    rule3 = (slopes == 3) and len(exps) == 1 and q%2 == 0 and coeffs[0]%2 == 1

    # if slope not type 0 and simple and slopes are even even odd
    rule4 = coeffs[0]%2 == 0 and coeffs[1]%2 == 0 and coeffs[2]%2 == 1 and slopes != 0 and len(exps) == 1

    # simple and coeffs follow odd even odd pattern
    rule5 = coeffs[0]%2 == 1 and coeffs[1]%2 == 0 and coeffs[2]%2 == 1 and len(exps) == 1
    
    # 2-4 splitting pattern and coeffs follow odd-even-odd and not ordinary
    rule6 = coeffs[0]%2 == 1 and coeffs[1]%2 == 0 and coeffs[2]%2 == 1 and len(exps) == 2 and slopes != 0

    # odd even odd pattern for 2-4 splitting, ordinary, and odd characteristic
    rule7 = coeffs[0]%2 == 1 and coeffs[1]%2 == 0 and coeffs[2]%2 == 1 and len(exps) == 2 and slopes == 0 and p != 2

    # if iso_class.has_jacobian == -1 or rule1 or rule2 or rule3:
    if rule1 or rule2 or rule3 or rule4 or rule5 or rule6 or rule7:
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
                    
                    false_positives_data[q][exps, slope].append((coeffs, coeffs_mod_2))
                elif not have_hyp and guess==True:
                    new_obstructions.append(polynom)

                    new_obstr_data[q][exps, slope].append((coeffs, coeffs_mod_2))
                
                
                else:
                    known_obstructions.append(polynom)
    
    # return new_obstructions, known_obstructions, jacobians
    print(f"Ratio false positive: {len(false_positives)/len(jacobians)}")
    print(f"Percent detected negatives: {(len(known_obstructions)/len(new_obstructions))/(len(known_obstructions)/len(new_obstructions) + 1)}")
    # print(false_positives_data)
    print([[(key, len(false_positives_data[f][key])) for key in false_positives_data[f]] for f in false_positives_data])
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


                

