from lmf import db
from sage.all import PolynomialRing, ZZ
from collections import defaultdict

slope_types = {('0A', '0B', '0C', '1A', '1B', '1C'): 0,
               ('1/2A', '1/2B', '1/2C', '1/2D', '1/2E', '1/2F'): 4,
               ('0A', '1/2A', '1/2B', '1/2C', '1/2D', '1A'): 2,
               ('0A', '0B', '1/2A', '1/2B', '1A', '1B'): 3,
               ('1/3A', '1/3B', '1/3C', '2/3A', '2/3B', '2/3C'): 1}

def get_data(field = None):
    R = PolynomialRing(ZZ, 'x')
    data = defaultdict(lambda: defaultdict(list))
    fetch_keys = {"g":3}
    if field:
        fetch_keys["q"] = field
    for rec in db.av_fq_isog.search(fetch_keys, ["label", "q", "hyp_count", "poly", "slopes"]):
        q = rec["q"]
        slopes = tuple(sorted(rec["slopes"]))
        have_hyp = (rec["hyp_count"] > 0)
        poly = R(rec["poly"])
        factors = poly.factor()
        exps = classify_factors(factors)
        coeffs = get_abc(factors, exps)
        data[q][exps, slope_types[slopes], have_hyp].append((rec["poly"], coeffs))
    return data

def write_data(data):
    for q, D in data.items():
        for (exps, slopes, have_hyp), polys in D.items():
            
            filename = f"temp/data_{have_hyp}_{q}_{exps}_{slopes}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    # _ = F.write(", ".join(str(c) for c in poly) + "\n")
                    remainders = [r%q for r in poly[1]]
                    F.write(f"{poly[1]}, {remainders} \n")


def our_jacobi_rules(iso_class, poly, slopes, q, factors, exps, coeffs):
    p = q.smallest_prime_factor()
    R = PolynomialRing(ZZ, 'x')
    
    # irreducible polynomial with NP slopes 1/2-...-1/2
    # non-jacobians are of the form x^6 + pqx^3 + q^3 or x^6 + pqx^3 + q^3
    rule1 = len(factors) == 1 and (slopes == 4) and (coeffs[0] == 0) and (coeffs[1] == 0) and ((coeffs[2] == -p*q) or (coeffs[2] == p*q))
    
    # if slope type 4 and characteristic is even and variety splits as quadratic and quartic
    rule2 = (slopes == 4) and (q % 2 == 0) and (exps == (2, 4))
    
    if iso_class.has_jacobian == -1 or rule1 or rule2:
        return False
    return True

def classify_factors(factors):
    fs = []
    for (fac, exp) in factors:
        if fac.degree() == 1:
            for _ in range(exp//2):
                fs.append(2)
        else:
            for _ in range(exp):
                fs.append(fac.degree()) 
    fs.sort()
    return(tuple(fs))

def get_abc(factors, classification):
    factors_fixed = []
    for (fac, exp) in factors:
        if fac.degree() == 1:
            for i in range(exp//2):
                factors_fixed.append((fac*fac, 2))
        else:
            for i in range(exp):
                factors_fixed.append((fac, fac.degree()))
    factors_fixed.sort(key=lambda x: x[1])
    
    if classification == (6, ):
        a = factors_fixed[0][0][1]
        b = factors_fixed[0][0][2]
        c = factors_fixed[0][0][3]
        
    elif classification == (2, 4):
        a = factors_fixed[0][0][1]
        b = factors_fixed[1][0][1]
        c = factors_fixed[1][0][2]

    elif classification == (2, 2, 2):
        a = factors_fixed[0][0][1]
        b = factors_fixed[1][0][1]
        c = factors_fixed[2][0][1]

    
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
                iso_class = IsogenyClass(poly=polynom)
                factors = R(polynom).factor()
                exps = classify_factors(factors)
                coeffs = get_abc(factors, exps)
                guess = our_jacobi_rules(iso_class, polynom, slope, q, factors, exps, coeffs)
                if have_hyp and guess == True:
                    jacobians.append(polynom)
                elif have_hyp and guess == False:
                    false_positives.append(polynom)
                    
                    false_positives_data[q][exps, slope].append((polynom, coeffs))
                elif not have_hyp and guess==True:
                    new_obstructions.append(polynom)

                    new_obstr_data[q][exps, slope].append((polynom, coeffs))
                
                
                else:
                    known_obstructions.append(polynom)
    
    # return new_obstructions, known_obstructions, jacobians
    print(f"Ratio false positive: {len(false_positives)/len(jacobians)}")
    print(f"Ratio detected negatives: {len(known_obstructions)/len(new_obstructions)}")
    return new_obstr_data, false_positives


def write_new_obstructions_data(data):
    for q, D in data[0].items():
        for (exps, slopes), polys in D.items():
            
            filename = f"temp/new_obs_{q}_{exps}_{slopes}.txt"
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
    
                    



                

