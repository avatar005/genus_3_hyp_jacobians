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
        poly = R(list(reversed(rec["poly"])))
        factors = poly.factor()
        exps = classify_factors(factors)
        data[q][exps, slope_types[slopes], have_hyp].append(rec["poly"])
    return data

def write_data(data):
    for q, D in data.items():
        for (exps, slopes, have_hyp), polys in D.items():
            
            filename = f"temp/data_{q}_{exps}_{slopes}_{have_hyp}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    _ = F.write(",".join(str(c) for c in poly) + "\n")


def our_jacobi_rules(iso_class, poly):
    if iso_class.has_jacobian == -1:
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
    new_obstr_data = defaultdict(lambda: defaultdict(list))
    for q, D in data.items():
        for (exps, slope, have_hyp), polys in D.items():
            for polynom in polys:
                iso_class = IsogenyClass(poly=polynom)
                guess = our_jacobi_rules(iso_class, polynom)
                if have_hyp:
                    jacobians.append(polynom)
                
                elif not have_hyp and guess==True:
                    new_obstructions.append(polynom)
                    factors = R(polynom).factor()

                    exps = classify_factors(factors)
                    coeffs = get_abc(factors, exps)
                    new_obstr_data[q][exps, slope].append((polynom, coeffs))
                
                
                else:
                    known_obstructions.append(polynom)
    
    # return new_obstructions, known_obstructions, jacobians
    return new_obstr_data


def write_new_obstructions_data(data):
    for q, D in data.items():
        for (exps, slopes), polys in D.items():
            
            filename = f"temp/new_obs_{q}_{exps}_{slopes}.txt"
            with open(filename, "w") as F:
                for poly in polys:
                    _ = F.write(", ".join(str(c) for c in poly) + "\n")
                    



                

