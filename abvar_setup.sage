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
        fs = []
        for (fac, exp) in factors:
            if fac.degree() == 1:
                for i in range(exp//2):
                    fs.append(2)
            else:
                for i in range(exp):
                    fs.append(fac.degree()) 
        fs.sort()
        exps = tuple(fs)
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

def sort_data(data):
    for q, D in data.items():
        for (exps, slopes, have_hyp), polys in D.items():
            for polynom in polys:
                iso_class = IsogenyClass(poly=polynom)
                guess = our_jacobi_rules
                
