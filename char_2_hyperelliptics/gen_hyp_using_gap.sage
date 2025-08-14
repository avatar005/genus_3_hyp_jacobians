from collections import defaultdict
from sage.all import GF, PolynomialRing, GL, VectorSpace, gcd
from sage.libs.gap.libgap import libgap as gap

def orbit_reps(n):
    '''Calls the gap function defined in "generate_v.g" in order to produce the list of
    q's used in the next step of the algorithm.'''

    # run the gap function
    gap.eval('Read("generate_v.g")')
    result = gap.eval(f'gen_vs({n})')

    # interpret the result as Sage objects
    result = list(result)
    F = GF(2**n, 'a')
    a = F.gen()
    R = PolynomialRing(F, 'x')
    x = R.gen()
    out = []
    for poly in result:
        coeffs = gap.CoefficientsOfUnivariatePolynomial(poly).sage()
        coeffs = [F(str(c)) for c in coeffs]
        f = sum(coeffs[i] * x**i for i in range(len(coeffs)))
        out.append(f)
    return(out)


def generate_hyperelliptics_char2(n, VS=None):
    """
    Generates all isomorphism classes of genus 3 hyperelliptic curves y^2 + v(x)y = u(x)
        over the field k = \mathbb{F}_{2^n}.
    Inputs are n, the degree of the extension over \mathbb{F}_2, and VS is a list of all candidate
        polynomials v(x), as produced by the preceeding function orbit_reps(n).
    The output/ is a generator that produces all pairs v, u determining the desired set of isomorphism classes.
    """
   
    F = GF(2**n, 'a')
    a = F.gen()
    R = PolynomialRing(F, 'x')
    x = R.gen()
    S = R.quotient(x**9, names='x')
    x = S.gen()
    G = GL(2, F)

    def g_action(level, g, P):
        """
        Computes the action of G = GL_2(k) on the ring k[x], as described in [Borodin, May 2025].
        Level is the degree of the action, g is an element of G, and P is an element of k[x].
        The output is the polynomial \psi_{level}(g)(P).
        """
        [[a,b], [c,d]] = g.list()
        A, B = a*x + b, c*x + d
        coeffs = P.list() + [0]*(level+1 - len(P.list()))
        result = sum(p_i * A**i * B**(level - i)
                for i, p_i in enumerate(coeffs) if p_i)
        return result

    def basechange(f):
        """
        Computes the representation of a polynomial f \in k[x] as an element of
            the \mathbb{F}_2 vector space used below.
        Input is a polynomial in the ring k[x]/(x^9).
        Output is a list of length 9*n, corresponding to the basis described below.
        """
        if f == 1:
            return [1] + [0] * (9*n - 1)
        coeffs = f.list()
        coeffs += [0] * (9 - len(coeffs))
        v = [0]*(9*n)
        for exp, c in enumerate(coeffs):
            rep = c.polynomial().list()
            rep += [0] * (n - len(rep))
            for i, bit in enumerate(rep):
                if bit:
                    v[i*9 + exp] = 1
        return v

    def smooth(u,v):
        """
        checks the conditions due to Xarles for whether y^2 + q(x) y = p(x) is hyperelliptic
        """
        if not isinstance(u, int):
            u_lift = S.lift(u)
            u_deriv = u_lift.derivative()
        else:
            u_lift = u
            u_deriv = 0
        return gcd(v, u_deriv**2 + (v.derivative()**2)*u_lift) == 1 and (v.degree() == 4 or u_lift.monomial_coefficient(S.lift(x**7))**2 != u_lift.monomial_coefficient(S.lift(x**8))*v.monomial_coefficient(S.lift(x**3))**2)


    if not VS:
        V_set = orbit_reps(n) # We compute, or use an existing, list of all possible candidate polynomials q(x).
    else:
        V_set = VS

    print("Got v's:", V_set)
    print("Number of v's:", len(V_set))

    basis = [a**i * S(x**j) for i in range(n) for j in range(9)] # this vector space is equivalent as a set to the ring R, but for later purposes
    V = VectorSpace(GF(2), 9*n)                                # it is advantageous to perform computations over \mathbb{F}_2

    # U_set = defaultdict(list) # stores the list of u(x) for each v(x), such that y^2 + u(x) y + v(x) is representative
                           # of a distinct isomorphism class of curves.

    count = 0
    for v in V_set: # iterate over each generated v(x).
        print("v value:", v)
        basis_U = [a**j * x**i * v + a**(2*j) * x**(2*i) for j in range(n) for i in range(5)] # the symbolic basis of the space of polynomials
                                                                                          # of the form r(x)^2 + q(x)r(x)
        new_basis_U = [basechange(f) for f in basis_U]  # computes the appropriate representation
        U = V.subspace([V(base) for base in new_basis_U]) # and here as a subspace of the space of polynomials
        Q = V.quotient(U) # computes the quotient
        candidates = set()
        for coset in Q:
            vec = tuple(Q.lift(coset)) # chooses a vector representative in the quotient space
            u = sum(bit * b for bit, b in zip(vec, basis)) # computes the associated polynomial in U according to the symbolic basis
            candidates.add(u)
        Rq = R(v)
        stab_v = [A for A in G if g_action(4, A, v) == v] # this is the stabilizer of q(x) in the general linear group
        while candidates:
            u = candidates.pop()
            if u == 0:
                u_deg = -1
            else:
                u_deg = max(i for i, c in enumerate(u.list()) if c!=0)
            m = max(2*v.degree(), u_deg)

            if not m in [7,8]:
                continue

            orbit = set()
            for A in stab_v: # for each candidate p(x) we choose one representative from the stab_v orbit
                uS = S(g_action(8, A, u))
                vec_u = basechange(uS)
                cos_u = Q(V(vec_u))
                lift_u = tuple(Q.lift(cos_u))
                u_cand = sum(bit * b for bit, b in zip(lift_u, basis))
                orbit.add(u_cand)
            candidates -= orbit
            if smooth(u, v):
                count += 1
                yield v, u
        

    # hyperelliptics = {}

    # for v in U_set:
    #     u_list = []
    #     for u in U_set[v]:
    #         if smooth(u,v):
    #             u_list.append(u)
    #     hyperelliptics[v] = u_list

    print("Total hyperellitpics:", count)
    # return hyperelliptics

