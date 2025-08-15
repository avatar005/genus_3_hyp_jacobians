gen_vs := function(n)
    local q, z, xg, polys, d, G, F, R, x, coeffs, poly, c, i, GAction, OrbitDomain, orbs, reps, coeff, pol;
    q := 2^n;
    F := GF(q);
    z := PrimitiveElement(F);
    x := Indeterminate(F);

    # generate all monic polynomials of degree up to 4
    polys := [UnivariatePolynomial(F, [1], x)];
    for d in [1..4] do
        coeffs := Cartesian(List([1..d], i -> Elements(F)));
        for c in coeffs do
            # add leading coeff
            Add(c, 1);
            pol := Zero(F);
            for i in [1..Length(c)] do
                pol := pol + x^(i-1)*c[i];
            od;
            Add(polys, pol);
        od;
    od;

    G := GeneralLinearGroup(2, F);

    # define the group action
    GAction := function(f, A)
        local a, b, c, d, t, s, result;

        a := A[1][1]; b := A[1][2]; c := A[2][1]; d := A[2][2];
        t := a*x + b; s := c*x + d;

        if f = One(F) then
            result := s ^ 4;
        else
            result := AsPolynomial((s) ^ 4 * Value(f, t / s));
        fi;
        result := result / LeadingCoefficient(result);
        return result;
    end;
    
    # compute the orbits
    orbs := OrbitsDomain(G, polys, GAction);
    reps := List(orbs, o -> o[1]);

    return reps;

end;

