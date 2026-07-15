from sage.all import *
from sage.structure.sequence import Sequence

class ExplicitMorphism:
    """
    Container for explicit morphism formulas.

    This is not yet a formal Sage morphism object. It stores the source,
    target, and the rational expressions defining the map.

    The map is represented by:

        source_coordinates -> target_coordinates

    For example:

        C -> E
        (z, Y) -> (x_expr, y_expr)
    """

    def __init__(self, source, target, x_map, y_map, source_coordinates=None, target_coordinates=None):
        self.source = source
        self.target = target
        self.x = x_map
        self.y = y_map
        self.source_coordinates = source_coordinates
        self.target_coordinates = target_coordinates

    def formulas(self):
        """
        Return the coordinate formulas of the morphism.
        """
        return self.x, self.y

    def as_dict(self):
        """
        Return the morphism data as a dictionary.
        """
        return {
            "source": self.source,
            "target": self.target,
            "x": self.x,
            "y": self.y,
            "source_coordinates": self.source_coordinates,
            "target_coordinates": self.target_coordinates,
        }

    def __repr__(self):
        return (
            f"ExplicitMorphism("
            f"source_coordinates={self.source_coordinates}, "
            f"target_coordinates={self.target_coordinates}, "
            f"x={self.x}, "
            f"y={self.y})"
        )


def genus2_cover_from_two_elliptic_curves(E, F, alpha_roots, beta_roots):
    """
    Construct the genus 2 curve C: y^2 = h(x) from Proposition 4.

    INPUT:

    - ``E`` -- elliptic curve in the form y^2 = f(x)
    - ``F`` -- elliptic curve in the form y^2 = g(x)
    - ``alpha_roots`` -- ordered list [alpha_1, alpha_2, alpha_3] of the roots of f
    - ``beta_roots`` -- ordered list [beta_1, beta_2, beta_3] of the roots of g

    OUTPUT:

    - ``C`` -- the genus 2 hyperelliptic curve
    - ``morphisms`` -- dictionary containing explicit formulas for C -> E and C -> F
    """

    

    K = Sequence(
        list(E.a_invariants())
        + list(F.a_invariants())
        + list(alpha_roots)
        + list(beta_roots)
    ).universe()
    
    if K.characteristic() == 2:
        raise ValueError("the base field must have characteristic different from 2")

    R = PolynomialRing(K, "x")
    x = R.gen()

    # Sage elliptic curve equation:
    # y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    a1_E, a2_E, a3_E, a4_E, a6_E = [K(c) for c in E.a_invariants()]
    a1_F, a2_F, a3_F, a4_F, a6_F = [K(c) for c in F.a_invariants()]

    # We only support curves of the form y^2 = cubic(x)
    # This means a1 = 0 and a3 = 0
    if a1_E != 0 or a3_E != 0:
        raise ValueError("E must be in the form y^2 = f(x)")
    
    if a1_F != 0 or a3_F != 0:
        raise ValueError("F must be in the form y^2 = g(x)")

    if len(alpha_roots) != 3 or len(beta_roots) != 3:
        raise ValueError("alpha_roots and beta_roots must each contain exactly 3 roots")


    f = x**3 + a2_E*x**2 + a4_E*x + a6_E
    g = x**3 + a2_F*x**2 + a4_F*x + a6_F

    if f.discriminant() == 0:
        raise ValueError("f must be separable")

    if g.discriminant() == 0:
        raise ValueError("g must be separable")

    alpha_roots = [K(alpha) for alpha in alpha_roots]
    beta_roots = [K(beta) for beta in beta_roots]

    for alpha in alpha_roots:
        if f(alpha) != 0:
            raise ValueError("alpha_roots must contain roots of f")

    for beta in beta_roots:
        if g(beta) != 0:
            raise ValueError("beta_roots must contain roots of g")

    alpha1, alpha2, alpha3 = alpha_roots
    beta1, beta2, beta3 = beta_roots

    # Quantities from Proposition 4
    a1 = (
        (alpha3 - alpha2)**2 / (beta3 - beta2)
        + (alpha2 - alpha1)**2 / (beta2 - beta1)
        + (alpha1 - alpha3)**2 / (beta1 - beta3)
    )

    b1 = (
        (beta3 - beta2)**2 / (alpha3 - alpha2)
        + (beta2 - beta1)**2 / (alpha2 - alpha1)
        + (beta1 - beta3)**2 / (alpha1 - alpha3)
    )

    a2 = (
        alpha1*(beta3 - beta2)
        + alpha2*(beta1 - beta3)
        + alpha3*(beta2 - beta1)
    )

    b2 = (
        beta1*(alpha3 - alpha2)
        + beta2*(alpha1 - alpha3)
        + beta3*(alpha2 - alpha1)
    )

    if a1 == 0 or b1 == 0 or a2 == 0 or b2 == 0:
        raise ValueError("a1, b1, a2 and b2 must be nonzero")

    Delta_f = f.discriminant()
    Delta_g = g.discriminant()

    A = Delta_g * a1 / a2
    B = Delta_f * b1 / b2

    if A == 0 or B == 0:
        raise ValueError("A and B must be nonzero")

    factor1 = (
        A*(alpha2 - alpha1)*(alpha1 - alpha3)*x**2
        + B*(beta2 - beta1)*(beta1 - beta3)
    )

    factor2 = (
        A*(alpha3 - alpha2)*(alpha2 - alpha1)*x**2
        + B*(beta3 - beta2)*(beta2 - beta1)
    )

    factor3 = (
        A*(alpha1 - alpha3)*(alpha3 - alpha2)*x**2
        + B*(beta1 - beta3)*(beta3 - beta2)
    )

    h = -(factor1 * factor2 * factor3)

    if h.degree() != 6:
        raise ValueError("the resulting polynomial h must have degree 6")

    if h.discriminant() == 0:
        raise ValueError("the resulting polynomial h must be separable")

    C = HyperellipticCurve(h)

    # ------------------------------------------------------------------
    # Morphism formulas from Proposition 4.
    #
    # We use C-coordinates:
    #
    #     C : yC^2 = h(xC)
    #
    # The map C -> F is:
    #
    #     x_F = t1*xC^2 + t2
    #     y_F = (Delta_f / B^3)*yC
    #
    # The map C -> E is:
    #
    #     x_E = s1/xC^2 + s2
    #     y_E = (Delta_g / A^3)*(yC/xC^3)
    #
    # For now these are returned as exact rational expressions.
    # ------------------------------------------------------------------

    t1 = -(A / B) * (b2 / b1)

    t2 = (
        beta1*(beta3 - beta2)**2 / (alpha3 - alpha2)
        + beta2*(beta1 - beta3)**2 / (alpha1 - alpha3)
        + beta3*(beta2 - beta1)**2 / (alpha2 - alpha1)
    ) / b1

    # s1 and s2 are obtained from t1 and t2 by exchanging alpha and beta.
    s1 = -(B / A) * (a2 / a1)

    s2 = (
        alpha1*(alpha3 - alpha2)**2 / (beta3 - beta2)
        + alpha2*(alpha1 - alpha3)**2 / (beta1 - beta3)
        + alpha3*(alpha2 - alpha1)**2 / (beta2 - beta1)
    ) / a1

    S = PolynomialRing(K, ("xC", "yC"))
    xC, yC = S.gens()
    L = S.fraction_field()

    xC = L(xC)
    yC = L(yC)

    C_to_F = ExplicitMorphism(
        source=C,
        target=F.change_ring(K),
        x_map=t1*xC**2 + t2,
        y_map=(Delta_f / B**3)*yC,
        source_coordinates=("xC", "yC"),
        target_coordinates=("xF", "yF"),
    )

    C_to_E = ExplicitMorphism(
        source=C,
        target=E.change_ring(K),
        x_map=s1 / xC**2 + s2,
        y_map=(Delta_g / A**3) * (yC / xC**3),
        source_coordinates=("xC", "yC"),
        target_coordinates=("xE", "yE"),
    )

    morphisms = {
        "C_to_E": C_to_E,
        "C_to_F": C_to_F,
        "parameters": {
            "A": A,
            "B": B,
            "t1": t1,
            "t2": t2,
            "s1": s1,
            "s2": s2,
        },
    }

    return C, morphisms


def genus2_cover_from_point(E, P):
    """
    Construct a genus 2 cover from one elliptic curve and one point.

    INPUT:

    - ``E`` -- an elliptic curve
    - ``P`` -- a finite point on ``E``

    OUTPUT:

    - ``C`` -- the genus 2 hyperelliptic curve
    - ``F`` -- the complementary genus 1 hyperelliptic curve
    - ``morphisms`` -- dictionary containing explicit formulas for C -> E and C -> F

    EXAMPLES::

        sage: from elliptic_covers import genus2_cover_from_point
        sage: E = EllipticCurve(QQ, [0, -8, 0, 8, 0])
        sage: P = E(1, 1)
        sage: C, F, morphisms = genus2_cover_from_point(E, P)
        sage: C.genus()
        2
        sage: F.genus()
        1
    """
    if P.curve() != E:
        raise ValueError("P must be a point on E")

    if P.is_zero():
        raise ValueError("P must be a finite point on E")

    a1, a2, a3, a4, a6 = E.a_invariants()

    K = Sequence(
        list(E.a_invariants()) + [P[0], P[1]]
    ).universe()

    if K.characteristic() == 2:
        raise ValueError("the base field must have characteristic different from 2")

    a1 = K(a1)
    a2 = K(a2)
    a3 = K(a3)
    a4 = K(a4)
    a6 = K(a6)

    xP = K(P[0])
    yP = K(P[1])

    E_K = E.change_ring(K)

    R = PolynomialRing(K, "x")
    x = R.gen()

    # ------------------------------------------------------------------
    # Step 1: Kill a1 and a3.
    #
    # Original equation:
    #
    #     y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    #
    # New coordinates:
    #
    #     X = x
    #     Y = y + (a1*x + a3)/2
    #
    # New equation:
    #
    #     Y^2 = X^3 + (b2/4)*X^2 + (b4/2)*X + b6/4
    # ------------------------------------------------------------------

    b2 = a1**2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3**2 + 4*a6

    completed_polynomial = (
        x**3
        + (b2 / 4)*x**2
        + (b4 / 2)*x
        + (b6 / 4)
    )

    XP = xP
    YP = yP + (a1*xP + a3) / 2

    E_completed = EllipticCurve(K, [
        0,
        b2 / 4,
        0,
        b4 / 2,
        b6 / 4,
    ])

    if E_completed.j_invariant() != E_K.j_invariant():
        raise ValueError("j-invariant changed after killing a1 and a3")

    # ------------------------------------------------------------------
    # Step 2: Translate P so that its x-coordinate becomes 0.
    #
    # New coordinate:
    #
    #     u = X - XP
    #
    # so:
    #
    #     X = u + XP
    #
    # In the code, the variable x represents u.
    # ------------------------------------------------------------------

    transformed_polynomial = completed_polynomial(x + XP)

    if transformed_polynomial(0) != YP**2:
        raise ValueError("coordinate translation did not move P correctly")

    if transformed_polynomial(0) == 0:
        raise ValueError("P must not become a 2-torsion point after the coordinate changes")

    c2 = transformed_polynomial[2]
    c1 = transformed_polynomial[1]
    c0 = transformed_polynomial[0]

    E_transformed = EllipticCurve(K, [
        0,
        c2,
        0,
        c1,
        c0,
    ])

    if E_transformed.j_invariant() != E_K.j_invariant():
        raise ValueError("j-invariant changed after translating the x-coordinate")

    # ------------------------------------------------------------------
    # Step 3: Construct the genus 2 curve.
    #
    # Transformed elliptic curve:
    #
    #     Y^2 = f(u)
    #
    # Substitute:
    #
    #     u = z^2
    #
    # Then:
    #
    #     C : Y^2 = f(z^2)
    #
    # In the code, the variable x represents z here.
    # ------------------------------------------------------------------

    H = transformed_polynomial(x**2)

    if H.degree() != 6:
        raise ValueError("the resulting polynomial H must have degree 6")

    if H.discriminant() == 0:
        raise ValueError("the resulting polynomial H must be separable")

    C = HyperellipticCurve(H)

    # ------------------------------------------------------------------
    # Step 4: Construct the complementary curve F.
    #
    # If:
    #
    #     f(u) = u^3 + c2*u^2 + c1*u + c0
    #
    # then the complementary curve is:
    #
    #     F : v^2 = 1 + c2*w + c1*w^2 + c0*w^3
    #
    # In the code, the variable x represents w here.
    # ------------------------------------------------------------------

    complementary_polynomial = 1 + c2*x + c1*x**2 + c0*x**3

    if complementary_polynomial.degree() != 3:
        raise ValueError("the complementary polynomial must have degree 3")

    if complementary_polynomial.discriminant() == 0:
        raise ValueError("the complementary polynomial must be separable")

    F = HyperellipticCurve(complementary_polynomial)

    # ------------------------------------------------------------------
    # Morphism formulas.
    #
    # C-coordinates:
    #
    #     C : Y^2 = H(z)
    #
    # Map C -> original E:
    #
    #     u = z^2
    #     X = u + XP = z^2 + XP
    #     x = X
    #     y = Y - (a1*X + a3)/2
    #
    # Therefore:
    #
    #     x_E = z^2 + XP
    #     y_E = Y - (a1*(z^2 + XP) + a3)/2
    #
    # Map C -> F:
    #
    #     w = 1/u = 1/z^2
    #     v = Y/z^3
    # ------------------------------------------------------------------

    S = PolynomialRing(K, ("z", "Y"))
    z, Y = S.gens()
    L = S.fraction_field()

    z = L(z)
    Y = L(Y)

    C_to_E = ExplicitMorphism(
        source=C,
        target=E_K,
        x_map=z**2 + XP,
        y_map=Y - (a1*(z**2 + XP) + a3) / 2,
        source_coordinates=("z", "Y"),
        target_coordinates=("x", "y"),
    )

    C_to_F = ExplicitMorphism(
        source=C,
        target=F,
        x_map=1 / z**2,
        y_map=Y / z**3,
        source_coordinates=("z", "Y"),
        target_coordinates=("w", "v"),
    )

    morphisms = {
        "C_to_E": C_to_E,
        "C_to_F": C_to_F,
        "parameters": {
            "XP": XP,
            "YP": YP,
            "a_invariants": (a1, a2, a3, a4, a6),
            "b_invariants": (b2, b4, b6),
            "completed_polynomial": completed_polynomial,
            "transformed_polynomial": transformed_polynomial,
            "cover_polynomial": H,
            "complementary_polynomial": complementary_polynomial,
        },
    }

    return C, F, morphisms

def genus2_cover_from_two_points(E, P, Q):
    """
    Construct a genus 2 cover from one elliptic curve and two points.

    INPUT:

    - ``E`` -- an elliptic curve
    - ``P`` -- a finite point on ``E``
    - ``Q`` -- a finite point on ``E``

    OUTPUT:

    - ``C`` -- the genus 2 hyperelliptic curve
    - ``F`` -- the complementary genus 1 hyperelliptic curve
    - ``morphisms`` -- dictionary containing explicit formulas for C -> E and C -> F

    The function translates P and Q by a point T such that the translated
    points become opposites.

    We look for T satisfying:

        2*T = -(P + Q)

    Then:

        P_new = P + T
        Q_new = Q + T

    so:

        Q_new = -P_new

    Finally, we apply the one-point construction to P_new.
    """
    if P.curve() != E:
        raise ValueError("P must be a point on E")

    if Q.curve() != E:
        raise ValueError("Q must be a point on E")

    if P.is_zero():
        raise ValueError("P must be a finite point on E")

    if Q.is_zero():
        raise ValueError("Q must be a finite point on E")

    K = Sequence(
        list(E.a_invariants()) + [P[0], P[1], Q[0], Q[1]]
    ).universe()

    if K.characteristic() == 2:
        raise ValueError("the base field must have characteristic different from 2")

    E_K = E.change_ring(K)
    P_K = E_K(P)
    Q_K = E_K(Q)

    # We want to translate both points by T so that:
    #
    #     P_new = P + T
    #     Q_new = Q + T
    #
    # and:
    #
    #     Q_new = -P_new
    #
    # This means:
    #
    #     P + T + Q + T = 0
    #
    # so:
    #
    #     2*T = -(P + Q)

    target = -(P_K + Q_K)

    half_points = target.division_points(2)

    if len(half_points) == 0:
        raise ValueError(
            "could not find a point T such that 2*T = -(P + Q) over the current base field"
        )

    # Choose one valid half-point.
    # Different choices of T may lead to different but related constructions.
    T = half_points[0]

    P_new = P_K + T
    Q_new = Q_K + T

    if P_new.is_zero():
        raise ValueError("after translation, P became the point at infinity")

    if Q_new.is_zero():
        raise ValueError("after translation, Q became the point at infinity")

    if P_new + Q_new != E_K(0):
        raise ValueError("translation failed: the two points did not become opposites")

    if Q_new != -P_new:
        raise ValueError("translation failed: Q_new is not the negative of P_new")

    # Apply the one-point construction to the translated point.
    C, F, morphisms = genus2_cover_from_point(E_K, P_new)

    # Add Phase 3 information to the morphism dictionary.
    # The actual maps C -> E and C -> F are already produced by
    # genus2_cover_from_point. Here we only record how the chosen point
    # P_new was obtained from the original two points P and Q.
    morphisms["phase_3_translation"] = {
        "original_points": {
            "P": P_K,
            "Q": Q_K,
        },
        "target": target,
        "translation_point": T,
        "translated_points": {
            "P_new": P_new,
            "Q_new": Q_new,
        },
        "relation": {
            "half_point_equation": "2*T = -(P + Q)",
            "opposite_points": "Q_new = -P_new",
        },
    }

    return C, F, morphisms

def _genus2_degree3_cover_from_parameters(a, b, c, d):
    """
    Internal Appendix A construction.

    This is the forward construction from parameters a,b,c,d.
    It is kept private because the public package function should take
    elliptic curves as input.
    """
    K = Sequence([a, b, c, d]).universe()

    if K.characteristic() in (2, 3):
        raise ValueError("the base field must have characteristic neither 2 nor 3")

    a = K(a)
    b = K(b)
    c = K(c)
    d = K(d)

    if 12*a*c + 16*b*d != 1:
        raise ValueError("parameters must satisfy 12*a*c + 16*b*d = 1")

    Delta1 = a**3 + b**2
    Delta2 = c**3 + d**2

    if Delta1 == 0:
        raise ValueError("a^3 + b^2 must be nonzero")

    if Delta2 == 0:
        raise ValueError("c^3 + d^2 must be nonzero")

    R = PolynomialRing(K, "x")
    x = R.gen()

    f = (
        (x**3 + 3*a*x + 2*b)
        * (2*d*x**3 + 3*c*x**2 + 1)
    )

    if f.degree() != 6:
        raise ValueError("the genus 2 polynomial must have degree 6")

    if f.discriminant() == 0:
        raise ValueError("the genus 2 polynomial must be separable")

    C = HyperellipticCurve(f)

    A1 = 12*(2*a**2*d - b*c)
    B1 = 12*(16*a*d**2 + 3*c**2)*Delta1
    C1 = 512*Delta1**2*d**3

    f1 = x**3 + A1*x**2 + B1*x + C1

    if f1.discriminant() == 0:
        raise ValueError("the first elliptic curve polynomial must be separable")

    E1 = EllipticCurve(K, [0, A1, 0, B1, C1])

    A2 = 12*(2*b*c**2 - a*d)
    B2 = 12*(16*b**2*c + 3*a**2)*Delta2
    C2 = 512*Delta2**2*b**3

    f2 = x**3 + A2*x**2 + B2*x + C2

    if f2.discriminant() == 0:
        raise ValueError("the second elliptic curve polynomial must be separable")

    E2 = EllipticCurve(K, [0, A2, 0, B2, C2])

    # Verify the morphism identities internally.
    L = R.fraction_field()
    xL = L(x)
    fL = L(f)

    denominator1 = xL**3 + 3*a*xL + 2*b
    denominator2 = 2*d*xL**3 + 3*c*xL**2 + 1

    u1_check = 12*Delta1*(-2*d*xL + c) / denominator1
    v1_check = Delta1*(16*d*xL**3 - 12*c*xL**2 - 1) / denominator1**2

    u2_check = 12*Delta2*xL**2*(a*xL - 2*b) / denominator2
    v2_check = Delta2*(xL**3 + 12*a*xL - 16*b) / denominator2**2

    check1 = (
        u1_check**3
        + A1*u1_check**2
        + B1*u1_check
        + C1
        - fL*v1_check**2
    )

    check2 = (
        u2_check**3
        + A2*u2_check**2
        + B2*u2_check
        + C2
        - fL*v2_check**2
    )

    if check1 != 0:
        raise ValueError(
            "the first degree-3 morphism formula does not satisfy the target curve equation"
        )

    if check2 != 0:
        raise ValueError(
            "the second degree-3 morphism formula does not satisfy the target curve equation"
        )

    S = PolynomialRing(K, ("xC", "yC"))
    xC, yC = S.gens()
    M = S.fraction_field()

    xC = M(xC)
    yC = M(yC)

    denominator1 = xC**3 + 3*a*xC + 2*b
    denominator2 = 2*d*xC**3 + 3*c*xC**2 + 1

    u1 = 12*Delta1*(-2*d*xC + c) / denominator1
    v1 = Delta1*(16*d*xC**3 - 12*c*xC**2 - 1) / denominator1**2

    u2 = 12*Delta2*xC**2*(a*xC - 2*b) / denominator2
    v2 = Delta2*(xC**3 + 12*a*xC - 16*b) / denominator2**2

    C_to_E1 = ExplicitMorphism(
        source=C,
        target=E1,
        x_map=u1,
        y_map=yC*v1,
        source_coordinates=("xC", "yC"),
        target_coordinates=("x1", "y1"),
    )

    C_to_E2 = ExplicitMorphism(
        source=C,
        target=E2,
        x_map=u2,
        y_map=yC*v2,
        source_coordinates=("xC", "yC"),
        target_coordinates=("x2", "y2"),
    )

    morphisms = {
        "C_to_E1": C_to_E1,
        "C_to_E2": C_to_E2,
        "parameters": {
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "Delta1": Delta1,
            "Delta2": Delta2,
            "genus2_polynomial": f,
            "elliptic_polynomial_1": f1,
            "elliptic_polynomial_2": f2,
            "A1": A1,
            "B1": B1,
            "C1": C1,
            "A2": A2,
            "B2": B2,
            "C2": C2,
            "u1": u1,
            "v1": v1,
            "u2": u2,
            "v2": v2,
        },
    }

    return C, E1, E2, morphisms


def _elliptic_curve_cubic_coefficients(E, K):
    """
    Return a2, a4, a6 for an elliptic curve of the form

        y^2 = x^3 + a2*x^2 + a4*x + a6.
    """
    a1, a2, a3, a4, a6 = E.a_invariants()

    a1 = K(a1)
    a2 = K(a2)
    a3 = K(a3)
    a4 = K(a4)
    a6 = K(a6)

    if a1 != 0 or a3 != 0:
        raise ValueError(
            "elliptic curves must be in the form y^2 = x^3 + a2*x^2 + a4*x + a6"
        )

    return a2, a4, a6


def _translated_cubic_coefficients(a2, a4, a6, beta):
    """
    Translate x_old = x_new + beta in

        x_old^3 + a2*x_old^2 + a4*x_old + a6.

    The result is

        x_new^3 + A*x_new^2 + B*x_new + C.
    """
    A = a2 + 3*beta
    B = a4 + 2*a2*beta + 3*beta**2
    C = a6 + a4*beta + a2*beta**2 + beta**3

    return A, B, C


def _solution_value(solution, variable, K):
    """
    Read a value from a Sage variety solution dictionary.
    """
    if variable in solution:
        return K(solution[variable])

    variable_name = str(variable)

    if variable_name in solution:
        return K(solution[variable_name])

    raise KeyError(f"could not find value for variable {variable}")


def _solve_degree3_parameters_from_translated_coefficients(A1, B1, C1, A2, B2, C2):
    """
    Solve Raghda's concrete Appendix A inverse system.

    Input curves after translation are assumed to be:

        E1 : y^2 = x^3 + A1*x^2 + B1*x + C1
        E2 : y^2 = x^3 + A2*x^2 + B2*x + C2

    The unknowns are a,b,c,d.
    """
    K = Sequence([A1, B1, C1, A2, B2, C2]).universe()

    A1 = K(A1)
    B1 = K(B1)
    C1 = K(C1)
    A2 = K(A2)
    B2 = K(B2)
    C2 = K(C2)

    if K.characteristic() in (2, 3):
        raise ValueError("the base field must have characteristic neither 2 nor 3")

    P = PolynomialRing(K, ("a", "b", "c", "d"))
    a, b, c, d = P.gens()

    Delta1 = a**3 + b**2
    Delta2 = c**3 + d**2

    equations = [
        12*(2*a**2*d - b*c) - A1,
        12*(16*a*d**2 + 3*c**2)*Delta1 - B1,
        512*Delta1**2*d**3 - C1,

        12*(2*b*c**2 - a*d) - A2,
        12*(16*b**2*c + 3*a**2)*Delta2 - B2,
        512*Delta2**2*b**3 - C2,

        12*a*c + 16*b*d - 1,
    ]

    I = P.ideal(equations)

    try:
        raw_solutions = I.variety(ring=K)
    except Exception as error:
        raise NotImplementedError(
            "The degree-3 parameter system could not be solved over the current base field."
        ) from error

    valid_solutions = []

    for solution in raw_solutions:
        aa = _solution_value(solution, a, K)
        bb = _solution_value(solution, b, K)
        cc = _solution_value(solution, c, K)
        dd = _solution_value(solution, d, K)

        if aa**3 + bb**2 == 0:
            continue

        if cc**3 + dd**2 == 0:
            continue

        if 12*aa*cc + 16*bb*dd != 1:
            continue

        valid_solutions.append((aa, bb, cc, dd))

    if len(valid_solutions) == 0:
        raise ValueError("no valid degree-3 parameters a, b, c, d were found")

    return valid_solutions


def genus2_degree3_cover_from_two_elliptic_curves(E1, E2, beta=0, tbeta=0, return_all=False):
    """
    Construct a genus 2 curve with degree-3 maps to two elliptic curves.

    This is the public degree-3 analogue of the Phase 1 function.

    INPUT:

    - ``E1`` -- elliptic curve in the form
          y^2 = x^3 + a2*x^2 + a4*x + a6

    - ``E2`` -- elliptic curve in the form
          y^2 = x^3 + ta2*x^2 + ta4*x + ta6

    - ``beta`` -- x-coordinate translation for E1

    - ``tbeta`` -- x-coordinate translation for E2

    - ``return_all`` -- if True, return all valid solutions

    OUTPUT:

    If ``return_all`` is False:

        C, morphisms

    If ``return_all`` is True:

        list of pairs (C, morphisms)

    The function translates the input curves, solves the Appendix A inverse
    equations for a,b,c,d, applies the internal forward construction, and then
    adjusts the morphism formulas so that the targets are the original input
    elliptic curves.
    """
    K = Sequence(
        list(E1.a_invariants())
        + list(E2.a_invariants())
        + [beta, tbeta]
    ).universe()

    if K.characteristic() in (2, 3):
        raise ValueError("the base field must have characteristic neither 2 nor 3")

    beta = K(beta)
    tbeta = K(tbeta)

    E1_K = E1.change_ring(K)
    E2_K = E2.change_ring(K)

    a2, a4, a6 = _elliptic_curve_cubic_coefficients(E1_K, K)
    ta2, ta4, ta6 = _elliptic_curve_cubic_coefficients(E2_K, K)

    A1, B1, C1 = _translated_cubic_coefficients(a2, a4, a6, beta)
    A2, B2, C2 = _translated_cubic_coefficients(ta2, ta4, ta6, tbeta)

    parameter_solutions = _solve_degree3_parameters_from_translated_coefficients(
        A1,
        B1,
        C1,
        A2,
        B2,
        C2,
    )

    results = []

    for a, b, c, d in parameter_solutions:
        C, appendix_E1, appendix_E2, appendix_morphisms = _genus2_degree3_cover_from_parameters(
            a,
            b,
            c,
            d,
        )

        C_to_translated_E1 = appendix_morphisms["C_to_E1"]
        C_to_translated_E2 = appendix_morphisms["C_to_E2"]

        C_to_E1 = ExplicitMorphism(
            source=C,
            target=E1_K,
            x_map=C_to_translated_E1.x + beta,
            y_map=C_to_translated_E1.y,
            source_coordinates=C_to_translated_E1.source_coordinates,
            target_coordinates=("x_E1", "y_E1"),
        )

        C_to_E2 = ExplicitMorphism(
            source=C,
            target=E2_K,
            x_map=C_to_translated_E2.x + tbeta,
            y_map=C_to_translated_E2.y,
            source_coordinates=C_to_translated_E2.source_coordinates,
            target_coordinates=("x_E2", "y_E2"),
        )

        morphisms = {
            "C_to_E1": C_to_E1,
            "C_to_E2": C_to_E2,
            "parameters": {
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "beta": beta,
                "tbeta": tbeta,
                "translated_coefficients_E1": (A1, B1, C1),
                "translated_coefficients_E2": (A2, B2, C2),
                "appendix_E1": appendix_E1,
                "appendix_E2": appendix_E2,
                "appendix_morphisms": appendix_morphisms,
            },
        }

        results.append((C, morphisms))

    if return_all:
        return results

    return results[0]