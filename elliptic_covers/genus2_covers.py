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