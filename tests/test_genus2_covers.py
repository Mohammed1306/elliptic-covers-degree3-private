from sage.all import *

from elliptic_covers import (
    genus2_cover_from_two_elliptic_curves,
    genus2_cover_from_point,
    genus2_cover_from_two_points,
    genus2_degree3_cover_from_two_elliptic_curves,
)

from elliptic_covers.genus2_covers import (
    ExplicitMorphism,
    _genus2_degree3_cover_from_parameters,
)


sqrt2 = QQbar(2).sqrt()

E = EllipticCurve(QQbar, [0, -8, 0, 8, 0])
F = EllipticCurve(QQbar, [0, 16, 0, 32, 0])

alpha_roots = [
    4 + 2*sqrt2,
    0,
    4 - 2*sqrt2,
]

beta_roots = [
    -8 + 4*sqrt2,
    0,
    -8 - 4*sqrt2,
]


def assert_valid_morphism_output(morphisms):
    assert "C_to_E" in morphisms
    assert "C_to_F" in morphisms
    assert "parameters" in morphisms

    C_to_E = morphisms["C_to_E"]
    C_to_F = morphisms["C_to_F"]

    assert isinstance(C_to_E, ExplicitMorphism)
    assert isinstance(C_to_F, ExplicitMorphism)

    assert C_to_E.source is not None
    assert C_to_E.target is not None
    assert C_to_E.x != 0
    assert C_to_E.y != 0

    assert C_to_F.source is not None
    assert C_to_F.target is not None
    assert C_to_F.x != 0
    assert C_to_F.y != 0


def assert_phase1_parameters(morphisms):
    assert "A" in morphisms["parameters"]
    assert "B" in morphisms["parameters"]
    assert "t1" in morphisms["parameters"]
    assert "t2" in morphisms["parameters"]
    assert "s1" in morphisms["parameters"]
    assert "s2" in morphisms["parameters"]

    assert morphisms["parameters"]["A"] != 0
    assert morphisms["parameters"]["B"] != 0
    assert morphisms["parameters"]["t1"] != 0
    assert morphisms["parameters"]["s1"] != 0


def test_valid_case1():
    C, morphisms = genus2_cover_from_two_elliptic_curves(
        E,
        F,
        alpha_roots,
        beta_roots,
    )

    assert C.genus() == 2
    assert_valid_morphism_output(morphisms)
    assert_phase1_parameters(morphisms)

    h, q = C.hyperelliptic_polynomials()

    assert q == 0
    assert h.degree() == 6
    assert h.discriminant() != 0

    R = h.parent()
    x = R.gen()

    expected_h = 2**53 * (2*x**2 - 1) * (4*x**4 + 12*x**2 + 1)

    assert h == expected_h


def test_valid_case2():
    sqrt5 = QQbar(5).sqrt()

    E2 = EllipticCurve(QQbar, [0, -20, 0, 20, 0])
    F2 = EllipticCurve(QQbar, [0, 40, 0, 320, 0])

    alpha_roots2 = [
        10 + 4*sqrt5,
        0,
        10 - 4*sqrt5,
    ]

    beta_roots2 = [
        -20 - 4*sqrt5,
        0,
        -20 + 4*sqrt5,
    ]

    C, morphisms = genus2_cover_from_two_elliptic_curves(
        E2,
        F2,
        alpha_roots2,
        beta_roots2,
    )

    assert C.genus() == 2
    assert_valid_morphism_output(morphisms)
    assert_phase1_parameters(morphisms)

    h, q = C.hyperelliptic_polynomials()

    assert q == 0
    assert h.degree() == 6
    assert h.discriminant() != 0

    R = h.parent()
    x = R.gen()

    expected_h = (
        -(2**44) * (5**12)
        * (2*x**2 + 2)
        * (16*x**4 + 12*x**2 + 1)
    )

    assert h == expected_h


def test_morphism_formulas_are_nonzero():
    C, morphisms = genus2_cover_from_two_elliptic_curves(
        E,
        F,
        alpha_roots,
        beta_roots,
    )

    assert_valid_morphism_output(morphisms)
    assert_phase1_parameters(morphisms)


def test_invalid_root_length():
    try:
        genus2_cover_from_two_elliptic_curves(
            E,
            F,
            [4 + 2*sqrt2, 0],
            beta_roots,
        )
    except ValueError as e:
        assert str(e) == "alpha_roots and beta_roots must each contain exactly 3 roots"
    else:
        raise AssertionError("Expected ValueError for invalid root list length")


def test_invalid_alpha_root():
    try:
        genus2_cover_from_two_elliptic_curves(
            E,
            F,
            [10, 0, 4 - 2*sqrt2],
            beta_roots,
        )
    except ValueError as e:
        assert str(e) == "alpha_roots must contain roots of f"
    else:
        raise AssertionError("Expected ValueError for invalid alpha root")


def test_invalid_beta_root():
    try:
        genus2_cover_from_two_elliptic_curves(
            E,
            F,
            alpha_roots,
            [10, 0, -8 - 4*sqrt2],
        )
    except ValueError as e:
        assert str(e) == "beta_roots must contain roots of g"
    else:
        raise AssertionError("Expected ValueError for invalid beta root")


def test_cover_from_point_valid_case():
    E_point = EllipticCurve(QQ, [0, -8, 0, 8, 0])
    P = E_point(1, 1)

    C, F_point, morphisms = genus2_cover_from_point(E_point, P)

    assert C.genus() == 2
    assert F_point.genus() == 1

    assert_valid_morphism_output(morphisms)

    H, qC = C.hyperelliptic_polynomials()
    fF, qF = F_point.hyperelliptic_polynomials()

    assert qC == 0
    assert qF == 0

    assert H.degree() == 6
    assert H.discriminant() != 0

    assert fF.degree() == 3
    assert fF.discriminant() != 0

    R = H.parent()
    x = R.gen()

    expected_H = x**6 - 5*x**4 - 5*x**2 + 1

    assert H == expected_H

    R_F = fF.parent()
    w = R_F.gen()

    expected_fF = w**3 - 5*w**2 - 5*w + 1

    assert fF == expected_fF


def test_cover_from_point_rejects_point_at_infinity():
    E_point = EllipticCurve(QQ, [0, -8, 0, 8, 0])
    P = E_point(0)

    try:
        genus2_cover_from_point(E_point, P)
    except ValueError as e:
        assert str(e) == "P must be a finite point on E"
    else:
        raise AssertionError("Expected ValueError for point at infinity")


def test_cover_from_point_rejects_point_on_different_curve():
    E_point = EllipticCurve(QQ, [0, -8, 0, 8, 0])
    E_other = EllipticCurve(QQ, [0, 0, 0, -1, 0])
    P = E_other(0, 0)

    try:
        genus2_cover_from_point(E_point, P)
    except ValueError as e:
        assert str(e) == "P must be a point on E"
    else:
        raise AssertionError("Expected ValueError for point on different curve")


def test_cover_from_point_rejects_two_torsion_point():
    E_point = EllipticCurve(QQ, [0, -8, 0, 8, 0])
    P = E_point(0, 0)

    try:
        genus2_cover_from_point(E_point, P)
    except ValueError as e:
        assert str(e) == "P must not become a 2-torsion point after the coordinate changes"
    else:
        raise AssertionError("Expected ValueError for 2-torsion point")


def test_cover_from_two_points_valid_case_already_opposite():
    E_points = EllipticCurve(QQ, [0, -8, 0, 8, 0])

    P = E_points(1, 1)
    Q = -P

    C, F_points, morphisms = genus2_cover_from_two_points(E_points, P, Q)

    assert C.genus() == 2
    assert F_points.genus() == 1

    assert_valid_morphism_output(morphisms)
    assert "phase_3_translation" in morphisms

    H, qC = C.hyperelliptic_polynomials()
    fF, qF = F_points.hyperelliptic_polynomials()

    assert qC == 0
    assert qF == 0

    assert H.degree() == 6
    assert H.discriminant() != 0

    assert fF.degree() == 3
    assert fF.discriminant() != 0


def test_cover_from_two_points_valid_case_nontrivial_translation():
    E_points = EllipticCurve(QQ, [0, -8, 0, 8, 0])

    T = E_points(1, 1)

    P = T
    Q = -3*T

    assert Q != -P

    C, F_points, morphisms = genus2_cover_from_two_points(E_points, P, Q)

    assert C.genus() == 2
    assert F_points.genus() == 1

    assert_valid_morphism_output(morphisms)
    assert "phase_3_translation" in morphisms

    translation_data = morphisms["phase_3_translation"]

    P_new = translation_data["translated_points"]["P_new"]
    Q_new = translation_data["translated_points"]["Q_new"]

    assert Q_new == -P_new

    H, qC = C.hyperelliptic_polynomials()
    fF, qF = F_points.hyperelliptic_polynomials()

    assert qC == 0
    assert qF == 0

    assert H.degree() == 6
    assert H.discriminant() != 0

    assert fF.degree() == 3
    assert fF.discriminant() != 0
    
def test_degree3_private_forward_construction_valid_case():
    a = QQ(1)
    b = QQ(1)
    c = QQ(-1)
    d = QQ(13) / 16

    C, E1, E2, morphisms = _genus2_degree3_cover_from_parameters(a, b, c, d)

    assert C.genus() == 2
    assert E1.genus() == 1
    assert E2.genus() == 1

    assert "C_to_E1" in morphisms
    assert "C_to_E2" in morphisms
    assert "parameters" in morphisms

    assert isinstance(morphisms["C_to_E1"], ExplicitMorphism)
    assert isinstance(morphisms["C_to_E2"], ExplicitMorphism)

    H, q = C.hyperelliptic_polynomials()

    assert q == 0
    assert H.degree() == 6
    assert H.discriminant() != 0


def test_degree3_curve_input_function_valid_case():
    a = QQ(1)
    b = QQ(1)
    c = QQ(-1)
    d = QQ(13) / 16

    C0, E1, E2, morphisms0 = _genus2_degree3_cover_from_parameters(a, b, c, d)

    C, morphisms = genus2_degree3_cover_from_two_elliptic_curves(
        E1,
        E2,
        beta=QQ(0),
        tbeta=QQ(0),
    )

    assert C.genus() == 2

    assert "C_to_E1" in morphisms
    assert "C_to_E2" in morphisms
    assert "parameters" in morphisms

    assert isinstance(morphisms["C_to_E1"], ExplicitMorphism)
    assert isinstance(morphisms["C_to_E2"], ExplicitMorphism)

    aa = morphisms["parameters"]["a"]
    bb = morphisms["parameters"]["b"]
    cc = morphisms["parameters"]["c"]
    dd = morphisms["parameters"]["d"]

    assert aa**3 + bb**2 != 0
    assert cc**3 + dd**2 != 0
    assert 12*aa*cc + 16*bb*dd == 1


def test_degree3_curve_input_return_all():
    a = QQ(1)
    b = QQ(1)
    c = QQ(-1)
    d = QQ(13) / 16

    C0, E1, E2, morphisms0 = _genus2_degree3_cover_from_parameters(a, b, c, d)

    results = genus2_degree3_cover_from_two_elliptic_curves(
        E1,
        E2,
        beta=QQ(0),
        tbeta=QQ(0),
        return_all=True,
    )

    assert len(results) >= 1

    for C, morphisms in results:
        assert C.genus() == 2
        assert "C_to_E1" in morphisms
        assert "C_to_E2" in morphisms
        assert isinstance(morphisms["C_to_E1"], ExplicitMorphism)
        assert isinstance(morphisms["C_to_E2"], ExplicitMorphism)
        

if __name__ == "__main__":
    test_valid_case1()
    test_valid_case2()
    test_morphism_formulas_are_nonzero()

    test_invalid_root_length()
    test_invalid_alpha_root()
    test_invalid_beta_root()

    test_cover_from_point_valid_case()
    test_cover_from_point_rejects_point_at_infinity()
    test_cover_from_point_rejects_point_on_different_curve()
    test_cover_from_point_rejects_two_torsion_point()

    test_cover_from_two_points_valid_case_already_opposite()
    test_cover_from_two_points_valid_case_nontrivial_translation()

    test_degree3_private_forward_construction_valid_case()
    test_degree3_curve_input_function_valid_case()
    test_degree3_curve_input_return_all()

    print("All genus 2 cover tests passed.")