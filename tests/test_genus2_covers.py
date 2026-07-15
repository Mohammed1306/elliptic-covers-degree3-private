from sage.all import *

from sage.all import *

from elliptic_covers import (
    ExplicitMorphism,
    genus2_degree3_cover_from_two_elliptic_curves,
)

from elliptic_covers.genus2_covers import (
    _genus2_degree3_cover_from_parameters,
)


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
    test_degree3_private_forward_construction_valid_case()
    test_degree3_curve_input_function_valid_case()
    test_degree3_curve_input_return_all()

    print("All genus 2 cover tests passed.")