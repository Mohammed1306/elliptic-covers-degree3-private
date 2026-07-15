# Elliptic Covers Sage

SageMath-compatible package for computing covers of elliptic curves.

This repository contains the official internship package development work for the project **Computing Covers of Elliptic Curves**.

The package currently focuses on constructing genus 2 covers of elliptic curves and returning explicit formulas for the corresponding maps.

---

## Requirements

This package is meant to be used with **SageMath**.

Check that SageMath is available:

```bash
sage --version
```

If using the local Sage development build:

```bash
~/Summer2026_internship/sage/sage --version
```

---

## Installation for Development

From the root of the repository, run:

```bash
sage -pip install -e .
```

or, if using the local Sage build:

```bash
~/Summer2026_internship/sage/sage -pip install -e .
```

This installs the package in editable mode, so changes in the source code are immediately reflected when importing the package.

---

## Package Structure

```text
elliptic-covers-sage/
├── elliptic_covers/
│   ├── __init__.py
│   └── genus2_covers.py
├── tests/
│   └── test_genus2_covers.py
├── experiments/
│   └── run_raghda_groebner.sage
├── pyproject.toml
├── README.md
└── .gitignore
```

The main implementation is in:

```text
elliptic_covers/genus2_covers.py
```

The main tests are in:

```text
tests/test_genus2_covers.py
```

Experimental scripts are stored in:

```text
experiments/
```

---

## Implemented Functions

### 1. Construction from two elliptic curves

```python
genus2_cover_from_two_elliptic_curves(E, F, alpha_roots, beta_roots)
```

This function implements the genus 2 construction from Proposition 4 of Howe--Leprévost--Poonen.

Input:

- `E`: elliptic curve in the form `y^2 = f(x)`
- `F`: elliptic curve in the form `y^2 = g(x)`
- `alpha_roots`: ordered roots of `f`
- `beta_roots`: ordered roots of `g`

The order of the roots determines the matching of the 2-torsion points.

Output:

```python
C, morphisms
```

where:

- `C` is the constructed genus 2 hyperelliptic curve
- `morphisms` contains explicit formulas for the maps:
  - `C_to_E`
  - `C_to_F`

Example:

```python
from sage.all import *
from elliptic_covers import genus2_cover_from_two_elliptic_curves

sqrt2 = QQbar(2).sqrt()

E = EllipticCurve(QQbar, [0, -8, 0, 8, 0])
F = EllipticCurve(QQbar, [0, 16, 0, 32, 0])

alpha_roots = [4 + 2*sqrt2, 0, 4 - 2*sqrt2]
beta_roots = [-8 + 4*sqrt2, 0, -8 - 4*sqrt2]

C, morphisms = genus2_cover_from_two_elliptic_curves(
    E,
    F,
    alpha_roots,
    beta_roots,
)

C.genus()
```

---

### 2. Construction from one elliptic curve and one point

```python
genus2_cover_from_point(E, P)
```

This function constructs a genus 2 cover starting from one elliptic curve `E` and one finite point `P` on `E`.

The function performs the following steps:

1. Removes the `a1` and `a3` terms from the elliptic curve equation.
2. Translates the x-coordinate so that the point `P` has x-coordinate `0`.
3. Constructs the genus 2 curve by substituting `x^2` into the transformed cubic.
4. Constructs the complementary genus 1 curve.
5. Returns explicit formulas for the maps.

Output:

```python
C, F, morphisms
```

where:

- `C` is the genus 2 hyperelliptic curve
- `F` is the complementary genus 1 hyperelliptic curve
- `morphisms` contains explicit formulas for:
  - `C_to_E`
  - `C_to_F`

Example:

```python
from sage.all import *
from elliptic_covers import genus2_cover_from_point

E = EllipticCurve(QQ, [0, -8, 0, 8, 0])
P = E(1, 1)

C, F, morphisms = genus2_cover_from_point(E, P)

C.genus()
F.genus()
```

---

### 3. Construction from one elliptic curve and two points

```python
genus2_cover_from_two_points(E, P, Q)
```

This function starts with one elliptic curve `E` and two finite points `P` and `Q` on `E`.

The goal is to translate the points so that they become opposites.

The function looks for a point `T` such that:

```text
2*T = -(P + Q)
```

Then it defines:

```text
P_new = P + T
Q_new = Q + T
```

so that:

```text
Q_new = -P_new
```

After this, it calls the one-point construction on `P_new`.

Output:

```python
C, F, morphisms
```

where:

- `C` is the genus 2 hyperelliptic curve
- `F` is the complementary genus 1 hyperelliptic curve
- `morphisms` contains the formulas inherited from the one-point construction, plus information about the Phase 3 translation.

Example:

```python
from sage.all import *
from elliptic_covers import genus2_cover_from_two_points

E = EllipticCurve(QQ, [0, -8, 0, 8, 0])

P = E(1, 1)
Q = -P

C, F, morphisms = genus2_cover_from_two_points(E, P, Q)

C.genus()
F.genus()
```

---

### 4. Degree-3 construction from two elliptic curves

```python
genus2_degree3_cover_from_two_elliptic_curves(E1, E2, beta=0, tbeta=0, return_all=False)
```

This is the public degree-3 construction.

The function takes two elliptic curves and coordinate translations as input, solves the Appendix A parameter equations for `a, b, c, d`, applies the internal Appendix A construction, and returns a genus 2 curve together with explicit degree-3 maps to the two original elliptic curves.

Input:

- `E1`: elliptic curve in the form `y^2 = x^3 + a2*x^2 + a4*x + a6`
- `E2`: elliptic curve in the form `y^2 = x^3 + ta2*x^2 + ta4*x + ta6`
- `beta`: x-coordinate translation for `E1`
- `tbeta`: x-coordinate translation for `E2`
- `return_all`: if `True`, return all valid parameter solutions found

Output if `return_all=False`:

```python
C, morphisms
```

Output if `return_all=True`:

```python
results
```

where `results` is a list of pairs:

```python
(C, morphisms)
```

The morphism dictionary contains:

- `C_to_E1`
- `C_to_E2`
- `parameters`

Example:

```python
from sage.all import *
from elliptic_covers import genus2_degree3_cover_from_two_elliptic_curves
from elliptic_covers.genus2_covers import _genus2_degree3_cover_from_parameters

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

C.genus()
```

The parameter-based Appendix A function is kept private as:

```python
_genus2_degree3_cover_from_parameters(a, b, c, d)
```

It is used internally by the public curve-input function.

---

## Morphism Output Format

The package returns maps using the `ExplicitMorphism` class.

This class stores:

- the source curve
- the target curve
- the x-coordinate formula
- the y-coordinate formula
- the source coordinate names
- the target coordinate names

Example structure:

```python
morphisms = {
    "C_to_E": ExplicitMorphism(...),
    "C_to_F": ExplicitMorphism(...),
    "parameters": {
        ...
    },
}
```

The coordinate formulas can be accessed by:

```python
x_formula, y_formula = morphisms["C_to_E"].formulas()
```

or directly by:

```python
morphisms["C_to_E"].x
morphisms["C_to_E"].y
```

The morphisms are explicit Sage rational expressions. They are not yet formal Sage scheme morphism objects.

For the two-point construction, the dictionary also contains:

```python
morphisms["phase_3_translation"]
```

This stores the translation point and the translated points used to reduce the two-point case to the one-point case.

For the degree-3 construction, the dictionary contains:

```python
morphisms["parameters"]
```

with the recovered Appendix A parameters and the translated elliptic curve coefficients.

---

## Running Tests

Run the tests with:

```bash
sage -python tests/test_genus2_covers.py
```

or, if using the local Sage build:

```bash
~/Summer2026_internship/sage/sage -python tests/test_genus2_covers.py
```

---

## Running Doctests

Run the doctests with:

```bash
sage -t elliptic_covers/genus2_covers.py
```

or, if using the local Sage build:

```bash
~/Summer2026_internship/sage/sage -t elliptic_covers/genus2_covers.py
```

---

## Experimental Groebner Computation

The repository also contains an experimental script based on the symbolic Groebner basis computation suggested by the supervisor:

```text
experiments/run_raghda_groebner.sage
```

Run it with a timeout:

```bash
timeout 600 ~/Summer2026_internship/sage/sage experiments/run_raghda_groebner.sage
```

On the tested machine, Sage created the ideal successfully but fell back to the very slow toy Groebner implementation, and the Groebner basis computation did not finish within the timeout.

Because of this, the package implementation uses a concrete-input approach: for given elliptic curves and translations `beta` and `tbeta`, the function sets up the Appendix A equations and solves for `a, b, c, d` over the current base field.

---

## Development Notes

The package currently supports the main genus 2 constructions over fields of characteristic different from `2`.

The degree-3 construction requires characteristic different from both `2` and `3`.

The code checks that:

- input points lie on the given elliptic curve
- points are finite
- elliptic curves are in the supported Weierstrass form when required
- relevant polynomials are separable
- the constructed genus 2 polynomial has degree `6`
- the complementary polynomial has degree `3`
- coordinate transformations preserve the j-invariant where appropriate
- degree-3 morphism formulas satisfy the target elliptic curve equations internally

The Phase 3 construction depends on Sage being able to find a point `T` satisfying:

```text
2*T = -(P + Q)
```

over the current base field.

The degree-3 curve-input construction depends on Sage being able to solve the Appendix A parameter system over the current base field.

---

## Notes on HLP Section 4

Section 4 of Howe--Leprévost--Poonen concerns genus 3 curves whose Jacobians are related to products of three elliptic curves.

This is a separate construction from the genus 2 degree-3 construction implemented in this package.

The package currently focuses on genus 2 constructions.

---

## Current Status

Implemented:

- genus 2 construction from two elliptic curves
- genus 2 construction from one elliptic curve and one point
- genus 2 construction from one elliptic curve and two points
- public degree-3 genus 2 construction from two elliptic curves and translations
- internal Appendix A parameter construction for degree-3 maps
- explicit morphism formula objects using `ExplicitMorphism`
- tests for valid and invalid cases
- tests for the degree-3 curve-input construction