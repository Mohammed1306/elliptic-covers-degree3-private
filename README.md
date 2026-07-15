# Elliptic Covers Degree 3 Private

Private experimental SageMath repository for the degree-3 genus 2 construction in the project **Computing Covers of Elliptic Curves**.

This repository contains the unfinished and experimental work related to the Appendix A degree-3 construction.

The main goal is to work towards a function that takes two elliptic curves as input and constructs a genus 2 curve with explicit degree-3 maps to them.

This repository is separate from the main stable package because the parameter-solving step is not fully resolved yet.

---

## Requirements

This repository is meant to be used with **SageMath**.

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

From the root of this repository, run:

```bash
sage -pip install -e .
```

or, if using the local Sage build:

```bash
~/Summer2026_internship/sage/sage -pip install -e .
```

This installs this repository in editable mode.

To make sure Python is importing from this private repo and not from the main package repo, run:

```bash
~/Summer2026_internship/sage/sage -python -c "import elliptic_covers; print(elliptic_covers.__file__)"
```

The output should point to:

```text
elliptic-covers-degree3-private/elliptic_covers/__init__.py
```

---

## Package Structure

```text
elliptic-covers-degree3-private/
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

### 1. Degree-3 construction from two elliptic curves

```python
genus2_degree3_cover_from_two_elliptic_curves(E1, E2, beta=0, tbeta=0, return_all=False)
```

This is the experimental public degree-3 function.

The function takes two elliptic curves and coordinate translations as input, sets up the Appendix A equations, tries to solve for the parameters `a, b, c, d`, and then constructs a genus 2 curve together with explicit degree-3 maps.

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

---

### 2. Internal Appendix A construction

```python
_genus2_degree3_cover_from_parameters(a, b, c, d)
```

This is the internal forward Appendix A construction.

It starts directly from the parameters `a, b, c, d` and constructs:

- the genus 2 curve `C`
- the first elliptic curve `E1`
- the second elliptic curve `E2`
- explicit degree-3 maps from `C` to `E1` and `E2`

This function is kept private because the final package interface should take elliptic curves as input, not raw parameters.

---

## Morphism Output Format

The repository returns maps using the `ExplicitMorphism` class.

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
    "C_to_E1": ExplicitMorphism(...),
    "C_to_E2": ExplicitMorphism(...),
    "parameters": {
        ...
    },
}
```

The coordinate formulas can be accessed by:

```python
x_formula, y_formula = morphisms["C_to_E1"].formulas()
```

or directly by:

```python
morphisms["C_to_E1"].x
morphisms["C_to_E1"].y
```

The morphisms are explicit Sage rational expressions. They are not yet formal Sage scheme morphism objects.

---

## Experimental Groebner Computation

The repository contains an experimental script based on the symbolic Groebner basis computation suggested by the supervisor:

```text
experiments/run_raghda_groebner.sage
```

Run it with a timeout:

```bash
timeout 600 ~/Summer2026_internship/sage/sage experiments/run_raghda_groebner.sage
```

On the tested machine, Sage created the ideal successfully but fell back to the very slow toy Groebner implementation, and the Groebner basis computation did not finish within the timeout.

This means the symbolic elimination step is still a blocker.

---

## Development Notes

The degree-3 construction requires characteristic different from both `2` and `3`.

The current curve-input function depends on Sage being able to solve the Appendix A parameter system over the current base field.

The translations `beta` and `tbeta` still have to be provided as input.

The main unresolved problem is to find a practical way to solve the polynomial system and recover the parameters `a, b, c, d` from the coefficients of the two elliptic curves.

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

## Current Status

Implemented:

- internal Appendix A parameter construction for degree-3 maps
- experimental curve-input degree-3 construction
- explicit morphism formula objects using `ExplicitMorphism`
- tests for the degree-3 parameter construction
- tests for the degree-3 curve-input construction
- experimental Groebner computation script

Current blockers:

- symbolic Groebner computation does not finish in Sage
- parameter recovery from general elliptic curve coefficients is not fully solved
- the current function only works when Sage can solve the concrete parameter system
- `beta` and `tbeta` are not found automatically yet

Future work:

- find a better way to solve or simplify the Appendix A polynomial system
- investigate whether the degree-16 equations can be solved in a more practical way
- determine how to recover `a, b, c, d` from the elliptic curve coefficients
- automate or reduce the need for choosing `beta` and `tbeta`
- decide which parts should eventually be moved back into the main package