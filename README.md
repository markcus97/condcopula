# condcopula

**Conditional sampling from an arbitrary copula in Python.**

`condcopula` provides a simple interface to conditionally sample a bivariate or a vine copula, with arbitrary conditioning sets. Copulas should be defined using  the [`pyvinecopulib`](https://github.com/vinecopulib/pyvinecopulib) library.

Sampling of a bivariate copula is implemented using the inverse Rosenblatt transform. Vine copula sampling is done using a Markov chain Monte Carlo approach as described and implemented [here](https://doi.org/10.1007/s11222-025-10652-4).

Only continuous variables are supported.



---

## Installation

You can install directly from GitHub using `pip`:

```bash
pip install git+https://github.com/<your-username>/condcopula.git
```

---

## Examples

```python
import pyvinecopulib as pv
import numpy as np
from condcopula import bicop_conditional_sample, vinecop_conditional_sample

## Bivariate copula sampling
bicop = pv.Bicop.from_family(
            family = pv.gaussian,
            rotation = 0,
            parameters = np.array([[0.9]])
        )
    
sample = bicop_conditional_sample(
            bicop,
            value_for_conditional = 0.5,
            conditional_var = 2,
            seed = 12345
        )

print(sample[:10])
[0.37225283 0.41768285 0.64161105 0.57899029 0.45204658 0.42529231
0.54321012 0.34903022 0.57733348 0.75313433]

## Vine copula sampling
# Specify vine copula
bicop11 = pv.Bicop.from_family(
            family = pv.bb1,
            rotation = 0,
            parameters = np.array([[1.0], [2.0]])
        )
bicop12 = pv.Bicop.from_family(
            family = pv.gumbel,
            rotation = 0,
            parameters = np.array([[3.0]])
        )
bicop13 = pv.Bicop.from_family(
            family = pv.gaussian,
            rotation = 0,
            parameters = np.array([[0.9]])
        )
bicop14 = pv.Bicop.from_family(
            family = pv.clayton,
            rotation = 0,
            parameters = np.array([[3.0]])
        )
bicop21 = pv.Bicop.from_family(
            family = pv.clayton,
            rotation = 0,
            parameters = np.array([[1.4]])
        )
bicop22 = pv.Bicop.from_family(
            family = pv.frank,
            rotation = 0,
            parameters = np.array([[2.8]])
        )
bicop23 = pv.Bicop.from_family(
            family = pv.bb7,
            rotation = 0,
            parameters = np.array([[1.0], [1.5]])
        )
bicop31 = pv.Bicop.from_family(
            family = pv.frank,
            rotation = 0,
            parameters = np.array([[1.8]])
        )
bicop32 = pv.Bicop.from_family(
            family = pv.gaussian,
            rotation = 0,
            parameters = np.array([[0.3]])
        )
bicop41 = pv.Bicop.from_family(
            family = pv.gaussian,
            rotation = 0,
            parameters = np.array([[0.05]])
        )
structure_matrix = np.array([
            [2,3,4,5,5],
            [3,4,5,4,0],
            [4,5,3,0,0],
            [5,2,0,0,0],
            [1,0,0,0,0]
        ])
vinecop = pv.Vinecop.from_structure(
            matrix = structure_matrix,
            pair_copulas = [
                [bicop11,bicop12,bicop13,bicop14],
                [bicop21,bicop22,bicop23],
                [bicop31,bicop32],
                [bicop41]
            ]
        )

# Specify values for conditional variables
# We will condition on variables 1, 3, and 5, taking on values
# 0.20, 0.45, 0.78
values_for_conditional = {
            1: 0.20,
            3: 0.45,
            5: 0.78
        }

# Sample from the copula
sample = vinecop_conditional_sample(
            vinecop = vinecop,
            values_for_conditional = values_for_conditional,
            seed = 12345
        )
print(sample)
[[0.31933428 0.24004681]
[0.4696318  0.3716008 ]
[0.18236888 0.38030853]
...
[0.34555428 0.38613813]
[0.51055103 0.34963525]
[0.24754671 0.47006729]]

```


## Attribution

Implementation of conditional vine sampling is a wrapper around the STAN code accompanying:

Hanebeck, A., Şahin, Ö., Havlíčková, P. et al. Sampling from Conditional Distributions of Simplified Vines. Stat Comput 35, 128 (2025). [https://doi.org/10.1007/s11222-025-10652-4](https://doi.org/10.1007/s11222-025-10652-4)

The wrapper itself is a translation into Python of the wrapper the authors [released](https://github.com/ArianeHanebeck/Sampling_Conditional_Vines/) in the R programming language.


