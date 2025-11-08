from condcopula import bicop_conditional_sample
import pyvinecopulib as pv
import numpy as np

def test_condcopula():
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

