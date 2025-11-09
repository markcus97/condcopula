# ======================================================
# Dependencies
# ======================================================
from typing import Union, Dict, Optional
from pathlib import Path

from pydantic import (
      BaseModel,
      StrictFloat,
      StrictInt,
      conint,
      model_validator
)
import pyvinecopulib as pv
import numpy as np

import cmdstanpy
from cmdstanpy import CmdStanModel, install_cmdstan
cmdstanpy.disable_logging()

this_dir = Path(__file__).resolve().parent
stan1_file = this_dir / "stan/STAN.stan"
stan1_file = stan1_file.resolve()
stan2_file = this_dir / "stan/STAN2.stan"
stan2_file = stan2_file.resolve()
try:
    stan1_model = CmdStanModel(stan_file=stan1_file)
    stan2_model = CmdStanModel(stan_file=stan2_file)
except ValueError as e:
    if str(e).startswith("No CmdStan installation found"):
        print("No CmdStan installation found, attempting to install it.")
        _ = install_cmdstan(compiler=True)
        stan1_model = CmdStanModel(stan_file=stan1_file)
        stan2_model = CmdStanModel(stan_file=stan2_file)
    else:
        print("The two STAN files are: \n"
              f"{stan1_file}\n"
              f"{stan2_file}\n"
              "You need to be able to run: import cmdstanpy; cmdstanpy.CmdStanModel(stan_file=stan_file_path).\n" \
              "Ensure this works and try importing the library again.")
        raise



# ======================================================
# Constants
# ======================================================
_NUM_CHAINS = 4
_ADAPT_DELTA = 0.8
_MAX_TREEDEPTH = 10



# ======================================================
# Input validation
# ======================================================
class _vinecop_validation(BaseModel,arbitrary_types_allowed=True):
    vinecop: pv.Vinecop
    values_for_conditional: Dict[StrictInt, Union[StrictInt,StrictFloat, np.integer, np.floating]]
    num_samples: conint(ge=1) = 1000 # type: ignore
    burnin: StrictInt = 1000
    thin: StrictInt = 10

    @model_validator(mode="after")
    def check_values_for_conditional(self):
        num_conditional_variables = len(self.values_for_conditional)
        copula_dim = self.vinecop.dim

        if num_conditional_variables == 0:
            raise ValueError("Dictionary values_for_conditional should have at least one entry.")

        if num_conditional_variables >= copula_dim:
            raise ValueError(
                f"The number of conditional variables (currently {num_conditional_variables})\n"
                f"should be smaller than the copula dimension {copula_dim}."
            )

        for key, value in self.values_for_conditional.items():
            if key <= 0:
                raise ValueError("Keys of values_for_conditional should be integers corresponding"
                                 "to indices of conditional variables. The lowest possible index is 1.")
            if key > copula_dim:
                raise ValueError("Keys of values_for_conditional should be integers corresponding"
                                 f"to indices of conditional variables. The index cannot be greater than copula dimension ({copula_dim}).")
            
            if (value < 0) or (value > 1):
                raise ValueError("Values in values_for_conditional should be numbers between 0 and 1.")

        return self

    @model_validator(mode="after")
    def check_var_types(self):
        if not all(x == "c" for x in self.vinecop.var_types):
            raise ValueError("Copula contains discrete variables - this is not supported.")
        
        return self



class _bicop_validation(BaseModel,arbitrary_types_allowed=True):
    bicop: pv.Bicop
    value_for_conditional: Union[StrictInt,StrictFloat, np.integer, np.floating]
    conditional_var: Union[StrictInt,StrictFloat, np.integer, np.floating]
    num_samples: conint(ge=1) = 1000 # type: ignore
    seed: Optional[StrictInt] = None

    @model_validator(mode="after")
    def check_values_for_conditional(self):
        if (self.conditional_var != 1) and (self.conditional_var != 2):
            raise ValueError(f"conditional_var should be 1 or 2 (currently {self.conditional_var}).")
        
        if (self.value_for_conditional < 0) or (self.value_for_conditional > 1):
            raise ValueError(f"value_for_conditional should be between 0 and 1 (currently {self.value_for_conditional}).")

        return self

    @model_validator(mode="after")
    def check_var_types(self):
        if not all(x == "c" for x in self.bicop.var_types):
            raise ValueError("Copula contains discrete variables - this is not supported.")
        
        return self


# ======================================================
# Supporting functions
# ======================================================
def _get_fam_par2(paircop):
    """Extracts and translates pair copula family, rotation, and parameters."""
    fam = paircop.family.name
    rot = paircop.rotation
    params = paircop.parameters

    r = 0
    z = 1
    f = 0
    p = 0
    p2 = 0

    if rot == 180:
      r= 1
      z = 1
    elif rot == 90:
      r = 3
      z = -1
    elif rot == 270:
      r = 2
      z = -1
    else:
      r = 0
      z = 1

    if fam == "indep":
        f = 0
        p = 0
        p2 = 0
    elif fam == "gaussian":
        f = 1
        p = params
        p2 = 0
    elif (fam == "student") or (fam == "t"):
        f = 2
        p = params[0]
        p2 = params[1]
    elif fam == "clayton":
        f = 3 + r*10
        p = params*z
        p2 = 0
    elif fam == "gumbel":
        f = 4 + r*10
        p = params*z
        p2 = 0
    elif fam == "frank":
        f = 5
        p = params*z
        p2 = 0
    elif fam == "joe":
        f = 6 + r*10
        p = params*z
        p2 = 0
    elif fam == "bb1":
        f = 7 + r*10
        p = params[0]*z
        p2 = params[1]*z
    elif fam == "bb6":
        f = 8 + r*10
        p = params[0]*z
        p2 = params[1]*z
    elif fam == "bb7":
        f = 9 + r*10
        p = params[0]*z
        p2 = params[1]*z
    elif fam == "bb8":
        f = 10 + r*10
        p = params[0]*z
        p2 = params[1]*z
    else:
        raise NotImplementedError(f"Copula family {fam!r} is not supported.")

    if isinstance(p,np.ndarray):
        p = p[0]
        if isinstance(p,np.ndarray):
            p = p[0]
    
    if isinstance(p2,np.ndarray):
        p2 = p2[0]
        if isinstance(p2,np.ndarray):
            p2 = p2[0]

    return f, p, p2



def _transformation_to_RVM_all(copula):
    """Transforms pyvinecopulib to R's VineCopula style."""
    d = copula.structure.dim
    order = copula.structure.order
    array = copula.structure.struct_array

    p = list(range(d - 1, 0, -1)) # check this

    Matrix = np.diag(order)
    Family = np.zeros((d,d))
    Par = np.zeros((d,d))
    Par2 = np.zeros((d,d))

    for i in range(d - 1, 0, -1):
        for j in range(1, i + 1):
            array_idx_i = p[i - 1] - 1
            array_idx_j = j - 1

            ij = array(array_idx_i,array_idx_j)
            Matrix[i, j - 1] = ij
            
            cop = copula.get_pair_copula(array_idx_i,array_idx_j)
            cop_f, cop_p, cop_p2 = _get_fam_par2(cop)

            Family[i, j - 1] = cop_f
            Par[i, j - 1] = cop_p
            Par2[i, j - 1] = cop_p2


    RVM = {
        "Matrix": Matrix,
        "Family": Family,
        "Par": Par,
        "Par2": Par2
    }

    return RVM



def _reorderRVineMatrix(Matrix, old_order=None):
    """Re-order RVine matrix in R's VineCopula style."""
    matrix_new = np.array(Matrix, copy=True)
    d = Matrix.shape[0]

    if old_order is None:
        old_order = np.diag(Matrix)

    for i in range(d):
        matrix_new[Matrix == old_order[i]] = d - i

    return matrix_new



def _normalizeRVineMatrix(Matrix):
    return _reorderRVineMatrix(Matrix)



def _createMaxMat(Matrix):
    """Create helper MaxMatr (translated exactly from R VineCopula library)"""

    if Matrix.shape[0] != Matrix.shape[1]:
        raise ValueError("Structure matrix has to be quadratic.")

    MaxMat = _reorderRVineMatrix(Matrix) # Can be simplified to MaxMat = np.array(Matrix, copy=True) because we'll be passing a re-ordered matrix

    n = Matrix.shape[0]

    for j in range(1,n):
        for i in range(n-1,j-1,-1):
            j_idx = j - 1
            i_idx = i - 1

            MaxMat[i_idx,j_idx] = np.max(MaxMat[i_idx:(i_idx+1+1),j_idx])

    tMaxMat = np.array(MaxMat, copy=True)
    tMaxMat[np.isnan(tMaxMat)] = 0

    oldSort = np.diag(Matrix)
    oldSort = oldSort[::-1]

    for i in range(1,n+1):
        i_idx = i - 1

        MaxMat[tMaxMat == i] = oldSort[i_idx]

    return MaxMat



def _neededCondDistr(Vine):
    """Create helper matrices (translated exactly from R VineCopula library)"""

    if Vine.shape[0] != Vine.shape[1]:
        raise ValueError("Structure matrix has to be quadratic.")

    Vine = _reorderRVineMatrix(Vine) # Unnecessary because we're already passing a re-ordered matrix

    MaxMat = _createMaxMat(Vine)

    d = MaxMat.shape[0]
    if d <= 2:
        raise ValueError(f"Dimension is {d} - did you pass a vine copula?")

    M = dict()

    M["direct"] = np.full((d,d), False)
    M["indirect"] = np.full((d,d), False)
    
    M["direct"][1:d,0] = True

    for i in range(2,d):
        v = d - i + 1

        i_idx = i - 1

        bw = MaxMat[i_idx:d,0:i_idx] == v

        direct = Vine[i_idx:d,0:i_idx] == v

        M["indirect"][i_idx:d,i_idx] = np.any(bw & ~direct,axis=1)
        
        M["direct"][i_idx:d,i_idx] = True

        M["direct"][i_idx,i_idx] = np.any(bw[0,:] & direct[0,:])

    return M



# ======================================================
# User-facing functions
# ======================================================
def bicop_conditional_sample(
        bicop: pv.Bicop,
        value_for_conditional: float,
        conditional_var: int,
        num_samples: int = 1000,
        seed: int | None = None
    ) -> np.ndarray:
    """Conditionally sample a bivariate copula.
    
    You pass the fit copula object, the index (1 or 2) of the conditional variable and its value.
    The function will sample the other variable conditional on this using the inverse Rosenblatt.
    
    Parameters
    ----------
    bicop : pv.Bicop
        Bivariate copula object to conditionally sample.
    value_for_conditional : float
        Value of the conditional variable (between 0 and 1).
    conditional_var : int
        Index of the conditional variable (1 or 2).
    num_samples : int, optional
        Number of samples to draw
    seed : int, optional
        Seed for reproducibility.

    Returns
    -------
    np.ndarray
        Conditional samples.

    Examples
    --------
    >>> import pyvinecopulib as pv
    >>> import numpy as np
    
    Specify the bivariate copula:
    >>> bicop = pv.Bicop.from_family(
        family = pv.gaussian,
        rotation = 0,
        parameters = np.array([[0.9]])
    )

    Sample it:
    >>> sample = bicop_conditional_sample(
        bicop,
        value_for_conditional = 0.5,
        conditional_var = 2,
        seed = 12345
    )
    >>> print(sample[:10])
    [0.37225283 0.41768285 0.64161105 0.57899029 0.45204658 0.42529231
    0.54321012 0.34903022 0.57733348 0.75313433]
    
    """
    # Validation
    _ = _bicop_validation(**locals())

    # Sample from uniform and prepare array
    if seed is None:
        uniform_samples = np.random.uniform(0,1,num_samples)
    else:
        rng = np.random.default_rng(seed=seed)
        uniform_samples = rng.uniform(0, 1, num_samples)
    
    # Clip for numerical safety
    eps = 1e-12
    uniform_samples = np.clip(uniform_samples, eps, 1.0 - eps)

    # Prepare arrays
    uniform_samples = np.asfortranarray(uniform_samples.reshape(-1, 1).astype(np.float64))
    value_for_conditional = float(value_for_conditional)
    value_for_conditional = min(max(value_for_conditional, eps), 1.0 - eps)
    value_for_conditional_array = np.full_like(uniform_samples,value_for_conditional)

    # Compute inverse Rosenblatt
    if conditional_var == 1:
        u = np.hstack([value_for_conditional_array,uniform_samples])
        u_other = bicop.hinv1(u=u)
    else:
        u = np.hstack([uniform_samples,value_for_conditional_array])
        u_other = bicop.hinv2(u=u)

    return u_other # type: ignore



def vinecop_conditional_sample(
        vinecop: pv.Vinecop,
        values_for_conditional: dict[int, float],
        num_samples: int = 1000,
        burnin: int = 1000,
        thin: int = 10,
        **kwargs
    ) -> np.ndarray:
    """Conditionally sample an arbitrary vine copula.
    
    The indices and values of the conditional variables are passed via values_for_conditional.
    The function will then sample the remaining variables in the copula, conditional on those values.

    Sampling is done using Hamiltonian Monte Carlo.

    This is a wrapper function, which calls Stan behind the scenes, and requires cmdstan to be
    installed (see Notes).

    You may pass "seed" via additional keyword arguments.
    
    Parameters
    ----------
    vinecop : pv.Vinecop
        Vine copula object to conditionally sample.
    values_for_conditional : dict
        Dictionary of values for the conditional variables.
        Keys should be indices of the conditional variables (1-index based).
    num_samples : int, optional
        Number of samples to draw
    burnin : int, optional
        Number of samples in the initial phase of MCMC simulation to discard.
        Default is 1,000.
    thin : int, optional
        Parameter governing thinning in MCMC. Default is 10.

    Returns
    -------
    np.ndarray
        Sample of the conditioned variables in order, e.g., if we have a copula of five
        variables and we specify 1, 3, and 5 as conditional variables, the first column
        in the output will be the values for variable 2 and the second column will
        be values of variable 4.
    
    Notes
    ----------
    Cmdstan may be installed using:
    >>> install_cmdstan()

    On Windows, you may need to set the compiler argument of the install_cmdstan function
    equal to True.
    
    References
    ----------
    .. [1] Hanebeck, A., Sahin, O., Havlickova, P. et al. Sampling from Conditional Distributions
    of Simplified Vines. Stat Comput 35, 128 (2025). https://doi.org/10.1007/s11222-025-10652-4
    
    Examples
    --------
    >>> import pyvinecopulib as pv
    >>> import numpy as np
    
    Specify the vine copula:

    >>> bicop11 = pv.Bicop.from_family(
        family = pv.bb1,
        rotation = 0,
        parameters = np.array([[1.0], [2.0]])
    )
    >>> bicop12 = pv.Bicop.from_family(
        family = pv.gumbel,
        rotation = 0,
        parameters = np.array([[3.0]])
    )
    >>> bicop13 = pv.Bicop.from_family(
        family = pv.gaussian,
        rotation = 0,
        parameters = np.array([[0.9]])
    )
    >>> bicop14 = pv.Bicop.from_family(
        family = pv.clayton,
        rotation = 0,
        parameters = np.array([[3.0]])
    )
    >>> bicop21 = pv.Bicop.from_family(
        family = pv.clayton,
        rotation = 0,
        parameters = np.array([[1.4]])
    )
    >>> bicop22 = pv.Bicop.from_family(
        family = pv.frank,
        rotation = 0,
        parameters = np.array([[2.8]])
    )
    >>> bicop23 = pv.Bicop.from_family(
        family = pv.bb7,
        rotation = 0,
        parameters = np.array([[1.0], [1.5]])
    )
    >>> bicop31 = pv.Bicop.from_family(
        family = pv.frank,
        rotation = 0,
        parameters = np.array([[1.8]])
    )
    >>> bicop32 = pv.Bicop.from_family(
        family = pv.gaussian,
        rotation = 0,
        parameters = np.array([[0.3]])
    )
    >>> bicop41 = pv.Bicop.from_family(
        family = pv.gaussian,
        rotation = 0,
        parameters = np.array([[0.05]])
    )
    >>> structure_matrix = np.array([
        [2,3,4,5,5],
        [3,4,5,4,0],
        [4,5,3,0,0],
        [5,2,0,0,0],
        [1,0,0,0,0]
    ])
    >>> vinecop = pv.Vinecop.from_structure(
        matrix = structure_matrix,
        pair_copulas = [
            [bicop11,bicop12,bicop13,bicop14],
            [bicop21,bicop22,bicop23],
            [bicop31,bicop32],
            [bicop41]
        ]
    )

    Specify conditional variables and their values:

    >>> values_for_conditional = {
        1: 0.2,
        3: 0.45,
        5: 0.78
    }

    Sample the copula:
    >>> sample = vinecop_conditional_sample(
            vinecop = vinecop,
            values_for_conditional = values_for_conditional,
            seed = 12345
        )
    >>> print(sample)
    [[0.3514476  0.36845809]
    [0.2122368  0.36747993]
    [0.96914275 0.47400057]
    ...
    [0.32418243 0.28573965]
    [0.2134473  0.35859307]
    [0.41961971 0.1817589 ]]
    
    """
    # Validate inputs
    _ = _vinecop_validation(
            vinecop=vinecop,
            values_for_conditional=values_for_conditional, # type: ignore
            num_samples=num_samples,
            burnin=burnin,
            thin=thin
        )
    seed = kwargs.get("seed", None)
    if seed is not None:
        if not isinstance(seed,(int,np.integer)):
            raise ValueError("seed should be an integer.")
        if isinstance(seed,np.integer):
            seed = int(seed)

    # Extract dimensions
    dim = vinecop.dim
    num_conditional_vars = len(values_for_conditional)
    num_conditioned_vars = dim - num_conditional_vars

    # Create a lists of indices of conditional and conditioned variables
    # and a list of values for conditional variables
    conditional_vars = list(values_for_conditional.keys())

    conditional_vars_idx = [i if i in conditional_vars else False for i in range(1,dim+1)]
    conditional_values = [values_for_conditional[i] if i in conditional_vars else False for i in range(1,dim+1)]

    # Clip conditional values for numerical stability
    eps = 1e-12
    conditional_values = [float(x) if x else False for x in conditional_values]
    conditional_values = [min(max(x,eps),1 - eps) if x else False for x in conditional_values]

    # Load Stan
    if num_conditioned_vars == 1:
        stan_to_use = stan1_model
    else:
        stan_to_use = stan2_model
    
    # Create copula representation in R's VineCopula style
    RVM = _transformation_to_RVM_all(vinecop)

    # Normalize RVineMatrix - needed in R VineCopula library to compute log-lik
    dataflip = 0
    vine_matrix_diag = np.diag(RVM["Matrix"])

    dim_to_1 = np.array(range(dim,0,-1))
    if (vine_matrix_diag != dim_to_1).any():
        RVM["Matrix"] = _normalizeRVineMatrix(RVM["Matrix"])
        dataflip = 1

    # From the VineCopula package: the different parameters/inputs we need 
    # to compute log-lik in Stan - inserted into Stan as data
    T = 1
    families = RVM["Family"].ravel(order="F")
    families[np.isnan(families)] = 0

    parameters = RVM["Par"].ravel(order="F")
    parameters[np.isnan(parameters)] = 0

    parameters2 = RVM["Par2"].ravel(order="F")
    parameters2[np.isnan(parameters2)] = 0

    # Prepare other inputs
    M = _neededCondDistr(RVM["Matrix"])
    
    condirect = (M["direct"]*1).ravel("F")
    condirect[np.isnan(condirect)] = 0

    conindirect = (M["indirect"]*1).ravel("F")
    conindirect[np.isnan(conindirect)] = 0

    maxmat = _createMaxMat(RVM["Matrix"]).ravel("F")
    maxmat[np.isnan(maxmat)] = 0
    
    matri = RVM["Matrix"].ravel("F")
    matri[np.isnan(matri)] = 0

    # Parameters for the Stan-program
    num_samples_original = num_samples
    num_samples = int(np.ceil(num_samples/_NUM_CHAINS)) # per-chain draws to reach total >= original

    It = num_samples*thin
    iter_u_2 = It+burnin
    burnin_u_2 = burnin
    chains_u_2 = _NUM_CHAINS
    adapt_delta_u_2 = _ADAPT_DELTA
    max_treedepth_u_2 = _MAX_TREEDEPTH

    # Data for Stan
    data_stan_u_2 = {
        "T": T,
        "dataflip": dataflip,
        "d": dim,
        "o": vine_matrix_diag,
        "d1": num_conditioned_vars,
        "d2": num_conditional_vars,
        "family": families.astype(int).tolist(),
        "maxmat": maxmat.tolist(),
        "matri": matri.tolist(),
        "condirect": condirect.tolist(),
        "conindirect": conindirect.tolist(),
        "par": parameters.astype(float).tolist(),
        "par2": parameters2.astype(float).tolist(),
        "indexcon": conditional_vars_idx,
        "ucon": conditional_values
    }

    repeats = np.repeat(0.5,num_conditioned_vars)
    init_list_u_2 = {
        1: {"ucalculate": repeats},
        2: {"ucalculate": repeats},
        3: {"ucalculate": repeats},
        4: {"ucalculate": repeats}
    }

    # Conditional sampling in Stan
    results = stan_to_use.sample(
        iter_sampling = iter_u_2-burnin_u_2,
        iter_warmup = burnin_u_2,
        chains = chains_u_2,
        data = data_stan_u_2,
        inits = init_list_u_2,
        adapt_delta = adapt_delta_u_2,
        max_treedepth = max_treedepth_u_2,
        thin = thin,
        show_progress=False,
        show_console=False,
        **kwargs
    )

    # Permutations
    seed = kwargs.get("seed", None)
    rng = np.random.default_rng(seed=seed)
    permuted_draws = rng.permutation(results.stan_variables()["ucalculate"])
    permuted_draws = permuted_draws[:num_samples_original]

    return permuted_draws



def copula_conditional_sample(
        copula: pv.Bicop | pv.Vinecop,
        values_for_conditional: dict[int, float],
        num_samples: int = 1000, # type: ignore
        **kwargs
    ) -> np.ndarray:
    """Conditionally sample an arbitrary copula.

    Wrapper around functions bicop_conditional_sample and vinecop_conditional_sample, allowing you to
    conditionally sample either a bivariate or a vine copula.
    
    The index (or indices) and value(s) of the conditional variables are passed via values_for_conditional.
    The function will then sample the remaining variable(s) in the copula, conditional on those value(s).

    You may pass "seed" via additional keyword arguments. See documentation for vinecop_conditional_sample
    and bicop_conditional_sample for what other arguments may be passed.
    
    Parameters
    ----------
    copula : pv.Bicop or pv.Vinecop
        Bivariate or vine copula object to conditionally sample.
    values_for_conditional : dict
        Dictionary of values for the conditional variables. If copula is bivariate, there should be only one entry.
        Keys should be indices of the conditional variables (1-index based; 1 or 2 for bivariate copula).
    num_samples : int, optional
        Number of samples to draw

    Returns
    -------
    np.ndarray
        Sample of the conditioned variables in order, e.g., if we have a copula of five
        variables and we specify 1, 3, and 5 as conditional variables, the first column
        in the output will be the values for variable 2 and the second column will
        be values of variable 4.
    """
    if not isinstance(values_for_conditional,dict):
        raise ValueError("values_for_conditional should be a dictionary.")

    if isinstance(copula,pv.Bicop):
        if len(values_for_conditional) != 1:
            raise ValueError("You passed bivariate copula: values_for_conditional dictionary should have exactly one entry.")
        
        (conditional_var, value_for_conditional), = values_for_conditional.items()

        result = bicop_conditional_sample(
            bicop = copula,
            value_for_conditional = value_for_conditional,
            conditional_var = conditional_var,
            num_samples = num_samples,
            **kwargs
        )
    elif isinstance(copula,pv.Vinecop):
        result = vinecop_conditional_sample(
            vinecop = copula,
            values_for_conditional = values_for_conditional,
            num_samples = num_samples,
            **kwargs
        )
    else:
        raise TypeError("copula must be a pv.Bicop or pv.Vinecop instance.")
    
    return result




# ======================================================
# Deprecated
# ======================================================
def _sample_from_conditional_orig(
        N,
        RVM,
        indexcon,
        ucon,
        burnin = 1000,
        thin = 10,
        **kwargs
    ):
    """Direct translation of original function"""
    # Extracting the dimensions
    d = len(indexcon)
    d1 = sum([1 for x in indexcon if not x])
    d2 = d - d1

    # Load Stan
    if d1 == 1:
        stan_to_use = stan1_model
    else:
        stan_to_use = stan2_model
    
    # Normalizing RVineMatrix - needed in VineCopula package to compute log-lik
    dataflip = 0
    o = np.diag(RVM["Matrix"])

    len_o_to_1 = np.array(range(len(o),0,-1))
    if (o != len_o_to_1).any():
        # oldRVM = dc(RVM)
        RVM["Matrix"] = _normalizeRVineMatrix(RVM["Matrix"])
        dataflip = 1
    
    # From the VineCopula package: the different parameters/inputs we need 
    # to compute log-lik in Stan - inserted into Stan as data
    T = 1
    w1 = RVM["Family"].ravel(order="F")
    w1[np.isnan(w1)] = 0

    th = RVM["Par"].ravel(order="F")
    th[np.isnan(th)] = 0

    th2 = RVM["Par2"].ravel(order="F")
    th2[np.isnan(th2)] = 0

    M = _neededCondDistr(RVM["Matrix"])
    condirect = (M["direct"]*1).ravel("F")
    conindirect = (M["indirect"]*1).ravel("F")
    maxmat = _createMaxMat(RVM["Matrix"]).ravel("F")
    matri = RVM["Matrix"].ravel("F")

    matri[np.isnan(matri)] = 0
    maxmat[np.isnan(maxmat)] = 0
    condirect[np.isnan(condirect)] = 0
    conindirect[np.isnan(conindirect)] = 0

    # The parameters for the Stan-program
    It = N*thin
    iter_u_2 = It+burnin
    burnin_u_2 = burnin
    chains_u_2 = 4
    adapt_delta_u_2 = 0.8
    max_treedepth_u_2 = 10

    # Sample from Stan
    data_stan_u_2 = {
        "T": T,
        "dataflip": dataflip,
        "d": d,
        "o": o.tolist(), # Is this right?
        "d1": d1,
        "d2": d2,
        "family": w1.astype(int).tolist(),
        "maxmat": maxmat.tolist(),
        "matri": matri.tolist(),
        "condirect": condirect.tolist(),
        "conindirect": conindirect.tolist(),
        "par": th.astype(float).tolist(),
        "par2": th2.astype(float).tolist(),
        "indexcon": indexcon,
        "ucon": ucon
    }
    
    repeats = np.repeat(0.5,d1)
    init_list_u_2 = {
        1: {"ucalculate": repeats},
        2: {"ucalculate": repeats},
        3: {"ucalculate": repeats},
        4: {"ucalculate": repeats}
    }

    results = stan_to_use.sample(
        iter_sampling = iter_u_2-burnin_u_2,
        iter_warmup = burnin_u_2,
        chains = chains_u_2,
        data = data_stan_u_2,
        inits = init_list_u_2,
        adapt_delta = adapt_delta_u_2,
        max_treedepth = max_treedepth_u_2,
        thin = thin,
        # seed = seed,
        show_progress=False,
        show_console=False,
        **kwargs
    )

    permuted_draws = np.random.permutation(results.stan_variables()["ucalculate"])

    return permuted_draws




