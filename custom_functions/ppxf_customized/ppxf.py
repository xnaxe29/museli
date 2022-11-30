"""
    Copyright (C) 2001-2020, Michele Cappellari

    E-mail: michele.cappellari_at_physics.ox.ac.uk

    Updated versions of the software are available from my web page
    http://purl.org/cappellari/software

    If you have found this software useful for your research,
    I would appreciate an acknowledgement to the use of the
    "Penalized Pixel-Fitting method by Cappellari & Emsellem (2004)
    as upgraded in Cappellari (2017)".

    http://ui.adsabs.harvard.edu/abs/2004PASP..116..138C
    http://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C

    This software is provided as is without any warranty whatsoever.
    Permission to use, for non-commercial purposes is granted.
    Permission to modify for personal or internal use is granted,
    provided this copyright and disclaimer are included unchanged
    at the beginning of the file. All other rights are reserved.
    In particular, redistribution of the code is not allowed.

"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import legendre, hermite
from scipy import optimize, linalg, special

try:
    import cvxopt
    cvxopt_installed = True
except ImportError:
    cvxopt_installed = False

from .capfit import capfit, lsq_lin, lsq_box, cov_err, lsq_lin_cvxopt

###############################################################################

def trigvander(x, deg):
    """
    Analogue to legendre.legvander(), but for a trigonometric
    series rather than Legendre polynomials:

    `deg` must be an even integer

    """
    assert deg % 2 == 0, "`degree` must be even with trig=True"

    u = np.pi*x[:, None]   # [-pi, pi] interval
    j = np.arange(1, deg//2 + 1)
    mat = np.ones((x.size, deg + 1))
    mat[:, 1:] = np.hstack([np.cos(j*u), np.sin(j*u)])

    return mat

################################################################################

def trigval(x, c):
    """
    Analogue to legendre.legval(), but for a trigonometric
    series rather than Legendre polynomials:

    Evaluate a trigonometric series with coefficients `c` at points `x`.

    """
    return trigvander(x, c.size - 1) @ c

################################################################################

def rebin(x, factor):
    """
    Rebin a vector, or the first dimension of an array,
    by averaging within groups of "factor" adjacent values.
    The dimensionality of the input is retained.
    In particular, in the 1-dim case, a row-vector
    or colum-vector are retained as such in output.

    """
    if factor > 1:
        n = x.ndim
        x = x.reshape(len(x)//factor, factor, -1).mean(1)
        if n == 1:
            x = x.squeeze()

    return x

################################################################################

def robust_sigma(y, zero=False):
    """
    Biweight estimate of the scale (standard deviation).
    Implements the approach described in
    "Understanding Robust and Exploratory Data Analysis"
    Hoaglin, Mosteller, Tukey ed., 1983, Chapter 12B, pg. 417

    """
    y = np.ravel(y)
    d = y if zero else y - np.median(y)

    mad = np.median(np.abs(d))
    u2 = (d/(9.0*mad))**2  # c = 9
    good = u2 < 1.0
    u1 = 1.0 - u2[good]
    num = y.size * ((d[good]*u1**2)**2).sum()
    den = (u1*(1.0 - 5.0*u2[good])).sum()
    sigma = np.sqrt(num/(den*(den - 1.0)))  # see note in above reference

    return sigma

################################################################################

def reddening_cal00(lam, ebv):
    """
    Reddening curve of `Calzetti et al. (2000)
    <http://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C>`_
    This is reliable between 0.12 and 2.2 micrometres.
    - LAMBDA is the restframe wavelength in Angstrom of each pixel in the
      input galaxy spectrum (1 Angstrom = 1e-4 micrometres)
    - EBV is the assumed E(B-V) colour excess to redden the spectrum.
      In output the vector FRAC gives the fraction by which the flux at each
      wavelength has to be multiplied, to model the dust reddening effect.

    """
    ilam = 1e4/lam  # Convert Angstrom to micrometres and take 1/lambda
    rv = 4.05  # C+00 equation (5)

    # C+00 equation (3) but extrapolate for lam > 2.2
    # C+00 equation (4) (into Horner form) but extrapolate for lam < 0.12
    k1 = rv + np.where(lam >= 6300, 2.76536*ilam - 4.93776,
                       ilam*((0.029249*ilam - 0.526482)*ilam + 4.01243) - 5.7328)
    fact = 10**(-0.4*ebv*k1.clip(0))  # Calzetti+00 equation (2) with opposite sign

    return fact # The model spectrum has to be multiplied by this vector

################################################################################

def losvd_rfft(pars, nspec, moments, nl, ncomp, vsyst, factor, sigma_diff):
    """
    Analytic Fourier Transform (of real input) of the Gauss-Hermite LOSVD.
    Equation (38) of `Cappellari (2017)
    <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C>`_

    """
    losvd_rfft = np.empty((nl, ncomp, nspec), dtype=complex)
    p = 0
    for j, mom in enumerate(moments):  # loop over kinematic components
        for k in range(nspec):  # nspec=2 for two-sided fitting, otherwise nspec=1
            s = 1 if k == 0 else -1  # s=+1 for left spectrum, s=-1 for right one
            #print (pars[1 + p])
            vel, sig = vsyst + s*pars[0 + p], pars[1 + p]
            #print (sig)
            a, b = [vel, sigma_diff]/sig
            w = np.linspace(0, np.pi*factor*sig, nl)
            losvd_rfft[:, j, k] = np.exp(1j*a*w - 0.5*(1 + b**2)*w**2)

            if mom > 2:
                n = np.arange(3, mom + 1)
                nrm = np.sqrt(special.factorial(n)*2**n)   # vdMF93 Normalization
                coeff = np.append([1, 0, 0], (s*1j)**n * pars[p - 1 + n]/nrm)
                poly = hermite.hermval(w, coeff)
                losvd_rfft[:, j, k] *= poly
        p += mom

    return np.conj(losvd_rfft)

################################################################################

def regularization(a, npoly, npix, nspec, reg_dim, reg_ord, regul):
    """
    Add first or second order 1D, 2D or 3D linear regularization.
    Equation (25) of `Cappellari (2017)
    <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C>`_

    """
    #print (npoly, a.shape[1])
    #print (npoly, len(a[:, npoly : npoly + np.prod(reg_dim)]))
    #print (np.prod(reg_dim))
    b = a[:, npoly : npoly + np.prod(reg_dim)].reshape(-1, *reg_dim)
    p = npix*nspec

    if reg_ord == 1:   # Minimize integral of (Grad[w] @ Grad[w])
        diff = np.array([1, -1])*regul
        if reg_dim.size == 1:
            for j in range(reg_dim[0] - 1):
                b[p, j : j + 2] = diff
                p += 1
        elif reg_dim.size == 2:
            for k in range(reg_dim[1]):
                for j in range(reg_dim[0]):
                    if j < reg_dim[0] - 1:
                        b[p, j : j + 2, k] = diff
                        p += 1
                    if k < reg_dim[1] - 1:
                        b[p, j, k : k + 2] = diff
                        p += 1
        elif reg_dim.size == 3:
            for q in range(reg_dim[2]):
                for k in range(reg_dim[1]):
                    for j in range(reg_dim[0]):
                        if j < reg_dim[0] - 1:
                            b[p, j : j + 2, k, q] = diff
                            p += 1
                        if k < reg_dim[1] - 1:
                            b[p, j, k : k + 2, q] = diff
                            p += 1
                        if q < reg_dim[2] - 1:
                            b[p, j, k, q : q + 2] = diff
                            p += 1
    elif reg_ord == 2:   # Minimize integral of Laplacian[w]**2
        diff = np.array([1, -2, 1])*regul
        if reg_dim.size == 1:
            for j in range(1, reg_dim[0] - 1):
                b[p, j - 1 : j + 2] = diff
                p += 1
        elif reg_dim.size == 2:
            for k in range(reg_dim[1]):
                for j in range(reg_dim[0]):
                    if 0 < j < reg_dim[0] - 1:
                        b[p, j - 1 : j + 2, k] = diff
                    if 0 < k < reg_dim[1] - 1:
                        b[p, j, k - 1 : k + 2] += diff
                    p += 1
        elif reg_dim.size == 3:
            for q in range(reg_dim[2]):
                for k in range(reg_dim[1]):
                    for j in range(reg_dim[0]):
                        if 0 < j < reg_dim[0] - 1:
                            b[p, j - 1 : j + 2, k, q] = diff
                        if 0 < k < reg_dim[1] - 1:
                            b[p, j, k - 1 : k + 2, q] += diff
                        if 0 < q < reg_dim[2] - 1:
                            b[p, j, k, q - 1 : q + 2] += diff
                        p += 1

################################################################################

class ppxf:
    """
    pPXF Purpose
    ------------

    Extract galaxy stellar kinematics ``(V, sigma, h3, h4, h5, h6,...)``
    or the stellar population and gas emission by fitting a template
    to an observed spectrum in pixel space, using the Penalized PiXel-Fitting
    (``pPXF``) method originally described in

    `Cappellari & Emsellem (2004) <https://ui.adsabs.harvard.edu/abs/2004PASP..116..138C>`_

    and substantially upgraded in subsequent years and particularly in

    `Cappellari (2017) <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C>`_

    The following key optional features are also available:

    1)  An optimal template, positive linear combination of different input
        templates, can be fitted together with the kinematics.
    2)  One can enforce smoothness on the template weights during the fit. This
        is useful to attach a physical meaning to the weights e.g. in terms of
        the star formation history of a galaxy.
    3)  One can fit multiple kinematic components for both the stars and the gas
        emission lines. Both the stellar and gas LOSVD can be penalized and can
        be described by a general Gauss-Hermite series.
    4)  Any parameter of the LOSVD (e.g. sigma) for any kinematic component can
        either be fitted or held fixed to a given value, while other parameters
        are fitted. Alternatively, parameters can be constrained to lie within
        given limits or even tied by simple relations to other parameters.
    5)  One can enforce linear equality/inequality constraints on either the
        template weights or the kinematic parameters.
    6)  Additive and/or multiplicative polynomials can be included to adjust the
        continuum shape of the template to the observed spectrum.
    7)  Iterative sigma clipping can be used to clean the spectrum.
    8)  It is possible to fit a mirror-symmetric LOSVD to two spectra at the
        same time. This is useful for spectra taken at point-symmetric spatial
        positions with respect to the center of an equilibrium stellar system.
    9)  One can include sky spectra in the fit, to deal with cases where the sky
        dominates the observed spectrum and an accurate sky subtraction is
        critical.
    10) One can derive an estimate of the reddening in the spectrum. This can be
        done independently for the stellar spectrum or the Balmer emission lines.
    11) The covariance matrix can be input instead of the error spectrum, to
        account for correlated errors in the spectral pixels.
    12) One can specify the weights fraction between two kinematics components,
        e.g. to model bulge and disk contributions.
    13) One can use templates with higher resolution than the galaxy, to
        improve the accuracy of the LOSVD extraction at low dispersion.


    Calling Sequence
    ----------------

    .. code-block:: python

        from ppxf.ppxf import ppxf

        pp = ppxf(self, templates, galaxy, noise, velscale, start,
                 bias=None, bounds=None, clean=False, component=0,
                 constr_templ={}, constr_kinem={}, degree=4, fixed=None,
                 fraction=None, ftol=1e-4, gas_component=None, gas_names=None,
                 gas_reddening=None, goodpixels=None, lam=None, linear=False,
                 linear_method='lsq_box', mask=None, method='capfit', mdegree=0,
                 moments=2, plot=False, quiet=False, reddening=None,
                 reddening_func=None, reg_dim=None, reg_ord=2, regul=0,
                 sigma_diff=0, sky=None, templates_rfft=None, tied=None,
                 trig=False, velscale_ratio=1, vsyst=0, x0=None)

        print(pp.sol)  # print best-fitting kinematics (V, sigma, h3, h4)
        pp.plot()      # Plot best fit and gas lines

    Input Parameters
    ----------------

    templates:
        Vector containing the spectrum of a single template star or more
        commonly an array of dimensions ``templates[nPixels, nTemplates]``
        containing different templates to be optimized during the fit of the
        kinematics. It has to be ``nPixels >= galaxy.size``.

        To apply linear regularization to the ``weights`` via the keyword
        ``regul``, ``templates`` should be an array of two
        ``templates[nPixels, nAge]``, three ``templates[nPixels, nAge, nMetal]``
        or four ``templates[nPixels, nAge, nMetal, nAlpha]`` dimensions,
        depending on the number of population variables one wants to study.
        This can be useful to try to attach a physical meaning to the output
        ``weights``, in term of the galaxy star formation history and chemical
        composition distribution.
        In that case the templates may represent single stellar population SSP
        models and should be arranged in sequence of increasing age,
        metallicity or alpha along the second, third or fourth dimension of the
        array respectively.
    galaxy:
        Vector containing the spectrum of the galaxy to be measured. The
        star and the galaxy spectra have to be logarithmically rebinned but the
        continuum should *not* be subtracted. The rebinning may be performed
        with the ``log_rebin`` routine in ``ppxf.ppxf_util``.

        For high redshift galaxies, one should bring the spectra close to the
        restframe wavelength, before doing the ``pPXF`` fit. This can be done
        by dividing the observed wavelength by ``(1 + z)``, where ``z`` is a
        rough estimate of the galaxy redshift, before the logarithmic
        rebinning. See Section 2.4 of `Cappellari (2017)`_ for details.

        ``galaxy`` can also be an array of dimensions ``galaxy[nGalPixels, 2]``
        containing two spectra to be fitted, at the same time, with a
        reflection-symmetric LOSVD. This is useful for spectra taken at
        point-symmetric spatial positions with respect to the center of an
        equilibrium stellar system.
        For a discussion of the usefulness of this two-sided fitting see e.g.
        Section 3.6 of `Rix & White (1992)
        <http://ui.adsabs.harvard.edu/abs/1992MNRAS.254..389R>`_.

        IMPORTANT: (1) For the two-sided fitting the ``vsyst`` keyword has to
        be used. (2) Make sure the spectra are rescaled to be not too many
        order of magnitude different from unity, to avoid over or underflow
        problems in the calculation. E.g. units of ``erg/(s cm^2 A)`` may cause
        problems!
    noise:
        Vector containing the ``1*sigma`` uncertainty (per pixel) in the galaxy
        spectrum, or covariance matrix describing the correlated uncertainties
        in the galaxy spectrum. Of course this vector/matrix must have the same
        units as the galaxy spectrum.

        If ``galaxy`` is a ``Nx2`` array, ``noise`` has to be an array with the
        same dimensions.

        When ``noise`` has dimensions ``NxN`` it is assumed to contain the
        covariance matrix with elements ``sigma(i, j)``. When the errors in the
        spectrum are uncorrelated it is mathematically equivalent to input in
        ``pPXF`` an error vector ``noise=errvec`` or a ``NxN`` diagonal matrix
        ``noise=np.diag(errvec**2)`` (note squared!).

        IMPORTANT: the penalty term of the ``pPXF`` method is based on the
        *relative* change of the fit residuals. For this reason, the penalty
        will work as expected even if no reliable estimate of the ``noise`` is
        available (see `Cappellari & Emsellem (2004)`_ for details). If no
        reliable noise is available this keyword can just be set to::

            noise = np.ones_like(galaxy)  # Same weight for all pixels

    velscale:
        Velocity scale of the spectra in km/s per pixel. It has to be the
        same for both the galaxy and the template spectra.
        An exception is when the ``velscale_ratio`` keyword is used, in which
        case one can input ``templates`` with smaller ``velscale`` than
        ``galaxy``.

        ``velscale`` is *defined* in ``pPXF`` by
        ``velscale = c*Delta[np.log(lambda)]``, which is approximately
        ``velscale ~ c*Delta(lambda)/lambda``.
        See Section 2.3 of `Cappellari (2017)`_ for details.
    start:
        Vector, or list/array of vectors ``[start1, start2, ...]``, with the
        initial estimate for the LOSVD parameters.

        When LOSVD parameters are not held fixed, each vector only needs to
        contain ``start = [velStart, sigmaStart]`` the initial guess for the
        velocity and the velocity dispersion in km/s. The starting values for
        h3-h6 (if they are fitted) are all set to zero by default.
        In other words, when ``moments=4``::

            start = [velStart, sigmaStart]

        is interpreted as::

            start = [velStart, sigmaStart, 0, 0]

        When the LOSVD for some kinematic components is held fixed (see
        ``fixed`` keyword), all values for ``[Vel, Sigma, h3, h4,...]`` can be
        provided.

        Unless a good initial guess is available, it is recommended to set the
        starting ``sigma >= 3*velscale`` in km/s (i.e. 3 pixels). In fact, when
        the sigma is very low, and far from the true solution, the ``chi^2`` of
        the fit becomes weakly sensitive to small variations in sigma (see
        ``pPXF`` paper). In some instances, the near-constancy of ``chi^2`` may
        cause premature convergence of the optimization.

        In the case of two-sided fitting a good starting value for the velocity
        is ``velStart = 0.0`` (in this case ``vsyst`` will generally be
        nonzero). Alternatively on should keep in mind that ``velStart`` refers
        to the first input galaxy spectrum, while the second will have velocity
        ``-velStart``.

        With multiple kinematic components ``start`` must be a list of starting
        values, one for each different component.

        EXAMPLE: We want to fit two kinematic components. We fit 4 moments for
        the first component and 2 moments for the second one as follows::

            component = [0, 0, ... 0, 1, 1, ... 1]
            moments = [4, 2]
            start = [[V1, sigma1], [V2, sigma2]]

    Optional Keywords
    -----------------

    bias:
        This parameter biases the ``(h3, h4, ...)`` measurements towards zero
        (Gaussian LOSVD) unless their inclusion significantly decreases the
        error in the fit. Set this to ``bias=0`` not to bias the fit: the
        solution (including ``[V, sigma]``) will be noisier in that case. The
        default ``bias`` should provide acceptable results in most cases, but
        it would be safe to test it with Monte Carlo simulations. This keyword
        precisely corresponds to the parameter ``lambda`` in the
        `Cappellari & Emsellem (2004)`_ paper.
        Note that the penalty depends on the *relative* change of the fit
        residuals, so it is insensitive to proper scaling of the ``noise``
        vector. A nonzero ``bias`` can be safely used even without a reliable
        ``noise`` spectrum, or with equal weighting for all pixels.
    bounds:
        Lower and upper bounds for every kinematic parameter. This is an array,
        or list of arrays, with the same dimensions as ``start``, except for
        the last dimension, which is 2. In practice, for every element of
        ``start`` one needs to specify a pair of values ``[lower, upper]``.

        EXAMPLE: We want to fit two kinematic components, with 4 moments for
        the first component and 2 for the second (e.g. stars and gas). In this
        case::

            moments = [4, 2]
            start_stars = [V1, sigma1, 0, 0]
            start_gas = [V2, sigma2]
            start = [start_stars, start_gas]

        then we can specify boundaries for each kinematic parameter as::

            bounds_stars = [[V1_lo, V1_up], [sigma1_lo, sigma1_up],
                            [-0.3, 0.3], [-0.3, 0.3]]
            bounds_gas = [[V2_lo, V2_up], [sigma2_lo, sigma2_up]]
            bounds = [bounds_stars, bounds_gas]

    component:
        When fitting more than one kinematic component, this keyword should
        contain the component number of each input template. In principle, every
        template can belong to a different kinematic component.

        EXAMPLE: We want to fit the first 50 templates to component 0 and the
        last 10 templates to component 1. In this case::

            component = [0]*50 + [1]*10

        which, in Python syntax, is equivalent to::

            component = [0, 0, ... 0, 1, 1, ... 1]

        This keyword is especially useful when fitting both emissions (gas) and
        absorption (stars) templates simultaneously (see the example for the
        ``moments`` keyword).
    constr_templ:
        It enforces linear constraints on the template weights during the fit.
        This is specified by the following dictionary, where ``A_ineq`` and
        ``A_eq`` are arrays (have ``A.ndim = 2``), while ``b_ineq`` and ``b_eq``
        are vectors (have ``b.ndim = 1``). Either the ``_eq`` or the ``_ineq``
        keys can be omitted if not needed::

            constr_templ = {"A_ineq": A_ineq, "b_ineq": b_ineq, "A_eq": A_eq, "b_eq": b_eq}

        The resulting pPXF solution will satisfy the following linear matrix
        inequalities and/or equalities::

            A_ineq @ pp.weights <= b_ineq
            A_eq @ pp.weights == b_eq

        Inequality can be used e.g. to constrain the fluxes of emission lines to
        lie within prescribed ranges.
        Equalities can be used e.g. to force the weights for different kinematic
        components to contain prescribed fractions of the total weights.

        EXAMPLES: We are fitting a spectrum using four templates, the first two
        templates belong to one kinematic component and the rest to the other.
        (NOTE: This 4-templates example is for illustration, but in real
        applications one will use many more than two templates per component!)
        This implies we have::

            component=[0, 0, 1, 1]

        then we can set the equality constraint that the sum of the weights of
        the first kinematic component is a given ``fraction`` of the total::

            pp.weights[component == 0].sum()/pp.weights.sum() == fraction

        as follows [see equation 30 of `Cappellari (2017)`_]::

            A_eq = [[fraction - 1, fraction - 1, fraction, fraction]]
            b_eq = [0]
            constr_templ = {"A_eq": A_eq, "b_eq": b_eq}

        An identical result can be obtained in this case using the legacy
        ``fraction`` keyword, but ``constr_templ`` additionally allows for
        general linear constraints for multiple kinematic components.

        Similarly, we can set the inequality constraint that the total weights
        of each of the two kinematic components is larger than ``fraction``::

            fraction <= pp.weights[component == 0].sum()/pp.weights.sum()
            fraction <= pp.weights[component == 1].sum()/pp.weights.sum()

        as follows::

            A_ineq = [[fraction - 1, fraction - 1, fraction, fraction],
                      [fraction, fraction, fraction - 1, fraction - 1]]
            b_ineq = [0, 0]
            constr_templ = {"A_ineq": A_ineq, "b_ineq": b_ineq}

        We can constrain the ratio of the first two templates weights to lie in
        the interval ``ratio_min <= w[0]/w[1] <= ratio_max`` as follows::

            A_ineq = [[-1, ratio_min, 0, 0],    # -w[0] + ratio_min*w[1] <= 0
                      [1, -ratio_max, 0, 0]]    # +w[0] - ratio_max*w[1] <= 0
            b_ineq = [0, 0]
            constr_templ = {"A_ineq": A_ineq, "b_ineq": b_ineq}

        If we have six templates for three kinematics components::

            component=[0, 0, 1, 1, 2, 2]

        we can set the fractions for the first two components to be ``fraction1``
        and ``fraction2`` (of the total weights) respectively as follows
        (the third components will be ``1 - fraction1 - fraction2``)::

            A_eq = [[fraction1 - 1, fraction1 - 1, fraction1, fraction1, fraction1, fraction1],
                    [fraction2, fraction2, fraction2 - 1, fraction2 - 1, fraction2, fraction2]]
            b_eq = [0, 0]
            constr_templ = {"A_eq": A_eq, "b_eq": b_eq}

    constr_kinem:
        It enforces linear constraints on the kinematic parameters during the fit.
        This is specified by the following dictionary, where ``A_ineq`` and
        ``A_eq`` are arrays (have ``A.ndim = 2``), while ``b_ineq`` and ``b_eq``
        are vectors (have ``b.ndim = 1``). Either the ``_eq`` or the ``_ineq``
        keys can be omitted if not needed::

            constr_kinem = {"A_ineq": A_ineq, "b_ineq": b_ineq, "A_eq": A_eq, "b_eq": b_eq}

        The resulting pPXF kinematics solution will satisfy the following linear
        matrix inequalities and/or equalities::

            params = np.ravel(pp.sol)  # Unravel for multiple components
            A_ineq @ params <= b_ineq
            A_eq @ params == b_eq

        IMPORTANT: the starting guess ``start`` must satisfy the constraints,
        or in other words, it must lie in the feasible region.

        Inequalities can be used e.g. to force one kinematic component to have
        larger velocity or dispersion than another one. This is useful e.g. when
        extracting two stellar kinematic components or when fitting both narrow
        and broad components of gas emission lines.

        EXAMPLES: We want to fit two kinematic components, with two moments for
        both the first and second component. In this case::

            moments = [2, 2]
            start = [[V1, sigma1], [V2, sigma2]]

        then we can set the constraint ``sigma1 >= 3*sigma2`` as follows::

            A_ineq = [[0, -1, 0, 3]]  # 0*V1 - 1*sigma1 + 0*V2 + 3*sigma2 <= 0
            b_ineq = [0]
            constr_kinem = {"A_ineq": A_ineq, "b_ineq": b_ineq}

        We can set the constraint ``sigma1 >= sigma2 + 2*velscale`` as follows::

            A_ineq = [[0, -1, 0, 1]]  # -sigma1 + sigma2 <= -2*velscale
            b_ineq = [-2]             # kinem. in pixels (-2 --> -2*velscale)!
            constr_kinem =  {"A_ineq": A_ineq, "b_ineq": b_ineq}

        We can set both the constraints ``V1 >= V2`` and
        ``sigma1 >= sigma2 + 2*velscale`` as follows::

            A_ineq = [[-1, 0, 1, 0],   # -V1 + V2 <= 0
                      [0, -1, 0, 1]]   # -sigma1 + sigma2 <= -2*velscale
            b_ineq = [0, -2]           # kinem. in pixels (-2 --> -2*velscale)!
            constr_kinem =  {"A_ineq": A_ineq, "b_ineq": b_ineq}

        We can constrain the velocity dispersion of the second kinematic
        component to differ less than 10% from that of the first component
        ``sigma1/1.1 <= sigma2 <= sigma1*1.1`` as follows::

            A_ineq = [[0, 1/1.1, 0, -1],   # +sigma1/1.1 - sigma2 <= 0
                      [0, -1.1,  0,  1]]   # -sigma1*1.1 + sigma2 <= 0
            b_ineq = [0, 0]
            constr_kinem =  {"A_ineq": A_ineq, "b_ineq": b_ineq}

        NOTE: When possible, it is more efficient to set equality constraints
        using the ``tied`` keyword, instead of setting ``A_eq`` and ``b_eq`` in
        ``constr_kinem``.
    clean:
        Set this keyword to use the iterative sigma clipping method described
        in Section 2.1 of `Cappellari et al. (2002)
        <http://ui.adsabs.harvard.edu/abs/2002ApJ...578..787C>`_.
        This is useful to remove from the fit unmasked bad pixels, residual gas
        emissions or cosmic rays.

        IMPORTANT: This is recommended *only* if a reliable estimate of the
        ``noise`` spectrum is available. See also note below for ``.chi2``.
    degree:
        Degree of the *additive* Legendre polynomial used to correct the
        template continuum shape during the fit (default: 4).
        Set ``degree=-1`` not to include any additive polynomial.
    fixed:
        Boolean vector set to ``True`` where a given kinematic parameter has to
        be held fixed with the value given in ``start``. This is an array, or
        list, with the same dimensions as ``start``.

        EXAMPLE: We want to fit two kinematic components, with 4 moments for
        the first component and 2 for the second. In this case::

            moments = [4, 2]
            start = [[V1, sigma1, h3, h4], [V2, sigma2]]

        then we can held fixed e.g. the sigma (only) of both components using::

            fixed = [[0, 1, 0, 0], [0, 1]]

        NOTE: Setting a negative ``moments`` for a kinematic component is
        entirely equivalent to setting ``fixed = 1`` for all parameters of the
        given kinematic component. In other words::

            moments = [-4, 2]

        is equivalent to::

            moments = [4, 2]
            fixed = [[1, 1, 1, 1], [0, 0]]

    fraction:
        This keyword allows one to fix the ratio between the first two
        kinematic components. This is a scalar defined as follows::

            fraction = np.sum(weights[component == 0])
                     / np.sum(weights[component < 2])

        This is useful e.g. to try to kinematically decompose bulge and disk.

        The remaining kinematic components (``component > 1``) are left free,
        and this allows, for example, to still include gas emission line
        components.
        More general linear constraints, for multiple kinematic components at
        the same time, can be specified using the more general and flexible
        ``constr_templ`` keyword.
    ftol:
        Fractional tolerance for stopping the non-linear minimization (default
        1e-4).
    gas_component:
        Boolean vector, of the same size as ``component``, set to ``True``
        where the given ``component`` describes a gas emission line. If given,
        ``pPXF`` provides the ``pp.gas_flux`` and ``pp.gas_flux_error`` in
        output.

        EXAMPLE: In the common situation where ``component = 0`` are stellar
        templates and the rest are gas emission lines, one will set::

            gas_component = component > 0

        This keyword is also used to plot the gas lines with a different color.
    gas_names:
        String array specifying the names of the emission lines (e.g.
        ``gas_names=["Hbeta", "[OIII]",...]``, one per gas line. The length of
        this vector must match the number of nonzero elements in
        ``gas_component``. This vector is only used by ``pPXF`` to print the
        line names on the console.
    gas_reddening:
        Set this keyword to an initial estimate of the gas reddening
        ``E(B-V) >= 0`` to fit a positive gas reddening together with the
        kinematics and the templates. This reddening is applied only to the gas
        templates, namely to the templates with the corresponding element of
        ``gas_component=True``. The fit assumes by default the extinction curve
        of `Calzetti et al. (2000)
        <http://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C>`_ but any other
        prescription can be passed via the ``reddening_func`` keyword.
    goodpixels:
        Integer vector containing the indices of the good pixels in the
        ``galaxy`` spectrum (in increasing order). Only these spectral pixels
        are included in the fit.

        IMPORTANT: in all likely situations this keyword *has* to be specified.
    lam:
        When the keyword ``reddening`` or ``gas_reddening`` are used, the user
        has to pass in this keyword a vector with the same dimensions of
        ``galaxy``, giving the restframe wavelength in Angstrom of every pixel
        in the input galaxy spectrum. If one uses my ``log_rebin`` routine to
        rebin the spectrum before the ``pPXF`` fit::

            from ppxf.ppxf_util import log_rebin
            specNew, logLam, velscale = log_rebin(lamRange, galaxy)

        the wavelength can be obtained as ``lam = np.exp(logLam)``.

        When ``lam`` is given, the wavelength is shown in the best-fitting
        plot, instead of the pixels.
    linear:
        Set to ``True`` to keep *all* nonlinear parameters fixed and *only*
        perform a linear fit for the templates and additive polynomials
        weights. The output solution is a copy of the input one and the errors
        are zero.
    linear_method: {'nnls', 'lsq_box', 'lsq_lin', 'cvxopt'} optional
        Method used for the solution of the linear least-squares subproblem
        to fit for the templates weights (default 'lsq_box' fast box-constrained).

        The computational speed of the four alternative linear methods depends
        on the size of the problem, with the default 'lsq_box' generally being
        the fastest without linear inequality constraints. Note that 'lsq_lin'
        is included in ``ppxf``, while 'cvxopt' is an optional external package.
        The 'nnls' option (the only one before v7.0) is generally slower and
        for this reason is now deprecated.

        The inequality constraints in ``constr_templ`` are only supported
        with ``linear_method='lsq_lin'`` or ``linear_method='cvxopt'``.
    mask:
        Boolean vector of length ``galaxy.size`` specifying with ``True`` the
        pixels that should be included in the fit. This keyword is just an
        alternative way of specifying the ``goodpixels``.
    method: {'capfit', 'trf', 'dogbox', 'lm'}, optional.
        Algorithm to perform the non-linear minimization step.
        The default 'capfit' is a novel linearly-constrained non-linear
        least-squares optimization program, which combines the Sequential
        Quadratic Programming and the Levenberg-Marquardt methods.
        For a description of the other methods ('trf', 'dogbox', 'lm'), see the
        documentation of ``scipy.optimize.least_squares``.

        The use of linear constraints with ``constr_kinem`` is only supported
        with the default ``method='capfit'``.
    mdegree:
        Degree of the *multiplicative* Legendre polynomial (with a mean of 1)
        used to correct the continuum shape during the fit (default: 0). The
        zero degree multiplicative polynomial is always included in the fit as
        it corresponds to the weights assigned to the templates. Note that the
        computation time is longer with multiplicative polynomials than with
        the same number of additive polynomials.

        IMPORTANT: Multiplicative polynomials cannot be used when the
        ``reddening`` keyword is set, as they are degenerate with the
        reddening.
    moments:
        Order of the Gauss-Hermite moments to fit. Set this keyword to 4 to
        fit ``[h3, h4]`` and to 6 to fit ``[h3, h4, h5, h6]``. Note that in all
        cases the G-H moments are fitted (non-linearly) *together* with
        ``[V, sigma]``.

        If ``moments=2`` or ``moments`` is not set then only ``[V, sigma]`` are
        fitted.

        If ``moments`` is negative then the kinematics of the given
        ``component`` are kept fixed to the input values.
        NOTE: Setting a negative ``moments`` for a kinematic component is
        entirely equivalent to setting ``fixed = 1`` for all parameters of the
        given kinematic component.

        EXAMPLE: We want to keep fixed ``component = 0``, which has a LOSVD
        described by ``[V, sigma, h3, h4]`` and is modelled with 100 spectral
        templates; At the same time, we fit ``[V, sigma]`` for 
        ``component = 1``, which is described by 5 templates (this situation
        may arise when fitting stellar templates with pre-determined stellar
        kinematics, while fitting the gas emission).
        We should give in input to ``pPXF`` the following parameters::

            component = [0]*100 + [1]*5   # --> [0, 0, ... 0, 1, 1, 1, 1, 1]
            moments = [-4, 2]
            start = [[V, sigma, h3, h4], [V, sigma]]

    plot:
        Set this keyword to plot the best fitting solution and the residuals
        at the end of the fit.

        One can also call the class function ``pp.plot()`` after the call to
        ``pp = ppxf(...)``.
    quiet:
        Set this keyword to suppress verbose output of the best fitting
        parameters at the end of the fit.
    reddening:
        Set this keyword to an initial estimate of the stellar reddening
        ``E(B-V) >= 0`` to fit a positive stellar reddening together with the
        kinematics and the templates. This reddening is applied only to the
        stellar templates, namely to the templates with the corresponding
        element of ``gas_component=False``, or to all templates, if
        ``gas_component`` is not set. The fit assumes by default the extinction
        curve of `Calzetti et al. (2000)`_ but any other prescription can be
        passed via the ``reddening_func`` keyword.

        IMPORTANT: The ``mdegree`` keyword cannot be used when ``reddening`` is
        set.
    regul:
        If this keyword is nonzero, the program applies first or second-order
        linear regularization to the ``weights`` during the ``pPXF`` fit.
        Regularization is done in one, two or three dimensions depending on
        whether the array of ``templates`` has two, three or four dimensions
        respectively.
        Large ``regul`` values correspond to smoother ``weights`` output. When
        this keyword is nonzero the solution will be a trade-off between the
        smoothness of ``weights`` and goodness of fit.

        Section 3.5 of `Cappellari (2017)`_ gives a description of
        regularization.

        When fitting multiple kinematic ``component`` the regularization is
        applied only to the first ``component = 0``, while additional
        components are not regularized. This is useful when fitting stellar
        population together with gas emission lines. In that case, the SSP
        spectral templates must be given first and the gas emission templates
        are given last. In this situation, one has to use the ``reg_dim``
        keyword (below), to give ``pPXF`` the dimensions of the population
        parameters (e.g. ``nAge``, ``nMetal``, ``nAlpha``). A usage example is
        given in the file ``ppxf_example_population_gas_sdss.py``.

        The effect of the regularization scheme is the following:

        * With ``reg_ord=1`` it enforces the numerical first derivatives
          between neighbouring weights (in the 1-dim case) to be equal to
          ``w[j] - w[j+1] = 0`` with an error ``Delta = 1/regul``.

        * With ``reg_ord=2`` it enforces the numerical second derivatives
          between neighboring weights (in the 1-dim case) to be equal to
          ``w[j-1] - 2*w[j] + w[j+1] = 0`` with an error ``Delta = 1/regul``.

        It may be helpful to define ``regul = 1/Delta`` and think of ``Delta``
        as the regularization error.

        IMPORTANT: ``Delta`` needs to be smaller but of the same order of
        magnitude of the typical ``weights`` to play an effect on the
        regularization. One quick way to achieve this is:

        1. Divide the full ``templates`` array by a scalar in such a way that
           the typical template has a median of one::

                templates /= np.median(templates)

        2. Do the same for the input galaxy spectrum::

                galaxy /= np.median(galaxy)

        3. In this situation, a sensible guess for ``Delta`` will be a few
           percent (e.g. ``0.01 --> regul=100``).

        Alternatively, for a more rigorous definition of the parameter
        ``regul``:

        A. Perform an un-regularized fit (``regul=0``) and then rescale the
           input ``noise`` spectrum so that::

                Chi^2/DOF = Chi^2/goodPixels.size = 1.

           This is achieved by rescaling the input ``noise`` spectrum as::

                noise = noise*sqrt(Chi**2/DOF) = noise*sqrt(pp.chi2);

        B. Increase ``regul`` and iteratively redo the ``pPXF`` fit until the
           ``Chi^2`` increases from the unregularized value
           ``Chi^2 = goodPixels.size`` by
           ``DeltaChi^2 = sqrt(2*goodPixels.size)``.

        The derived regularization corresponds to the maximum one still
        consistent with the observations and the derived star formation history
        will be the smoothest (minimum curvature or minimum variation) that is
        still consistent with the observations.
    reg_dim:
        When using regularization with more than one kinematic component (using
        the ``component`` keyword), the regularization is only applied to the
        first one (``component=0``). This is useful to fit the stellar
        population and gas emission together.

        In this situation, one has to use the ``reg_dim`` keyword, to give
        ``pPXF`` the dimensions of the population parameters (e.g. ``nAge``,
        ``nMetal``, ``nAlpha``). One should creates the initial array of
        population templates like e.g.
        ``templates[nPixels, nAge, nMetal, nAlpha]`` and define::

            reg_dim = templates.shape[1:]   # = [nAge, nMetal, nAlpha]

        The array of stellar templates is then reshaped into a 2-dim array as::

            templates = templates.reshape(templates.shape[0], -1)

        and the gas emission templates are appended as extra columns at the
        end. An usage example is given in
        ``ppxf_example_population_gas_sdss.py``.

        When using regularization with a single component (the ``component``
        keyword is not used, or contains identical values), the number of
        population templates along different dimensions (e.g. ``nAge``,
        ``nMetal``, ``nAlpha``) is inferred from the dimensions of the
        ``templates`` array and this keyword is not necessary.
    reg_ord:
        Order of the derivative that is minimized by the regularization.
        The following two rotationally-symmetric estimators are supported:

        * ``reg_ord=1``: minimizes the integral over the weights of the squared
          gradient::

            Grad[w] @ Grad[w].

        * ``reg_ord=2``: minimizes the integral over the weights of the squared
          curvature::

            Laplacian[w]**2.

    sigma_diff:
        Quadratic difference in km/s defined as::

            sigma_diff**2 = sigma_inst**2 - sigma_temp**2

        between the instrumental dispersion of the galaxy spectrum and the
        instrumental dispersion of the template spectra.

        This keyword is useful when the templates have higher resolution than
        the galaxy and they were not convolved to match the instrumental
        dispersion of the galaxy spectrum. In this situation, the convolution
        is done by ``pPXF`` with increased accuracy, using an analytic Fourier
        Transform.
    sky:
        vector containing the spectrum of the sky to be included in the fit, or
        array of dimensions ``sky[nPixels, nSky]`` containing different sky
        spectra to add to the model of the observed ``galaxy`` spectrum. The
        ``sky`` has to be log-rebinned as the ``galaxy`` spectrum and needs to
        have the same number of pixels.

        The sky is generally subtracted from the data before the ``pPXF`` fit.
        However, for observations very heavily dominated by the sky spectrum,
        where a very accurate sky subtraction is critical, it may be useful
        *not* to subtract the sky from the spectrum, but to include it in the
        fit using this keyword.
    templates_rfft:
        When calling ``pPXF`` many times with an identical set of templates,
        one can use this keyword to pass the real FFT of the templates,
        computed in a previous ``pPXF`` call, stored in the
        ``pp.templates_rfft`` attribute. This keyword mainly exists to show
        that there is no need for it...

        IMPORTANT: Use this keyword only if you understand what you are doing!
    tied:
        A list of string expressions. Each expression "ties" the parameter to
        other free or fixed parameters.  Any expression involving constants and
        the parameter array ``p[j]`` are permitted. Since they are totally
        constrained, tied parameters are considered to be fixed; no errors are
        computed for them.

        This is an array, or list of arrays, with the same dimensions as
        ``start``. In practice, for every element of ``start`` one needs to
        specify either an empty string ``''`` implying that the parameter is
        free, or a string expression involving some of the variables ``p[j]``,
        where ``j`` represents the index of the flattened list of kinematic 
        parameters.

        EXAMPLE: We want to fit three kinematic components, with 4 moments for
        the first component and 2 moments for the second and third (e.g. stars
        and two gas components). In this case::

            moments = [4, 2, 2]
            start = [[V1, sigma1, 0, 0], [V2, sigma2], [V3, sigma3]]

        then we can force the equality constraint ``V2 = V3`` as follows::

            tied = [['', '', '', ''], ['', ''], ['p[4]', '']]

        or we can force the equality constraint ``sigma2 = sigma3`` as
        follows::

            tied = [['', '', '', ''], ['', ''], ['', 'p[5]']]

        NOTE: One could in principle use the ``tied`` keyword to completely tie
        the LOSVD of two kinematic components. However, this same effect is more
        efficient achieved by assigning them to the same kinematic component
        using the ``component`` keyword.
    trig:
        Set ``trig=True`` to use trigonometric series as an alternative to
        Legendre polynomials, for both the additive and multiplicative
        polynomials. When ``trig=True`` the fitted series below has
        ``N = degree/2`` or ``N = mdegree/2``::

            poly = A_0 + sum_{n=1}^{N} [A_n*cos(n*th) + B_n*sin(n*th)]

        IMPORTANT: The trigonometric series has periodic boundary conditions.
        This is sometimes a desirable property, but this expansion is not as
        flexible as the Legendre polynomials.
    velscale_ratio:
        Integer. Gives the integer ``ratio >= 1`` between the ``velscale`` of
        the ``galaxy`` and the ``templates``. When this keyword is used, the
        templates are convolved by the LOSVD at their native resolution, and
        only subsequently are integrated over the pixels and fitted to
        ``galaxy``. This keyword is generally unnecessary and mostly useful for
        testing.

        Note that in realistic situations the uncertainty in the knowledge and
        variations of the intrinsic line-spread function becomes the limiting
        factor in recovering the LOSVD well below ``velscale``.
    vsyst:
        Reference velocity in ``km/s`` (default 0). The input initial guess and
        the output velocities are measured with respect to this velocity. This
        keyword is generally used to account for the difference in the starting
        wavelength of the templates and the galaxy spectrum as follows::

            vsyst = c*np.log(wave_temp[0]/wave_gal[0])

        The value assigned to this keyword is *crucial* for the two-sided
        fitting. In this case ``vsyst`` can be determined from a previous
        normal one-sided fit to the galaxy velocity profile. After that initial
        fit, ``vsyst`` can be defined as the measured velocity at the galaxy
        center. More accurately ``vsyst`` is the value which has to be
        subtracted to obtain a nearly anti-symmetric velocity profile at the
        two opposite sides of the galaxy nucleus.

        IMPORTANT: this value is generally *different* from the systemic
        velocity one can get from the literature. Do not try to use that!

    Output Parameters
    -----------------

    Stored as attributes of the ``pPXF`` class:

    .apoly:
        Vector with the best fitting additive polynomial.
    .bestfit:
        Vector with the best fitting model for the galaxy spectrum.
        This is a linear combination of the templates, convolved with the best
        fitting LOSVD, multiplied by the multiplicative polynomials and
        with subsequently added polynomial continuum terms or sky components.
    .chi2:
        The reduced ``chi^2`` (namely ``chi^2/DOF``) of the fit
        (where ``DOF ~ pp.goodpixels.size``).

        IMPORTANT: if ``Chi^2/DOF`` is not ~1 it means that the errors are not
        properly estimated, or that the template is bad and it is *not* safe to
        set the ``clean`` keyword.
    .error:
        This variable contains a vector of *formal* uncertainty (``1*sigma``)
        for the fitted parameters in the output vector ``sol``. This option can
        be used when speed is essential, to obtain an order of magnitude
        estimate of the uncertainties, but we *strongly* recommend to run
        bootstrapping simulations to obtain more reliable errors. In fact, these
        errors can be severely underestimated in the region where the penalty
        effect is most important (``sigma < 2*velscale``).

        These errors are meaningless unless ``Chi^2/DOF ~ 1``. However if one
        *assumes* that the fit is good, a corrected estimate of the errors is::

            error_corr = error*sqrt(chi^2/DOF) = pp.error*sqrt(pp.chi2).

        IMPORTANT: when running Monte Carlo simulations to determine the error,
        the penalty (``bias``) should be set to zero, or better to a very small
        value. See Section 3.4 of `Cappellari & Emsellem (2004)`_ for an
        explanation.
    .gas_bestfit:
        If ``gas_component`` is not ``None``, this attribute returns the
        best-fitting gas spectrum alone. The stellar spectrum alone can be
        computed as ``stellar_spectrum = pp.bestfit - pp.gas_bestfit``
    .gas_flux:
        Vector with the integrated flux (in counts) of all lines set as
        ``True`` in the input ``gas_component`` keyword. This is the flux of
        individual gas templates, which may include multiple lines.
        This implies that, if a gas template describes a doublet, the flux is
        that of both lines. If the Balmer series is input as a single template,
        this is the flux of the entire series.

        The returned fluxes are not corrected in any way and in particular, no
        reddening correction is applied. In other words, the returned
        ``.gas_flux`` should be unchanged, within the errors, regardless of
        whether reddening or multiplicative polynomials were fitted by ``pPXF``
        or not. The fluxes are just raw values as one could measure (with lower
        accuracy) by summing the pixels, within the given gas lines, on the
        continuum-subtracted input galaxy spectrum.

        IMPORTANT: ``pPXF`` makes no assumptions about the input flux units:
        The returned ``.gas_flux`` has the same units as the quantity one would
        obtain by just summing the values of the spectral pixels within the gas
        emission. This implies that, if the spectrum is in units of 
        ``erg/(cm^2 s A)``, the gas flux returned by ``pPXF`` should be
        multiplied by the pixel size in Angstrom at the line wavelength to
        obtain the integrated line flux in units of ``erg/(cm^2 s)``.

        NOTE: If there is no gas reddening and each input gas templates was
        normalized to ``sum = 1``, then 
        ``pp.gas_flux = pp.weights[pp.gas_component]``.

        When a gas template is identically zero within the fitted region, then
        ``pp.gas_flux = pp.gas_flux_error = np.nan``. The corresponding
        components of ``pp.gas_zero_template`` are set to ``True``. These
        ``np.nan`` values are set at the end of the calculation to flag the
        undefined values. These flags generally indicate that some of the gas
        templates passed to ``pPXF`` consist of gas emission lines that fall
        outside the fitted wavelength range or within a masked spectral region.
        These ``np.nan`` do *not* indicate numerical issues with the actual
        ``pPXF`` calculation and the rest of the ``pPXF`` output is reliable.
    .gas_flux_error:
        *Formal* uncertainty (``1*sigma``) for the quantity ``pp.gas_flux``, in
        the same units as the gas fluxes.

        This error is approximate as it ignores the covariance between the gas
        flux and any non-linear parameter. Bootstrapping can be used for more
        accurate errors.

        These errors are meaningless unless ``Chi^2/DOF ~ 1``. However if one
        *assumes* that the fit is good, a corrected estimate of the errors is::

            gas_flux_error_corr = gas_flux_error*sqrt(chi^2/DOF)
                                = pp.gas_flux_error*sqrt(pp.chi2).

    .gas_mpoly:
        vector with the best-fitting gas reddening curve.
    .gas_reddening:
        Best fitting ``E(B-V)`` value if the ``gas_reddening`` keyword is set.
        This is especially useful when the Balmer series is input as a single
        template with an assumed theoretically predicted decrement e.g. using
        ``emission_lines(..., tie_balmer=True)`` in ``ppxf.ppxf_util`` to
        compute the gas templates.
    .gas_zero_template:
        vector of size ``gas_component.sum()`` set to ``True`` where
        the gas template was identically zero within the fitted region.
        For those gas components ``pp.gas_flux = pp.gas_flux_error = np.nan``.
        These flags generally indicate that some of the gas templates passed to
        ``pPXF`` consist of gas emission lines that fall outside the fitted
        wavelength range or within a masked spectral region.
    .goodpixels:
        Integer vector containing the indices of the good pixels in the fit.
        This vector is a copy of the input ``goodpixels`` if ``clean = False``
        otherwise it will be updated by removing the detected outliers.
    .matrix:
        Prediction ``matrix[nPixels, degree+nTemplates]`` of the linear system.

        ``pp.matrix[nPixels, :degree]`` contains the additive polynomials if
        ``degree >= 0``.

        ``pp.matrix[nPixels, degree:]`` contains the templates convolved by the
        LOSVD, and multiplied by the multiplicative polynomials if 
        ``mdegree > 0``.
    .mpoly:
        Best fitting multiplicative polynomial (or reddening curve when
        ``reddening`` is set).
    .mpolyweights:
        This is largely superseded by the ``.mpoly`` attribute above.

        When ``mdegree > 0`` this contains in output the coefficients of the
        multiplicative Legendre polynomials of order ``1, 2,... mdegree``.
        The polynomial can be explicitly evaluated as::

            from numpy.polynomial import legendre
            x = np.linspace(-1, 1, len(galaxy))
            mpoly = legendre.legval(x, np.append(1, pp.mpolyweights))

        When ``trig = True`` the polynomial is evaluated as::

            mpoly = pp.trigval(x, np.append(1, pp.mpolyweights))

    .polyweights:
        This is largely superseded by the ``.apoly`` attribute above.

        When ``degree >= 0`` contains the weights of the additive Legendre
        polynomials of order ``0, 1,... degree``. The best fitting additive
        polynomial can be explicitly evaluated as::

            from numpy.polynomial import legendre
            x = np.linspace(-1, 1, len(galaxy))
            apoly = legendre.legval(x, pp.polyweights)

        When ``trig=True`` the polynomial is evaluated as::

            apoly = pp.trigval(x, pp.polyweights)

        When doing a two-sided fitting (see help for ``galaxy`` parameter), the
        additive polynomials are allowed to be different for the left and right
        spectrum. In that case, the output weights of the additive polynomials
        alternate between the first (left) spectrum and the second (right)
        spectrum.
    .reddening:
        Best fitting ``E(B-V)`` value if the ``reddening`` keyword is set.
    .sol:
        Vector containing in output the parameters of the kinematics.

        * If ``moments=2`` this contains ``[Vel, Sigma]``
        * If ``moments=4`` this contains ``[Vel, Sigma, h3, h4]``
        * If ``moments=N`` this contains ``[Vel, Sigma, h3,... hN]``

        When fitting multiple kinematic ``component``, ``pp.sol`` contains a
        list with the solution for all different components, one after the
        other, sorted by ``component``: ``pp.sol = [sol1, sol2,...]``.

        ``Vel`` is the velocity, ``Sigma`` is the velocity dispersion, 
        ``h3 - h6`` are the Gauss-Hermite coefficients. The model parameters
        are fitted simultaneously.

        IMPORTANT: The precise relation between the output ``pPXF`` velocity
        and redshift is ``Vel = c*np.log(1 + z)``. See Section 2.3 of
        `Cappellari (2017)`_ for a detailed explanation.

        These are the default safety limits on the fitting parameters
        (they can be changed using the ``bounds`` keyword):

        * Vel is constrained to be +/-2000 km/s from the first input guess
        * ``velscale/100 < Sigma < 1000`` km/s
        * ``-0.3 < [h3, h4, ...] < 0.3``  (extreme value for real galaxies)

        In the case of two-sided LOSVD fitting the output values refer to the
        first input galaxy spectrum, while the second spectrum will have by
        construction kinematics parameters ``[-Vel, Sigma, -h3, h4, -h5, h6]``.
        If ``vsyst`` is nonzero (as required for two-sided fitting), then the
        output velocity is measured with respect to ``vsyst``.
    .status:
        Contains the output status of the optimization. Positive values
        generally represent success (the status is defined as in
        ``scipy.optimize.least_squares``).
    .weights:
        Receives the value of the weights by which each template was
        multiplied to best fit the galaxy spectrum. The optimal template can be
        computed with an array-vector multiplication::

            bestemp = templates @ weights

        These weights do not include the weights of the additive polynomials
        which are separately stored in ``pp.polyweights``.

        When the ``sky`` keyword is used ``weights[:nTemplates]`` contains the
        weights for the templates, while ``weights[nTemplates:]`` gives the
        ones for the sky. In that case the best fitting galaxy template and sky
        are given by::

            bestemp = templates @ weights[:nTemplates]
            bestsky = sky @ weights[nTemplates:]

        When doing a two-sided fitting (see help for ``galaxy`` parameter)
        *together* with the ``sky`` keyword, the sky weights are allowed to be
        different for the left and right spectrum. In that case the output sky
        weights alternate between the first (left) spectrum and the second
        (right) spectrum.

    How to Set Regularization
    -------------------------

    The ``pPXF`` routine can give sensible quick results with the default
    ``bias`` parameter, however, like in any penalized/filtered/regularized
    method, the optimal amount of penalization generally depends on the problem
    under study.

    The general rule here is that the penalty should leave the line-of-sight
    velocity-distribution (LOSVD) virtually unaffected, when it is well sampled
    and the signal-to-noise ratio (``S/N``) is sufficiently high.

    EXAMPLE: If you expect a LOSVD with up to a high ``h4 ~ 0.2`` and your
    adopted penalty (``bias``) biases the solution towards a much lower 
    ``h4 ~ 0.1``, even when the measured ``sigma > 3*velscale`` and the S/N is
    high, then you are *misusing* the ``pPXF`` method!

    THE RECIPE: The following is a simple practical recipe for a sensible
    determination of the penalty in ``pPXF``:

    1. Choose a minimum ``(S/N)_min`` level for your kinematics extraction and
       spatially bin your data so that there are no spectra below
       ``(S/N)_min``;
    2. Perform a fit of your kinematics *without* penalty (keyword ``bias=0``).
       The solution will be noisy and may be affected by spurious solutions,
       however, this step will allow you to check the expected mean ranges in
       the Gauss-Hermite parameters ``[h3, h4]`` for the galaxy under study;
    3. Perform a Monte Carlo simulation of your spectra, following e.g. the
       included ``ppxf_example_simulation.py`` routine. Adopt as ``S/N`` in the
       simulation the chosen value ``(S/N)_min`` and as input ``[h3, h4]`` the
       maximum representative values measured in the non-penalized ``pPXF`` fit
       of the previous step;
    4. Choose as the penalty (``bias``) the *largest* value such that, for
       ``sigma > 3*velscale``, the mean difference delta between the output 
       ``[h3, h4]`` and the input ``[h3, h4]`` is well within (e.g. 
       ``delta ~ rms/3``) the rms scatter of the simulated values (see an
       example in Fig. 2 of `Emsellem et al. 2004
       <http://ui.adsabs.harvard.edu/abs/2004MNRAS.352..721E>`_).

    Problems with Your First Fit?
    -----------------------------

    Common problems with your first ``pPXF`` fit are caused by incorrect
    wavelength ranges or different velocity scales between galaxy and
    templates. To quickly detect these problems try to overplot the (log
    rebinned) galaxy and the template just before calling the ``pPXF``
    procedure.

    You can use something like the following Python lines while adjusting the
    smoothing window and the pixels shift. If you cannot get a rough match
    by eye it means something is wrong and it is unlikely that ``pPXF``
    (or any other program) will find a good match:

    .. code-block:: python

        import numpy as np
        import matplotlib.pyplot as plt
        import scipy.ndimage

        sigma = 2       # Velocity dispersion in pixels
        shift = -20     # Velocity shift in pixels
        template = np.roll(ndimage.gaussian_filter1d(template, sigma), shift)
        plt.plot(galaxy, 'k')
        plt.plot(template*np.median(galaxy)/np.median(template), 'r')

    ###########################################################################
    """
    def __init__(self, templates, galaxy, noise, velscale, start,
                 bias=None, bounds=None, clean=False, component=0,
                 constr_templ={}, constr_kinem={}, degree=4, fixed=None,
                 fraction=None, ftol=1e-4, gas_component=None, gas_names=None,
                 gas_reddening=None, goodpixels=None, lam=None, linear=False,
                 linear_method='lsq_box', mask=None, method='capfit', mdegree=0,
                 moments=2, plot=False, quiet=False, reddening=None,
                 reddening_func=None, reg_dim=None, reg_ord=2, regul=0,
                 sigma_diff=0, sky=None, templates_rfft=None, tied=None,
                 trig=False, velscale_ratio=1, vsyst=0, x0=None):

        # Do extensive checking of possible input errors
        self.galaxy = galaxy
        self.nspec = galaxy.ndim     # nspec=2 for reflection-symmetric LOSVD
        self.npix = galaxy.shape[0]  # total pixels in the galaxy spectrum
        self.noise = noise
        self.clean = clean
        self.fraction = fraction
        self.ftol = ftol
        self.gas_reddening = gas_reddening
        self.degree = max(degree, -1)
        self.mdegree = max(mdegree, 0)
        self.method = method
        self.quiet = quiet
        self.sky = sky
        self.vsyst = vsyst/velscale
        self.regul = regul
        self.lam = lam
        self.nfev = 0
        self.reddening = reddening
        self.reg_dim = np.asarray(reg_dim)
        self.reg_ord = reg_ord
        self.templates = templates.reshape(templates.shape[0], -1)
        self.npix_temp, self.ntemp = self.templates.shape
        self.sigma_diff = sigma_diff/velscale
        self.status = 0   # Initialize status as failed
        self.velscale = velscale
        self.velscale_ratio = velscale_ratio
        self.tied = tied
        self.gas_flux = self.gas_flux_error = self.gas_bestfit = None
        self.linear_method = linear_method
        self.x0 = x0  # Initialization for linear solution

        if reddening_func is None:
            self.reddening_func = reddening_cal00
        else:
            assert callable(reddening_func), "`reddening_func` must be callable"
            self.reddening_func = reddening_func

        if method != 'capfit':
            assert method in ['trf', 'dogbox', 'lm'], \
                "`method` must be 'capfit', 'trf', 'dogbox' or 'lm'"
            assert tied is None, "Parameters can only be tied with method='capfit'"
            assert fixed is None, "Parameters can only be fixed with method='capfit'"
            if method == 'lm':
                assert bounds is None, "Bounds not supported with method='lm'"

        assert linear_method in ['nnls', 'lsq_box', 'lsq_lin', 'cvxopt'], \
            "`linear_method` must be 'nnls', 'lsq_box', 'lsq_lin' or 'cvxopt'"

        if (linear_method == 'cvxopt') or (constr_kinem != {}):
            assert cvxopt_installed, \
                "To use `linear_method`='cvxopt' or `constr_kinem` " \
                "the cvxopt package must be installed"

        if trig:
            assert degree < 0 or degree % 2 == 0, \
                "`degree` must be even with trig=True"
            assert mdegree < 0 or mdegree % 2 == 0, \
                "`mdegree` must be even with trig=True"
            self.polyval = trigval
            self.polyvander = trigvander
        else:
            self.polyval = legendre.legval
            self.polyvander = legendre.legvander

        assert isinstance(velscale_ratio, int), "VELSCALE_RATIO must be an integer"

        component = np.atleast_1d(component)
        assert np.issubdtype(component.dtype, np.integer), "COMPONENT must be integers"

        if component.size == 1 and self.ntemp > 1:  # component is a scalar
            # all templates have the same LOSVD
            self.component = np.zeros(self.ntemp, dtype=int)
        else:
            assert component.size == self.ntemp, \
                "There must be one kinematic COMPONENT per template"
            self.component = component

        if gas_component is None:
            self.gas_component = np.zeros(self.component.size, dtype=bool)
            self.gas_any = False
        else:
            self.gas_component = np.asarray(gas_component)
            assert self.gas_component.dtype == bool, \
                "`gas_component` must be boolean"
            assert self.gas_component.size == component.size, \
                "`gas_component` and `component` must have the same size"
            assert np.any(gas_component), "`gas_component` must be nonzero"
            if gas_names is None:
                self.gas_names = np.full(np.sum(gas_component), "Unknown")
            else:
                assert self.gas_component.sum() == len(gas_names), \
                    "There must be one name per gas emission line template"
                self.gas_names = gas_names
            self.gas_any = True

        tmp = np.unique(component)
        self.ncomp = tmp.size
        assert np.array_equal(tmp, np.arange(self.ncomp)), \
            "COMPONENT must range from 0 to NCOMP-1"

        if fraction is not None:
            assert 0 < fraction < 1, "Must be `0 < fraction < 1`"
            assert self.ncomp >= 2, \
                "At least 2 COMPONENTs are needed with FRACTION keyword"

        if regul > 0 and reg_dim is None:
            assert self.ncomp == 1, \
                "REG_DIM must be specified with more than one component"
            self.reg_dim = np.asarray(templates.shape[1:])
            
        assert reg_ord in [1, 2], "`reg_ord` must be 1 or 2"

        moments = np.atleast_1d(moments)
        if moments.size == 1:
            # moments is scalar: all LOSVDs have same number of G-H moments
            moments = np.full(self.ncomp, moments, dtype=int)

        self.fixall = moments < 0  # negative moments --> keep entire LOSVD fixed
        self.moments = np.abs(moments)

        assert tmp.size == self.moments.size, \
            "MOMENTS must be an array of length NCOMP"

        if sky is not None:
            assert sky.shape[0] == galaxy.shape[0], \
                "GALAXY and SKY must have the same size"
            self.sky = sky.reshape(sky.shape[0], -1)

        assert galaxy.ndim < 3 and noise.ndim < 3, \
            "Wrong GALAXY or NOISE input dimensions"

        if noise.ndim == 2 and noise.shape[0] == noise.shape[1]:
            # NOISE is a 2-dim covariance matrix
            assert noise.shape[0] == galaxy.shape[0], \
                "Covariance Matrix must have size npix*npix"
            # Cholesky factor of symmetric, positive-definite covariance matrix
            noise = linalg.cholesky(noise, lower=1)
            # Invert Cholesky factor
            self.noise = linalg.solve_triangular(noise, np.eye(noise.shape[0]), lower=1)
        else:   # NOISE is an error spectrum
            assert galaxy.shape == noise.shape, \
                "GALAXY and NOISE must have the same size"
            assert np.all((noise > 0) & np.isfinite(noise)), \
                "NOISE must be a positive vector"
            if self.nspec == 2:   # reflection-symmetric LOSVD
                self.noise = self.noise.T.reshape(-1)
                self.galaxy = self.galaxy.T.reshape(-1)

        assert np.all(np.isfinite(galaxy)), 'GALAXY must be finite'

        assert self.npix_temp >= galaxy.shape[0]*velscale_ratio, \
            "TEMPLATES length cannot be smaller than GALAXY"

        if reddening is not None:
            assert lam is not None, "LAM must be given with REDDENING keyword"
            assert mdegree < 1, "MDEGREE cannot be used with REDDENING keyword"

        if gas_reddening is not None:
            assert lam is not None, "LAM must be given with GAS_REDDENING keyword"
            assert self.gas_any, "GAS_COMPONENT must be nonzero with GAS_REDDENING keyword"

        if lam is not None:
            assert lam.shape == galaxy.shape, \
                "GALAXY and LAM must have the same size"

        if mask is not None:
            assert mask.dtype == bool, "MASK must be a boolean vector"
            assert mask.shape == galaxy.shape, \
                "GALAXY and MASK must have the same size"
            assert goodpixels is None, \
                "GOODPIXELS and MASK cannot be used together"
            goodpixels = np.flatnonzero(mask)

        if goodpixels is None:
            self.goodpixels = np.arange(galaxy.shape[0])
        else:
            assert np.all(np.diff(goodpixels) > 0), \
                "GOODPIXELS is not monotonic or contains duplicated values"
            assert goodpixels[0] >= 0 and goodpixels[-1] < galaxy.shape[0], \
                "GOODPIXELS are outside the data range"
            self.goodpixels = goodpixels

        if bias is None:
            # Cappellari & Emsellem (2004) pg.144 left
            self.bias = 0.7*np.sqrt(500./self.goodpixels.size)
        else:
            self.bias = bias

        start1 = [start] if self.ncomp == 1 else list(start)
        assert np.all([hasattr(a, "__len__") and 2 <= len(a) <= b
                       for (a, b) in zip(start1, self.moments)]), \
            "START must be a list/array of vectors [start1, start2,...] " \
            "with each vector satisfying 2 <= len(START[j]) <= MOMENTS[j]"
        assert len(start1) == self.ncomp, \
            "There must be one START per COMPONENT"

        # The lines below deal with the possibility for the input
        # gas templates to be identically zero in the fitted region
        self.gas_any_zero = False
        if self.gas_any:
            vmed = np.median([a[0] for a in start1])
            dx = int(np.round((vsyst + vmed)/velscale))  # Approximate velocity shift
            gas_templates = self.templates[:, self.gas_component]
            tmp = rebin(gas_templates, velscale_ratio)
            gas_peak = np.max(np.abs(tmp), axis=0)
            tmp = np.roll(tmp, dx, axis=0)
            good_peak = np.max(np.abs(tmp[self.goodpixels, :]), axis=0)
            self.gas_zero_template = good_peak <= gas_peak/1e3
            if np.any(self.gas_zero_template):
                self.gas_any_zero = True
                gas_ind = np.flatnonzero(self.gas_component)
                self.gas_zero_ind = gas_ind[self.gas_zero_template]
                if not quiet:
                    print("Warning: Some gas templates are identically zero in "
                          "the fitted range, the gas emissions likely fall in a "
                          "masked region or outside the fitted wavelength range")
                    print(self.gas_zero_template)

        # Pad with zeros when `start[j]` has fewer elements than `moments[j]`
        for j, (st, mo) in enumerate(zip(start1, self.moments)):
            st = np.asarray(st, dtype=float)   # Make sure starting guess is float
            start1[j] = np.pad(st, [0, mo - len(st)])

        if bounds is not None:
            if self.ncomp == 1:
                bounds = [bounds]
            assert list(map(len, bounds)) == list(map(len, start1)), \
                "BOUNDS and START must have the same shape"
            assert np.all([hasattr(c, "__len__") and len(c) == 2
                           for a in bounds for c in a]), \
                "All BOUNDS must have two elements [lb, ub]"

        if fixed is not None:
            if self.ncomp == 1:
                fixed = [fixed]
            assert list(map(len, fixed)) == list(map(len, start1)), \
                "FIXED and START must have the same shape"

        if tied is not None:
            if self.ncomp == 1:
                tied = [tied]
            assert list(map(len, tied)) == list(map(len, start1)), \
                "TIED and START must have the same shape"

        self.set_linear_constraints(constr_templ, constr_kinem, method)

        if galaxy.ndim == 2:
            # two-sided fitting of LOSVD
            assert vsyst != 0, "VSYST must be defined for two-sided fitting"
            self.goodpixels = np.append(self.goodpixels, galaxy.shape[0] + self.goodpixels)

        self.npad = 2**int(np.ceil(np.log2(self.templates.shape[0])))
        if templates_rfft is None:
            # Pre-compute FFT of real input of all templates
            self.templates_rfft = np.fft.rfft(self.templates, self.npad, axis=0)
        else:
            self.templates_rfft = templates_rfft

        # Convert velocity from km/s to pixels
        for j, st in enumerate(start1):
            start1[j][:2] = st[:2]/velscale

        if (np.all(moments < 0) or np.all(fixed)) and (mdegree <= 0) \
                and (reddening is None) and (gas_reddening is None):
            linear = True

        if linear:
            assert mdegree <= 0, "Must be `mdegree` <= 0 with `linear`=True"
            params = np.concatenate(start1)  # Flatten list
            if reddening is not None:
                params = np.append(params, reddening)
            if gas_reddening is not None:
                params = np.append(params, gas_reddening)
            perror = np.zeros_like(params)
            self.method = 'linear'
            self.status = 1   # Status irrelevant for linear fit
            self.njev = 0     # Jacobian is not evaluated
        else:
            params, perror = self.nonlinear_fit(start1, bounds, fixed, tied, clean)

        self.bias = 0   # Evaluate residuals without bias
        err = self.linear_fit(params)
        self.ndof = err.size - (perror > 0).sum()
        self.chi2 = np.sum(err**2)/self.ndof   # Chi**2/DOF
        self.format_output(params, perror)
        if plot:   # Plot final data-model comparison if required.
            self.plot()

################################################################################

    def set_linear_constraints(self, constr_templ, constr_kinem, method):

        npoly = (self.degree + 1)*self.nspec
        self.nsky = 0 if self.sky is None else self.sky.shape[1]*self.nspec
        ncols = npoly + self.ntemp + self.nsky
        self.ngh = self.moments.sum()    # Parameters of the LOSVD only
        self.npars = self.ngh + self.mdegree*self.nspec
        if self.reddening is not None:
            self.npars += 1
        if self.gas_reddening is not None:
            self.npars += 1

        # See Equation (30) of Cappellari (2017)
        self.A_eq_templ = self.b_eq_templ = None
        if self.fraction is not None:
            self.A_eq_templ = np.zeros((1, ncols))
            ff = self.A_eq_templ[0, npoly: npoly + self.ntemp]
            ff[self.component == 0] = self.fraction - 1
            ff[self.component == 1] = self.fraction
            self.b_eq_templ = np.zeros(1)

        # Constrain identically-zero gas templates to zero weight
        if self.gas_any_zero:
            nz = self.gas_zero_ind.size
            A_gas_zero = np.zeros((nz, ncols))
            ff = A_gas_zero[:, npoly: npoly + self.ntemp]
            ff[np.arange(nz), self.gas_zero_ind] = 1
            b_gas_zero = np.zeros(nz)
            if self.A_eq_templ is None:
                self.A_eq_templ, self.b_eq_templ = A_gas_zero, b_gas_zero
            else:
                self.A_eq_templ = np.vstack([self.A_eq_templ, A_gas_zero])
                self.b_eq_templ = np.append(self.b_eq_templ, b_gas_zero)

        pos = self.ntemp + self.nsky
        self.A_ineq_templ = np.pad(-np.eye(pos), [(0, 0), (npoly, 0)])
        self.b_ineq_templ = np.zeros(pos)   # Positivity
        if constr_templ != {}:
            assert isinstance(constr_templ, dict), "`constr_templ` must be a dictionary"
            if "A_ineq" in constr_templ:
                assert self.linear_method not in ['lsq_box', 'nnls'], \
                    "`A_ineq` not supported with linear_method='lsq_box' or 'nnls'"
                A_ineq = np.pad(constr_templ["A_ineq"], [(0, 0), (npoly, self.nsky)])
                self.A_ineq_templ = np.vstack([self.A_ineq_templ, A_ineq])
                self.b_ineq_templ = np.append(self.b_ineq_templ, constr_templ["b_ineq"])
            if "A_eq" in constr_templ:
                A_eq = np.pad(constr_templ["A_eq"], [(0, 0), (npoly, self.nsky)])
                b_eq = constr_templ["b_eq"]
                if self.A_eq_templ is None:
                    self.A_eq_templ, self.b_eq_templ = A_eq, b_eq
                else:
                    self.A_eq_templ = np.vstack([self.A_eq_templ, A_eq])
                    self.b_eq_templ = np.append(self.b_eq_templ, b_eq)

        self.A_ineq_kinem = self.b_ineq_kinem = self.A_eq_kinem = self.b_eq_kinem = None
        if constr_kinem != {}:
            assert method == 'capfit', "Linear constraints on kinematics require method='capfit'"
            assert isinstance(constr_kinem, dict), "`constr_kinem` must be a dictionary"
            if "A_ineq" in constr_kinem:
                self.A_ineq_kinem = np.pad(constr_kinem["A_ineq"], [(0, 0), (0, self.npars - self.ngh)])
                self.b_ineq_kinem = constr_kinem["b_ineq"]
            if "A_eq" in constr_kinem:
                self.A_eq_kinem = np.pad(constr_kinem["A_eq"], [(0, 0), (0, self.npars - self.ngh)])
                self.b_eq_kinem = constr_kinem["b_eq"]

################################################################################

    def solve_linear(self, A, b, npoly):

        m, n = A.shape
        if n == 1:  # A is a column vector, not an array
            soluz = ((b @ A)/(A.T @ A))[0]
        elif n == npoly + 1:  # Fitting a single template with polynomials
            soluz = linalg.lstsq(A, b)[0]
        else:  # Fitting multiple templates
            #print (A, b)
            #A = np.nan_to_num(A)
            if self.linear_method == 'lsq_lin':
                soluz = lsq_lin(A, b, self.A_ineq_templ, self.b_ineq_templ,
                               self.A_eq_templ, self.b_eq_templ, x=self.x0).x
                self.x0 = soluz
            elif self.linear_method == 'cvxopt':
                res = lsq_lin_cvxopt(A, b, self.A_ineq_templ, self.b_ineq_templ,
                                    self.A_eq_templ, self.b_eq_templ, initvals=self.x0)
                self.x0 = res.initvals
                soluz = res.x
            else:   # linear_method='lsq_box' or 'nnls'
                if self.A_eq_templ is not None:
                    scale = 1e-4*linalg.norm(self.A_eq_templ)/linalg.norm(A)
                    A = np.vstack([self.A_eq_templ/scale, A])
                    b = np.append(self.b_eq_templ/scale, b)
                if self.linear_method == 'lsq_box':
                    lb = np.zeros(n)
                    lb[:npoly] = -np.inf
                    soluz = lsq_box(A, b, [lb, np.inf], self.x0).x
                    self.x0 = soluz
                else:  # linear_method='nnls'
                    AA = np.hstack([A, -A[:, :npoly]])
                    x = optimize.nnls(AA, b)[0]
                    x[:npoly] -= x[n:]
                    soluz = x[:n]

        return soluz

################################################################################

    def nonlinear_fit(self, start0, bounds0, fixed0, tied0, clean):
        """
        This function implements the procedure described in
        Section 3.4 of Cappellari M., 2017, MNRAS, 466, 798
        http://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C

        """
        # Explicitly specify the step for the numerical derivatives
        # and force safety limits on the fitting parameters.
        #
        # Set [h3, h4, ...] and mult. polynomials to zero as initial guess
        # and constrain -0.3 < [h3, h4, ...] < 0.3
        start = np.zeros(self.npars)
        fixed = np.full(self.npars, False)
        tied = np.full(self.npars, '', dtype=object)
        step = np.full(self.npars, 0.001)
        bounds = np.tile([-0.3, 0.3], (self.npars, 1))

        p = 0
        for j, st in enumerate(start0):
            if bounds0 is None:
                bn = [st[0] + np.array([-2e3, 2e3])/self.velscale,  # V bounds
                      [0.01, 1e3/self.velscale]]                # sigma bounds
            else:
                bn = np.array(bounds0[j][:2], dtype=float)/self.velscale
            for k in range(self.moments[j]):
                if self.fixall[j]:  # Negative moment --> keep entire LOSVD fixed
                    fixed[p + k] = True
                elif fixed0 is not None:  # Keep individual LOSVD parameters fixed
                    fixed[p + k] = fixed0[j][k]
                if tied0 is not None:
                    tied[p + k] = tied0[j][k]
                if k < 2:
                    start[p + k] = st[k].clip(*bn[k])
                    bounds[p + k] = bn[k]
                    step[p + k] = 0.01
                else:
                    start[p + k] = st[k]
                    if bounds0 is not None:
                        bounds[p + k] = bounds0[j][k]
            p += self.moments[j]

        if self.mdegree > 0:
            for q in range(self.mdegree*self.nspec):
                bounds[p + q] = [-1, 1]  # Force <100% corrections
        elif self.reddening is not None:
            start[p] = self.reddening
            bounds[p] = [0., 10.]  # Force positive E(B-V) < 10 mag

        if self.gas_reddening is not None:
            start[-1] = self.gas_reddening
            bounds[-1] = [0., 10.]  # Force positive E(B-V) < 10 mag

        if self.method == 'lm':
            step = 0.01  # only a scalar is supported
            bounds = np.array([-np.inf, np.inf])   # No bounds

        # Here the actual calculation starts.
        # If required, once the minimum is found, clean the pixels deviating
        # more than 3*sigma from the best fit and repeat the minimization
        # until the set of cleaned pixels does not change any more.
        good = self.goodpixels.copy()
        for j in range(5):  # Do at most five cleaning iterations
            self.clean = False  # No cleaning during chi2 optimization
            if self.method == 'capfit':
                res = capfit(
                    self.linear_fit, start, ftol=self.ftol, bounds=bounds.T,
                    abs_step=step, x_scale='jac', tied=tied, fixed=fixed,
                    A_ineq=self.A_ineq_kinem, b_ineq=self.b_ineq_kinem,
                    A_eq=self.A_eq_kinem, b_eq=self.b_eq_kinem)
                perror = res.x_err
            else:
                res = optimize.least_squares(
                    self.linear_fit, start, ftol=self.ftol, bounds=bounds.T,
                    diff_step=step, x_scale='jac', method=self.method)
                perror = cov_err(res.jac)[1]
            params = res.x
            if not clean:
                break
            good_old = self.goodpixels.copy()
            self.goodpixels = good.copy()  # Reset goodpixels
            self.clean = True  # Do cleaning during linear fit
            self.linear_fit(params)
            if np.array_equal(good_old, self.goodpixels):
                break

        self.status = res.status
        self.njev = res.njev

        return params, perror

################################################################################

    def linear_fit(self, pars, in_jac=False):
        """
        This function implements the procedure described in
        Sec.3.3 of Cappellari M., 2017, MNRAS, 466, 798
        http://ui.adsabs.harvard.edu/abs/2017MNRAS.466..798C

        """
        # pars = [vel_1, sigma_1, h3_1, h4_1, ... # Velocities are in pixels.
        #         ...                             # For all kinematic components
        #         vel_n, sigma_n, h3_n, h4_n, ...
        #         m_1, m_2,... m_mdegree          # Multiplicative polynomials
        #         E(B-V)_stars, E(B-V)_gas]       # Reddening

        nspec, npix, ngh = self.nspec, self.npix, self.ngh
        lvd_rfft = losvd_rfft(pars, nspec, self.moments, self.templates_rfft.shape[0],
                              self.ncomp, self.vsyst, self.velscale_ratio, self.sigma_diff)

        # The zeroth order multiplicative term is already included in the
        # linear fit of the templates. The polynomial below has mean of 1.
        # X needs to be within [-1, 1] for Legendre Polynomials
        x = np.linspace(-1, 1, npix)
        if self.mdegree > 0:
            pars_mpoly = pars[ngh:] if self.gas_reddening is None else pars[ngh:-1]
            if nspec == 2:  # Different multiplicative poly for left/right spectra
                mpoly1 = self.polyval(x, np.append(1.0, pars_mpoly[::2]))
                mpoly2 = self.polyval(x, np.append(1.0, pars_mpoly[1::2]))
                mpoly = np.append(mpoly1, mpoly2).clip(0.1)
            else:
                mpoly = self.polyval(x, np.append(1.0, pars_mpoly)).clip(0.1)
        elif self.reddening is not None:
            # Multiplicative polynomials do not make sense when fitting stellar reddening.
            # In that case one has to assume the spectrum is well calibrated.
            mpoly = self.reddening_func(self.lam, pars[ngh])
        else:
            mpoly = None

        if self.gas_reddening is not None:
            gas_mpoly = self.reddening_func(self.lam, pars[-1])
        else:
            gas_mpoly = None

        # This array is used for estimating predictions
        npoly = (self.degree + 1)*nspec  # Number of additive polynomials in fit
        #print (npoly, self.ntemp, self.nsky)
        ncols = npoly + self.ntemp + self.nsky
        c = np.zeros((npix*nspec, ncols))

        if self.degree >= 0:  # Fill first columns of the Design Matrix
            vand = self.polyvander(x, self.degree)
            c[: npix, : npoly//nspec] = vand
            if nspec == 2:
                c[npix :, npoly//nspec : npoly] = vand  # poly for right spectrum

        tmp = np.empty((nspec, self.npix))
        for j, template_rfft in enumerate(self.templates_rfft.T):  # columns loop
            for k in range(nspec):
                tt = np.fft.irfft(template_rfft*lvd_rfft[:, self.component[j], k], self.npad)
                tmp[k, :] = rebin(tt[:self.npix*self.velscale_ratio], self.velscale_ratio)
            c[:, npoly + j] = tmp.ravel()
            if self.gas_component[j]:
                if gas_mpoly is not None:
                    c[:, npoly + j] *= gas_mpoly
            elif mpoly is not None:
                c[:, npoly + j] *= mpoly

        if self.nsky > 0:
            k = npoly + self.ntemp
            c[: npix, k : k + self.nsky//nspec] = self.sky
            if nspec == 2:
                c[npix :, k + self.nsky//nspec : k + self.nsky] = self.sky  # Sky for right spectrum

        #print (self.reg_dim.size, np.prod(self.reg_dim))
        if self.regul > 0:
            if self.reg_ord == 1:
                nr = self.reg_dim.size
                nreg = nr*np.prod(self.reg_dim)
            elif self.reg_ord == 2:
                nreg = np.prod(self.reg_dim)
        else:
            nreg = 0

        # This array is used for the system solution
        #print (npix, nspec, nreg)
        nrows = npix*nspec + nreg
        a = np.zeros((nrows, ncols))

        if self.noise.ndim > 1 and self.noise.shape[0] == self.noise.shape[1]:
            # input NOISE is a npix*npix covariance matrix
            a[:npix*nspec, :] = self.noise @ c
            b = self.noise @ self.galaxy
        else:
            # input NOISE is a 1sigma error vector
            a[:npix*nspec, :] = c/self.noise[:, None] # Weight columns with errors
            b = self.galaxy/self.noise

        if self.regul > 0:
            regularization(a, npoly, npix, nspec, self.reg_dim, self.reg_ord, self.regul)

        # Select the spectral region to fit and solve the over-conditioned system
        # using SVD/BVLS. Use unweighted array for estimating bestfit predictions.
        # Iterate to exclude pixels deviating >3*sigma if clean=True.
        m = 1
        while m > 0:
            if nreg > 0:
                aa = a[np.append(self.goodpixels, np.arange(npix*nspec, nrows)), :]
                bb = np.append(b[self.goodpixels], np.zeros(nreg))
            else:
                aa = a[self.goodpixels, :]
                bb = b[self.goodpixels]
            #print (aa, bb)
            self.weights = self.solve_linear(aa, bb, npoly)
            self.bestfit = c @ self.weights
            if self.noise.ndim > 1 and self.noise.shape[0] == self.noise.shape[1]:
                # input NOISE is a npix*npix covariance matrix
                err = (self.noise @ (self.galaxy - self.bestfit))[self.goodpixels]
            else:
                # input NOISE is a 1sigma error vector
                err = ((self.galaxy - self.bestfit)/self.noise)[self.goodpixels]
            if self.clean:
                w = np.abs(err) < 3  # select residuals smaller than 3*sigma
                m = err.size - w.sum()
                if m > 0:
                    self.goodpixels = self.goodpixels[w]
                    if not self.quiet:
                        print('Outliers:', m)
            else:
                break

        self.matrix = c          # Return LOSVD-convolved templates matrix
        self.mpoly = mpoly
        self.gas_mpoly = gas_mpoly

        # Penalize the solution towards (h3, h4, ...) = 0 if the inclusion of
        # these additional terms does not significantly decrease the error.
        # The lines below implement eq.(8)-(9) in Cappellari & Emsellem (2004)
        if np.any(self.moments > 2) and self.bias > 0:
            D2 = p = 0
            for mom in self.moments:  # loop over kinematic components
                if mom > 2:
                    D2 += np.sum(pars[2 + p : mom + p]**2)  # eq.(8) CE04
                p += mom
            err += self.bias*robust_sigma(err, zero=True)*np.sqrt(D2)  # eq.(9) CE04

        self.nfev += 1

        return err

################################################################################

    def plot(self):
        """
        Produces a plot of the pPXF best fit.
        One can call pp.plot() after pPXF terminates.

        """
        if self.lam is None:
            plt.xlabel("Pixels")
            x = np.arange(self.galaxy.size)
        else:
            plt.xlabel(r"Wavelength (nm)")
            x = self.lam/10

        ll, rr = np.min(x), np.max(x)
        resid = self.galaxy - self.bestfit
        mn = np.min(self.bestfit[self.goodpixels])
        mn -= np.percentile(np.abs(resid[self.goodpixels]), 99)
        mx = np.max(self.bestfit[self.goodpixels])
        resid += mn   # Offset residuals to avoid overlap
        mn1 = np.min(resid[self.goodpixels])
        plt.ylabel("Relative Flux")
        plt.xlim([ll, rr] + np.array([-0.02, 0.02])*(rr - ll))
        plt.ylim([mn1, mx] + np.array([-0.05, 0.05])*(mx - mn1))
        plt.plot(x, self.galaxy, 'k')
        plt.plot(x[self.goodpixels], resid[self.goodpixels], 'd',
                 color='LimeGreen', mec='LimeGreen', ms=4)
        w = np.flatnonzero(np.diff(self.goodpixels) > 1)
        for wj in w:
            a, b = self.goodpixels[wj : wj + 2]
            plt.axvspan(x[a], x[b], facecolor='lightgray')
            plt.plot(x[a : b + 1], resid[a : b + 1], 'b')
        for k in self.goodpixels[[0, -1]]:
            plt.plot(x[[k, k]], [mn, self.bestfit[k]], 'lightgray')

        if self.gas_any:
            stars_spectrum = self.bestfit - self.gas_bestfit
            plt.plot(x, self.gas_bestfit + mn, c='magenta', linewidth=2)
            plt.plot(x, self.bestfit, c='orange', linewidth=2)
            plt.plot(x, stars_spectrum, 'r', linewidth=2)
        else:
            plt.plot(x, self.bestfit, 'r', linewidth=2)
            plt.plot(x[self.goodpixels], self.goodpixels*0 + mn, '.k', ms=1)

################################################################################

    def format_output(self, params, perror):
        """
        Store the best fitting parameters in the output solution
        and print the results on the console if quiet=False

        """
        p = 0
        self.sol = []
        self.error = []
        for mom in self.moments:
            params[p : p + 2] *= self.velscale  # Bring velocity scale back to km/s
            self.sol.append(params[p : p + mom])
            perror[p : p + 2] *= self.velscale  # Bring velocity scale back to km/s
            self.error.append(perror[p : p + mom])
            p += mom
        if self.mdegree > 0:
            self.mpolyweights = params[p : p + self.mdegree*self.nspec]
        elif self.reddening is not None:
            self.reddening = params[p]  # Replace input with best fit
            self.mpolyweights = None
        if self.gas_reddening is not None:
            self.gas_reddening = params[-1]  # Replace input with best fit
        npoly = (self.degree + 1)*self.nspec
        if self.degree >= 0:
            # output weights for the additive polynomials
            self.polyweights = self.weights[:npoly]
            self.apoly = self.matrix[:, :npoly] @ self.polyweights
        else:
            self.polyweights = self.apoly = None
        # output weights for the templates (and sky) only
        self.weights = self.weights[npoly:]

        if not self.quiet:
            nmom = np.max(self.moments)
            txt = ["Vel", "sigma"] + [f"h{j}" for j in range(3, nmom + 1)]
            print((" Best Fit:" + "{:>10}"*nmom).format(*txt))
            for j, (sol, mom) in enumerate(zip(self.sol, self.moments)):
                print((" comp. {:2d}:" + "{:10.0f}"*2 + "{:10.3f}"*(mom - 2)).format(j, *sol))
            if self.reddening is not None:
                print(f"Stars Reddening E(B-V): {self.reddening:.3f}")
            if self.gas_reddening is not None:
                print(f"Gas Reddening E(B-V): {self.gas_reddening:.3f}")
            print(f"chi2/DOF: {self.chi2:#.4g}; degree = {self.degree}; mdegree = {self.mdegree}")
            print(f"method = {self.method}; Jac calls: {self.njev}; "
                  f"Func calls: {self.nfev}; Status: {self.status}")
            nw = self.weights.size
            print(f"linear_method = {self.linear_method}; Nonzero Templates (>0.1%): ",
                  f"{np.sum(self.weights > 0.001*self.weights.sum())} / {nw}")
            if self.weights.size <= 20:
                print('Templates weights:')
                print(("{:10.3g}"*self.weights.size).format(*self.weights))
            if self.tied is not None:
                tied = np.concatenate(self.tied)
                w = np.flatnonzero(tied != '')
                print("Tied parameters:")
                for j in w:
                    print(f' p[{j}] = {tied[j]}')

        if self.fraction is not None:
            fracFit = np.sum(self.weights[self.component == 0])\
                    / np.sum(self.weights[self.component < 2])
            if not self.quiet:
                print(f"Weights Fraction w[0]/w[0+1]: {fracFit:#.3g}")

        if self.gas_any:
            #print (self.gas_component.shape, npoly, self.matrix[:, npoly:].shape, self.weights.shape)
            #print (self.sky)
            #Changing code for including sky lines
            #Changelog 14/07/2022 Adarsh Ranjan, KASI, Daejeon, ROK
            gas = self.gas_component
            #spectra = self.matrix[:, npoly+2:]    # Remove polynomials
            if self.sky is not None:
                spectra = self.matrix[:, :-(self.sky.shape[1])][:, npoly:]    # Remove polynomials
            else:
                spectra = self.matrix[:, npoly:]    # Remove polynomials

            integ = abs(spectra[:, gas].sum(0))
            #self.gas_flux = integ*self.weights[2:][gas]
            if self.sky is not None:
                self.gas_flux = integ*self.weights[:-(self.sky.shape[1])][gas]
            else:
                self.gas_flux = integ*self.weights[gas]

            design_matrix = spectra[:, gas]/self.noise[:, None]
            self.gas_flux_error = integ*cov_err(design_matrix)[1]
            #self.gas_bestfit = spectra[:, gas] @ self.weights[:-2][gas]
            if self.sky is not None:
                self.gas_bestfit = spectra[:, gas] @ self.weights[:-(self.sky.shape[1])][gas]
            else:
                self.gas_bestfit = spectra[:, gas] @ self.weights[gas]

            if self.gas_any_zero:
                self.gas_flux[self.gas_zero_template] = np.nan
                self.gas_flux_error[self.gas_zero_template] = np.nan

            if not self.quiet:
                print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                print('gas_component   name       flux       err      V     sig')
                print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
                for j, comp in enumerate(self.component[gas]):
                    print("Comp: %2d  %12s %#10.4g  %#8.2g  %6.0f  %4.0f" %
                          (comp, self.gas_names[j], self.gas_flux[j],
                           self.gas_flux_error[j], *self.sol[comp][:2]))
                print('---------------------------------------------------------')

        if self.ncomp == 1:
            self.sol = self.sol[0]
            self.error = self.error[0]

################################################################################