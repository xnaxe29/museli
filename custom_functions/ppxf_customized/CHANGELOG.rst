Changelog
---------

V7.4.4: MC, Oxford, 10 February 2021
++++++++++++++++++++++++++++++++++++

- ``ppxf``: More robust matrix scaling when using linear equality constraints
  in ``constr_templ`` with ``linear_method='lsq_box'``. Thanks to Shravan Shetty
  (pku.edu.cn) for a detailed report and for testing my fix.

V7.4.3: MC, Oxford, 21 December 2020
++++++++++++++++++++++++++++++++++++

- ``capfit``: New ``linear_method`` keyword to select between ``cvxopt`` or
  ``lsq_lin``, when using linear constraints, for cases where the latter stops.
  The ``cvxopt`` package must be installed when setting that option.
- ``ppxf``: Adapted to use ``capfit`` with ``linear_method='cvxopt'`` when
  enforcing linear constraints on the kinematics with ``constr_kinem``.
- ``ppxf``: Included NOTE in the documentation of ``constr_kinem``.
  All changes above were after detailed reports by Kyle Westfall (ucolick.org).

V7.4.2: MC, Oxford, 9 October 2020
++++++++++++++++++++++++++++++++++

- ``ppxf``: Corrected typo in example in the documentation of ``constr_templ``.
- ``ppxf``: Ckeck that ``constr_templ`` and ``constr_kinem`` are dictionaries.
  Thanks to Davide Bevacqua (unibo.it) for the feedback.

V7.4.1: MC, Oxford, 11 September 2020
+++++++++++++++++++++++++++++++++++++

- ``capfit``: Fixed possible infinite loop in ``lsq_box`` and ``lsq_lin``.
  Thanks to Shravan Shetty (pku.edu.cn) for the detailed report and to both
  him and Kyle Westfall (ucolick.org) for testing the fix.
- ``capfit``: Use Numpy rather than Scipy version of ``linalg.lstsq`` to avoid
  a current Scipy bug in the default criterion for rank deficiency.
- ``capfit``: Renamed ``cond`` keyword to ``rcond`` for consistency with Numpy.
- ``capfit``: Passed ``rcond`` keyword to ``cov_err`` function.
- ``ppxf``: removed ``rcond`` keyword in ``capfit`` call. Use default instead.

V7.4.0: MC, Oxford, 20 August 2020
++++++++++++++++++++++++++++++++++

- ``capfit``: New function ``lsq_lin`` implementing a linear least-squares
  linearly-constrained algorithm supporting rank-deficient matrices and allowing
  for a starting guess.
- ``capfit``: Removed the ``lsqlin`` procedure which is superseded by ``lsq_lin``.
- ``capfit``: Renamed ``lsqbox`` to ``lsq_box`` and revised its interface.
- ``ppxf``: Modified to use the new ``lsq_lin`` and the updated ``lsq_box`` functions.
- ``ppxf``: More examples for the ``constr_templ`` and ``constr_kinem`` keywords.
- Set redshift ``z = 0`` when one uncomments the lines to bring the spectrum to
  the rest-frame in ``ppxf_example_kinematics_sdss.py``. Thanks to
  Vaidehi S. Paliya (desy.de) for pointing out the inconsistency in my example.

V7.3.0: MC, Oxford, 10 July 2020
++++++++++++++++++++++++++++++++

- ``capfit``: New function ``lsqbox`` implementing a fast linear least-squares
  box-constrained (bounds) algorithm which allows for a starting guess.
  This also fixed a major mistake in the current implementation of
  ``scipy.optimize.lsq_linear`` (my fix will be included in Scipy).
- ``ppxf``: The new ``linear_method='lsqbox'`` and ``linear_method='cvxopt'``
  now use an initial guess for the solution, which significantly speeds up the
  kinematic fit with multiple templates. As an example, my procedure
  ``ppxf_example_population_gas_sdss`` is now about 4 times faster with the new
  ``linear_method='lsqbox'`` than with the legacy ``linear_method='nnls'``.
- ``ppxf``: Added support for linear equality constraints on the templates
  ``constr_templ`` and for using the keyword ``fraction`` with both
  ``linear_method='lsqbox'`` and ``linear_method='nnls'``.
- Print ``degree`` and ``mdegree`` with the final results.
- Set ``linear=True`` automatically if the fit has no free non-linear parameters,
  to avoid a program stop. Thanks to Shravan Shetty (pku.edu.cn) for the report.

V7.2.1: MC, Oxford, 12 June 2020
++++++++++++++++++++++++++++++++

- ``capfit``: New input keyword ``cond`` for Jacobian rank tolerance.
- ``capfit``: Use ``bvls`` to solve quadratic subproblem with only ``bounds``.
- ``ppxf``: Set ``cond=1e-7`` in ``capfit`` call, when using linear constraints. 
  The ``capfit`` related changes were due to detailed feedback by Kyle Westfall 
  (ucolick.org), to deal with situations with degenerate Jacobians, like when
  there is no stellar continuum and one uses multiplicative polynomials.
- ``ppxf``: Clarified documentation for ``.gas_zero_template`` and the 
  corresponding warning message, after feedback by Laura Salo (umn.edu).

V7.2.0: MC, Oxford, 4 May 2020
++++++++++++++++++++++++++++++

- Allow for ``linear_method='cvxopt'`` when the optional ``cvxopt`` package
  is installed.

V7.1.0: MC, Oxford, 30 April 2020
+++++++++++++++++++++++++++++++++

- Introduced new ``ppxf`` keyword ``linear_method``, and corresponding changes 
  in the code, to select between the old ('nnls') and the new ('lsqlin') 
  approach to the solution of the linear least-squares subproblem in ``ppxf``. 
  Thanks to Sam Vaughan (sydney.edu.au) for a convincing minimal example 
  illustrating the usefulness of this keyword.

V7.0.1: MC, Oxford, 8 April 2020
++++++++++++++++++++++++++++++++

- Support ``.gas_zero_template`` and ``fraction`` together with other 
  equality constraints.
- Included ``np.pad(...mode='constant')`` for backward compatibility with 
  Numpy 1.16. Thanks to Shravan Shetty (KIAA-PKU) for the suggestion.
- Fix ``rebin()`` not retaining the dimensionality of an input column-vector. 
  This resulted in a program stop with a single gas template and 
  ``velscale_ratio > 1``. Thanks to Zhiyuan Ji (astro.umass.edu) for a clear 
  example reproducing the bug.
- ``capfit``: New keyword ``cond`` for ``lsqlin``.
- ``capfit``: Relaxed assertion for inconsistent inequalities in ``lsqlin`` 
  to avoid false positives. Thanks to Kyle Westfall (UCO Lick) for a detailed 
  bug report.

V7.0.0: MC, Oxford, 10 January 2020
+++++++++++++++++++++++++++++++++++

- ``capfit``: New general linear least-squares optimization function
  ``lsqlin`` which is now used to solve the quadratic subproblem.
- ``capfit``: Allow for linear inequality/equality constraints
  ``A_ineq``, ``b_ineq`` and  ``A_eq``, ``b_eq``.
- ``ppxf``: Use (faster) ``capfit.lsqlin`` for the linear fit.
- ``ppxf``: Use updated ``capfit.capfit`` for the non-linear optimization.
- ``ppxf``: Allow for linear equalities/inequalities for both the template
  weights and the kinematic parameters with the ``constr_templ`` and
  ``constr_kinem`` optional keywords.
- ``ppxf``: New ``set_linear_constraints`` function.
- ``ppxf``: Updated documentation.

V6.7.17: MC, Oxford, 14 November 2019
+++++++++++++++++++++++++++++++++++++

- ``capfit``: Written complete documentation.
- ``capfit``: Improved print formatting.
- ``capfit``: Return ``.message`` attribute.
- ``capfit``: Improved ``xtol`` convergence test.
- ``capfit``: Only accept final move if ``chi2`` decreased.
- ``capfit``: Strictly satisfy bounds during Jacobian computation.

V6.7.16: MC, Oxford, 12 June 2019
+++++++++++++++++++++++++++++++++

- ``capfit``: Use only free parameters for ``xtol`` convergence test.
- ``capfit``: Describe in words convergence status with nonzero ``verbose``.
- ``capfit``: Fixed program stop when ``abs_step`` is undefined.
- ``capfit``: Fixed ignoring optional ``max_nfev``.

V6.7.15: MC, Oxford, 7 February 2019
++++++++++++++++++++++++++++++++++++
- Removed unused ``re`` import.
- Removed Scipy's ``next_fast_len`` usage due to an issue with odd padding size.
  Thanks to Eric Emsellem (ESO) for a clear example illustrating this rare and
  subtle bug.

V6.7.14: MC, Oxford, 27 November 2018
++++++++++++++++++++++++++++++++++++++
- Print used ``tied`` parameters equalities, if any.
- Return ``.ndof`` attribute.
- Do not remove ``fixed`` or ``tied`` parameters from the DOF calculation.
  Thanks to Joanna Woo (Univ. of Victoria) for the correction.
- Replaced ``normalize``, ``min_age``, ``max_age`` and ``metal`` keywords with
  ``norm_range``, ``age_range`` and ``metal_range`` in ``ppxf.miles_util.miles``.
- Fixed ``clock`` ``DeprecationWarning`` in Python 3.7.

V6.7.13: MC, Oxford, 20 September 2018
++++++++++++++++++++++++++++++++++++++
- Expanded documentation of ``reddening`` and ``gas_reddening``.
  Thanks to Nick Boardman (Univ. Utah) for the feedback.
- ``capfit`` now raises an error if one tries to tie parameters to themselves.
  Thanks to Kyle Westfall (Univ. Santa Cruz) for the suggestion.
- ``capfit`` uses Python 3.6 f-strings.

V6.7.12: MC, Oxford, 9 July 2018
++++++++++++++++++++++++++++++++
- Allow for ``velscale`` and ``vsyst`` to be Numpy arrays rather than scalars.
- Improved criterion for when the Balmer series is within the fitted wavelength
  range in ``ppxf.ppxf_util.emission_lines``. Thanks to Sam Vaughan
  (Univ. of Oxford) for the feedback.
- Included ``width`` keyword in ``ppxf.ppxf_util.determine_goodpixels``.
  Thanks to George Privon (Univ. of Florida) for the suggestion.
- Expanded ``.gas_flux`` documentation.

V6.7.11: MC, Oxford, 5 June 2018
++++++++++++++++++++++++++++++++

- Formatted ``ppxf.py`` docstring in reStructuredText.
- Removed CHANGELOG from the code and placed in a separate file.
- Modified ``setup.py`` to show help and CHANGELOG on PyPi page.
- Included ``ppxf.__version__``.

V6.7.8: MC, Oxford, 21 May 2018
+++++++++++++++++++++++++++++++

- Moved package to the Python Package Index (PyPi).
- Dropped legacy Python 2.7 support.

V6.7.6: MC, Oxford, 16 April 2018
+++++++++++++++++++++++++++++++++

- Changed imports for the conversion of ppxf to a package.
  Thanks to Joe Burchett (Santa Cruz) for the suggestion.

V6.7.5: MC, Oxford, 10 April 2018
+++++++++++++++++++++++++++++++++

- Fixed syntax error under Python 2.7.

V6.7.4: MC, Oxford, 16 February 2018
++++++++++++++++++++++++++++++++++++

- Fixed bug in ``reddening_cal00()``. It only affected NIR lam > 1000 nm.

V6.7.3: MC, Oxford, 8 February 2018
+++++++++++++++++++++++++++++++++++

- Plot wavelength in nm instead of Angstrom, following IAU rules.
- Ensures each element of ``start`` is not longer than its ``moments``.
- Removed underscore from internal function names.
- Included ``ftol`` keyword.

V6.7.2: MC, Oxford, 30 January 2018
+++++++++++++++++++++++++++++++++++

- Included dunder names as suggested by Peter Weilbacher (Potsdam).
- Fixed wrong ``.gas_reddening`` when ``mdegree > 0``.
- Improved formatting of documentation.

V6.7.1: MC, Oxford, 29 November 2017
++++++++++++++++++++++++++++++++++++

- Removed import of ``misc.factorial``, deprecated in Scipy 1.0.

V6.7.0: MC, Oxford, 6 November 2017
+++++++++++++++++++++++++++++++++++

- Allow users to input identically-zero gas templates while still
  producing a stable NNLS solution. In this case, warn the user and set
  the .gas_zero_template attribute. This situation can indicate an input
  bug or a gas line which entirely falls within a masked region.
- Corrected ``gas_flux_error`` normalization, when input not normalized.
- Return ``.gas_bestfit``, ``.gas_mpoly``, ``.mpoly`` and ``.apoly`` attributes.
- Do not multiply gas emission lines by polynomials, instead allow for
  ``gas_reddening`` (useful with tied Balmer emission lines).
- Use ``axvspan`` to visualize masked regions in plot.
- Fixed program stop with ``linear`` keyword.
- Introduced ``reddening_func`` keyword.

V6.6.4: MC, Oxford, 5 October 2017
++++++++++++++++++++++++++++++++++

- Check for NaN in ``galaxy`` and check all ``bounds`` have two elements.
- Allow ``start`` to be either a list or an array or vectors.

V6.6.3: MC, Oxford, 25 September 2017
+++++++++++++++++++++++++++++++++++++

- Reduced bounds on multiplicative polynomials and clipped to positive
  values. Thanks to Xihan Ji (Tsinghua University) for providing an
  example of slightly negative gas emission lines, when the spectrum
  contains essentially just noise.
- Improved visualization of masked pixels.

V6.6.2: MC, Oxford, 15 September 2017
+++++++++++++++++++++++++++++++++++++

- Fixed program stop with a 2-dim templates array and regularization.
  Thanks to Adriano Poci (Macquarie University) for the clear report and
  the fix.

V6.6.1: MC, Oxford, 4 August 2017
+++++++++++++++++++++++++++++++++

- Included note on ``.gas_flux`` output units. Thanks to Xihan Ji
  (Tsinghua University) for the feedback.

V6.6.0: MC, Oxford, 27 June 2017
++++++++++++++++++++++++++++++++

- Print and return gas fluxes and errors, if requested, with the new
  ``gas_component`` and ``gas_names`` keywords.

V6.5.0: MC, Oxford, 23 June 2017
++++++++++++++++++++++++++++++++

- Replaced ``MPFIT`` with ``capfit`` for a Levenberg-Marquardt method with
  fixed or tied variables, which rigorously accounts for box constraints.

V6.4.2: MC, Oxford, 2 June 2017
+++++++++++++++++++++++++++++++

- Fixed removal of bounds in solution, introduced in V6.4.1.
  Thanks to Kyle Westfall (Univ. Santa Cruz) for reporting this.
- Included ``method`` keyword to use Scipy's ``least_squares()``
  as alternative to MPFIT.
- Force float division in pixel conversion of ``start`` and ``bounds``.

V6.4.1: MC, Oxford, 25 May 2017
+++++++++++++++++++++++++++++++

- ``linear_fit()`` does not return unused status any more, for
  consistency with the correspinding change to ``cap_mpfit``.

V6.4.0: MC, Oxford, 12 May 2017
+++++++++++++++++++++++++++++++

- Introduced ``tied`` keyword to tie parameters during fitting.
- Included discussion of formal errors of ``.weights``.

V6.3.2: MC, Oxford, 4 May 2017
++++++++++++++++++++++++++++++

- Fixed possible program stop introduced in V6.0.7 and consequently
  removed unnecessary function ``_templates_rfft()``. Many thanks to
  Jesus Falcon-Barroso for a very clear and useful bug report!

V6.3.1: MC, Oxford, 13 April 2017
+++++++++++++++++++++++++++++++++

- Fixed program stop when fitting two galaxy spectra with
  reflection-symmetric LOSVD.

V6.3.0: MC, Oxford, 30 March 2017
+++++++++++++++++++++++++++++++++

- Included ``reg_ord`` keyword to allow for both first and second order
  regularization.

V6.2.0: MC, Oxford, 27 March 2017
+++++++++++++++++++++++++++++++++

- Improved curvature criterion for regularization when ``dim > 1``.

V6.1.0: MC, Oxford, 15 March 2017
+++++++++++++++++++++++++++++++++

- Introduced ``trig`` keyword to use a trigonometric series as
  alternative to Legendre polynomials.

V6.0.7: MC, Oxford, 13 March 2017
+++++++++++++++++++++++++++++++++

- Use ``next_fast_len()`` for optimal ``rfft()`` zero padding.
- Included keyword ``gas_component`` in the ``.plot()`` method, to
  distinguish gas emission lines in best-fitting plots.
- Improved plot of residuals for noisy spectra.
- Simplified regularization implementation.

V6.0.6: MC, Oxford, 23 February 2017
++++++++++++++++++++++++++++++++++++

- Added ``linear_fit()`` and ``nonlinear_fit()`` functions to better
  clarify the code structure. Included ``templates_rfft`` keyword.
- Updated documentation. Some code simplifications.

V6.0.5: MC, Oxford, 21 February 2017
++++++++++++++++++++++++++++++++++++

- Consistently use new format_output() function both with/without
  the ``linear`` keyword. Added ``.status`` attribute. Changes suggested by
  Kyle Westfall (Univ. Santa Cruz).

V6.0.4: MC, Oxford, 30 January 2017
+++++++++++++++++++++++++++++++++++

- Re-introduced ``linear`` keyword to only perform a linear fit and
  skip the non-linear optimization.

V6.0.3: MC, Oxford, 1 December 2016
+++++++++++++++++++++++++++++++++++

- Return usual ``Chi**2/DOF`` instead of Biweight estimate.

V6.0.2: MC, Oxford, 15 August 2016
++++++++++++++++++++++++++++++++++

- Improved formatting of printed output.

V6.0.1: MC, Oxford, 10 August 2016
++++++++++++++++++++++++++++++++++

- Allow ``moments`` to be an arbitrary integer.
- Allow for scalar ``moments`` with multiple kinematic components.

V6.0.0: MC, Oxford, 28 July 2016
++++++++++++++++++++++++++++++++

- Compute the Fourier Transform of the LOSVD analytically:
- Major improvement in velocity accuracy when ``sigma < velscale``.
- Removed ``oversample`` keyword, which is now unnecessary.
- Removed limit on velocity shift of templates.
- Simplified FFT zero padding. Updated documentation.

V5.3.3: MC, Oxford 24 May 2016
++++++++++++++++++++++++++++++

- Fixed Python 2 compatibility. Thanks to Masato Onodera (NAOJ).

V5.3.2: MC, Oxford, 22 May 2016
+++++++++++++++++++++++++++++++

- Backward compatibility change: allow ``start`` to be smaller than
  ``moments``. After feedback by Masato Onodera (NAOJ).
- Updated documentation of ``bounds`` and ``fixed``.

V5.3.1: MC, Oxford, 18 May 2016
+++++++++++++++++++++++++++++++

- Use wavelength in plot when available. Make ``plot()`` a class function.
  Changes suggested and provided by Johann Cohen-Tanugi (LUPM).

V5.3.0: MC, Oxford, 9 May 2016
++++++++++++++++++++++++++++++

- Included ``velscale_ratio`` keyword to pass a set of templates with
  higher resolution than the galaxy spectrum.
- Changed ``oversample`` keyword to require integers not Booleans.

V5.2.0: MC, Baltimore, 26 April 2016
++++++++++++++++++++++++++++++++++++

- Included ``bounds``, ``fixed`` and ``fraction`` keywords.

V5.1.18: MC, Oxford, 20 April 2016
++++++++++++++++++++++++++++++++++

- Fixed deprecation warning in Numpy 1.11. Changed order from 1 to 3
  during oversampling. Warn if sigma is under-sampled.

V5.1.17: MC, Oxford, 21 January 2016
++++++++++++++++++++++++++++++++++++

- Expanded explanation of the relation between output velocity and redshift.

V5.1.16: MC, Oxford, 9 November 2015
++++++++++++++++++++++++++++++++++++

- Fixed potentially misleading typo in documentation of ``moments``.

V5.1.15: MC, Oxford, 22 October 2015
++++++++++++++++++++++++++++++++++++

- Updated documentation. Thanks to Peter Weilbacher (Potsdam) for
  corrections.

V5.1.14: MC, Oxford, 19 October 2015
++++++++++++++++++++++++++++++++++++

- Fixed deprecation warning in Numpy 1.10.

V5.1.13: MC, Oxford, 24 April 2015
++++++++++++++++++++++++++++++++++

- Updated documentation.

V5.1.12: MC, Oxford, 25 February 2015
+++++++++++++++++++++++++++++++++++++

- Use ``color=`` instead of ``c=`` to avoid new Matplotlib 1.4 bug.

V5.1.11: MC, Sydney, 5 February 2015
++++++++++++++++++++++++++++++++++++

- Reverted change introduced in V5.1.2. Thanks to Nora Lu"tzgendorf
  for reporting problems with ``oversample``.

V5.1.10: MC, Oxford, 14 October 2014
++++++++++++++++++++++++++++++++++++

- Fixed bug in saving output introduced in previous version.

V5.1.9: MC, Las Vegas Airport, 13 September 2014
++++++++++++++++++++++++++++++++++++++++++++++++

- Pre-compute FFT and oversampling of templates. This speeds up the
  calculation for very long or highly-oversampled spectra. Thanks to
  Remco van den Bosch for reporting situations where this optimization
  may be useful.

V5.1.8: MC, Utah, 10 September 2014
+++++++++++++++++++++++++++++++++++

- Fixed program stop with ``reddening`` keyword. Thanks to Masatao
  Onodera for reporting the problem.

V5.1.7: MC, Oxford, 3 September 2014
++++++++++++++++++++++++++++++++++++

- Relaxed requirement on input maximum velocity shift.
- Minor reorganization of the code structure.

V5.1.6: MC, Oxford, 6 August 2014
+++++++++++++++++++++++++++++++++

- Catch an additional input error. Updated documentation for Python.
  Included templates ``matrix`` in output. Modified plotting colours.

V5.1.5: MC, Oxford, 21 June 2014
++++++++++++++++++++++++++++++++

- Fixed deprecation warning.

V5.1.4: MC, Oxford, 25 May 2014
+++++++++++++++++++++++++++++++

- Support both Python 2.7 and Python 3.

V5.1.3: MC, Oxford, 7 May 2014
++++++++++++++++++++++++++++++

- Allow for an input covariance matrix instead of an error spectrum.

V5.1.2: MC, Oxford, 6 May 2014
++++++++++++++++++++++++++++++

- Replaced REBIN with INTERPOLATE + /OVERSAMPLE keyword. This is
  to account for the fact that the Line Spread Function of the observed
  galaxy spectrum already includes pixel convolution. Thanks to Mike
  Blanton for the suggestion.

V5.1.1: MC, Dallas Airport, 9 February 2014
+++++++++++++++++++++++++++++++++++++++++++

- Fixed typo in the documentation of ``nnls_flags``.

V5.1.0: MC, Oxford, 9 January 2014
++++++++++++++++++++++++++++++++++

- Allow for a different LOSVD for each template. Templates can be
  stellar or can be gas emission lines. A ppxf version adapted for
  multiple kinematic components existed for years. It was updated in
  JAN/2012 for the paper by Johnston et al. (2013, MNRAS). This version
  merges those changes with the public ppxf version, making sure that all
  previous ppxf options are still supported.

V5.0.1: MC, Oxford, 12 December 2013
++++++++++++++++++++++++++++++++++++

- Minor cleaning and corrections.

V5.0.0: MC, Oxford, 6 December 2013
+++++++++++++++++++++++++++++++++++

- Translated from IDL into Python and tested against the original version.

V4.6.6: MC, Paranal, 8 November 2013
++++++++++++++++++++++++++++++++++++

- Uses CAP_RANGE to avoid potential naming conflicts.

V4.6.5: MC, Oxford, 15 November 2012
++++++++++++++++++++++++++++++++++++

- Expanded documentation of REGUL keyword.

V4.6.4: MC, Oxford, 9 December 2011
+++++++++++++++++++++++++++++++++++

- Increased oversampling factor to 30x, when the /OVERSAMPLE keyword
  is used. Updated corresponding documentation. Thanks to Nora
  Lu"tzgendorf for test cases illustrating errors in the recovered
  velocity when the sigma is severely undersampled.

V4.6.3: MC, Oxford 25 October 2011
++++++++++++++++++++++++++++++++++

- Do not change TEMPLATES array in output when REGUL is nonzero.
  From feedback of Richard McDermid.

V4.6.2: MC, Oxford, 17 October 2011
+++++++++++++++++++++++++++++++++++

- Included option for 3D regularization and updated documentation of
  REGUL keyword.

V4.6.1: MC, Oxford, 29 July 2011
++++++++++++++++++++++++++++++++

- Use Coyote Graphics (http://www.idlcoyote.com/) by David W. Fanning.
  The required routines are now included in NASA IDL Astronomy Library.

V4.6.0: MC, Oxford, 12 April 2011
+++++++++++++++++++++++++++++++++

- Important fix to /CLEAN procedure: bad pixels are now properly
  updated during the 3sigma iterations.

V4.5.0: MC, Oxford, 13 April 2010
+++++++++++++++++++++++++++++++++

- Dramatic speed up in the convolution of long spectra.

V4.4.0: MC, Oxford, 18 September 2009
+++++++++++++++++++++++++++++++++++++

- Introduced Calzetti et al. (2000) ppxf_REDDENING_CURVE function to
  estimate the reddening from the fit.

V4.3.0: MC, Oxford, 4 Mach 2009
+++++++++++++++++++++++++++++++

- Introduced REGUL keyword to perform linear regularization of WEIGHTS
  in one or two dimensions.

V4.2.3: MC, Oxford, 27 November 2008
++++++++++++++++++++++++++++++++++++

- Corrected error message for too big velocity shift.

V4.2.2: MC, Windhoek, 3 July 2008
+++++++++++++++++++++++++++++++++

- Added keyword POLYWEIGHTS.

V4.2.1: MC, Oxford, 17 May 2008
+++++++++++++++++++++++++++++++

- Use LA_LEAST_SQUARES (IDL 5.6) instead of SVDC when fitting a single
  template. Please let me know if you need to use ppxf with an older IDL
  version.

V4.2.0: MC, Oxford, 15 March 2008
+++++++++++++++++++++++++++++++++

- Introduced optional fitting of SKY spectrum. Many thanks to
  Anne-Marie Weijmans for testing.

V4.1.7: MC, Oxford, 6 October 2007
++++++++++++++++++++++++++++++++++

- Updated documentation with important note on penalty determination.

V4.1.6: MC, Leiden, 20 January 2006
+++++++++++++++++++++++++++++++++++

- Print number of nonzero templates. Do not print outliers in /QUIET mode.

V4.1.5: MC, Leiden, 10 February 2005
++++++++++++++++++++++++++++++++++++

- Verify that GOODPIXELS is monotonic and does not contain duplicated
  values. After feedback from Richard McDermid.

V4.1.4: MC, Leiden, 12 January 2005
+++++++++++++++++++++++++++++++++++

- Make sure input NOISE is a positive vector.

V4.1.3: MC, Vicenza, 30 December 2004
+++++++++++++++++++++++++++++++++++++

- Updated documentation.

V4.1.2: MC, Leiden, 11 November 2004
++++++++++++++++++++++++++++++++++++

- Handle special case where a single template without additive
  polynomials is fitted to the galaxy.

V4.1.1: MC, Leiden, 21 September 2004
+++++++++++++++++++++++++++++++++++++

- Increased maximum number of iterations ITMAX in BVLS. Thanks to
  Jesus Falcon-Barroso for reporting problems.
- Introduced error message when velocity shift is too big.
- Corrected output when MOMENTS=0.

V4.1.0: MC, Leiden, 3 September 2004
++++++++++++++++++++++++++++++++++++

- Corrected implementation of two-sided fitting of the LOSVD. Thanks
  to Stefan van Dongen for reporting problems.

V4.0.0: MC, Vicenza, 16 August 2004
+++++++++++++++++++++++++++++++++++

- Introduced optional two-sided fitting assuming a reflection
  symmetric LOSVD for two input spectra.

V3.7.3: MC, Leiden, 7 August 2004
+++++++++++++++++++++++++++++++++

- Corrected bug: keyword ERROR was returned in pixels instead of km/s.
- Decreased lower limit on fitted dispersion. Thanks to Igor V. Chilingarian.

V3.7.2: MC, Leiden, 28 April 2004
+++++++++++++++++++++++++++++++++

- Corrected program stop after fit when MOMENTS=2. Bug was introduced in V3.7.0.

V3.7.1: MC, Leiden, 31 March 2004
+++++++++++++++++++++++++++++++++

- Updated documentation.

V3.7.0: MC, Leiden, 23 March 2004
+++++++++++++++++++++++++++++++++

- Revised implementation of MDEGREE option. Nonlinear implementation:
  straightforward, robust, but slower.

V3.6.0: MC, Leiden, 19 March 2004
+++++++++++++++++++++++++++++++++

- Added MDEGREE option for multiplicative polynomials. Linear implementation:
  fast, works well in most cases, but can fail in certain cases.

V3.5.0: MC, Leiden, 11 December 2003
++++++++++++++++++++++++++++++++++++

- Included /OVERSAMPLE option.

V3.4.7: MC, Leiden, 8 December 2003
+++++++++++++++++++++++++++++++++++

- First released version.

V1.0.0: Leiden, 10 October 2001
+++++++++++++++++++++++++++++++

- Created by Michele Cappellari.

