import numpy as np
import sys
import os
from os import path
import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
ppxf_dir = path.dirname(path.realpath(ppxf_package.__file__))
import numpy
#from get_data_from_file import *
import get_data_from_file as gdff
import basic_functions as bf
from astropy.constants import c as c_ms
c_kms = float(c_ms.value / 1e3)
quite_val = False


##################################PPXF_FITTING_FUNCTIONS##################################
##################################PPXF_MAIN_FUNCTION##################################


def ppxf_function(data_file, **kwargs):
	file1 = data_file
	file_type = kwargs.get('data_file_type', 'direct')  # Is this an SDSS data file?
	z = kwargs.get('redshift', 0.0)  # Set the redshift
	require_air_to_vaccum = kwargs.get('air_to_vac', True)  # Does the observed wavelength need air to vacuum correction?
	noise = kwargs.get('flux_err', False)  # Set the default flux uncertainty
	pathname_default = ppxf_dir + '/miles_models/*.fits'  # Initialize a default stellar template
	pathname = kwargs.get('template_address', pathname_default)  # Set the address for Stellar template
	regul_err = kwargs.get('regul_err_value', 0.013)  # Initial value for regularization error
	tie_balmer_val = kwargs.get('tie_balmer', False)  # To Fit with tied Balmer lines
	limit_doublets = kwargs.get('limit_doublets', False)  # To Fit with limited [SII] doublet
	include_reddening = kwargs.get('include_reddening', False)  # To include reddening in the fit
	reddening_val = kwargs.get('reddening_val', 0.0)  # Initial value for reddening
	gas_reddening_val = kwargs.get('gas_reddening_val', 0.0)  # Initial value for gas reddening
	extra_redshift = kwargs.get('extra_redshift', 0.0) # Required for high redshift objects
	kin_lim_1 = kwargs.get('first_order_kinematic_limits', 10000.0) # Required for setting limits for sigma
	kin_lim_2 = kwargs.get('second_order_kinematic_limits', 0.3) # Required for setting limits for h3 and h4
	redshift_lim_lower = kwargs.get('lower_limit_for_redshift', 0.0) # Lower limit for redshift (for limits on V)
	redshift_lim_upper = kwargs.get('upper_limit_for_redshift', 1.0) # Upper limit for redshift (for limits on V)
	adeg_val = kwargs.get('adeg', 4) # Required for setting additive polynomial
	mdeg_val = kwargs.get('mdeg', 0) # Required for setting multiplicative polynomial
	method_of_fit = kwargs.get('fitting_method', 'capfit') # Required for setting method of ppxf fit
	method_of_linear_fit = kwargs.get('linear_fitting_method', 'lsq_box') # Required for setting method of ppxf fit
	sky_spectrum = kwargs.get('include_sky', None) # Required for including sky spectrum for the fit
	sky_flux_init = kwargs.get('sky_flux_init', None) # Required for including sky spectrum for the fit
	number_of_stellar_components = kwargs.get('number_of_stellar_components', 1) # Number of stellar components to fit
	number_of_gas_components = kwargs.get('number_of_gas_components', 0) # Number of gas components to fit
	fix_certain_component = kwargs.get('fix_certain_component', None) # Fix values for a certain component
	fixed_stellar_par = kwargs.get('fixed_stellar_par', [False, False, False, False]) # Fix a specific stellar parameter
	fixed_gas_par = kwargs.get('fixed_gas_par', [False, False]) # Fix a specific gas parameter
	bias_val = kwargs.get('bias_val', None) # Bias
	clean_val = kwargs.get('clean_val', None) # Clean
	fraction_val = kwargs.get('fraction_val', None) # Fraction
	ftol_val = kwargs.get('ftol_val', 1e-4) # Ftol
	good_pixels_array = kwargs.get('good_pixels_array', None) # Goodpixels
	linear_fit_bool = kwargs.get('linear_fit_bool', False) # Linear Fit
	mask_array = kwargs.get('mask_array', None) # Mask
	plot_bool = kwargs.get('plot_bool', False) # Plot?
	reddening_func_val = kwargs.get('reddening_func_val', None) # reddening_func
	sigma_diff_val = kwargs.get('sigma_diff_val', 0) # sigma_diff
	templates_rfft_val = kwargs.get('templates_rfft_val', None) # templates_rfft
	tied_val = kwargs.get('tied_val', None) # tied?
	trig_bool = kwargs.get('trig_bool', False) # trig?
	velscale_ratio_val = kwargs.get('velscale_ratio_val', 1) # Velscale
	x0_val = kwargs.get('x0_val', None) # x0
	quiet_val = kwargs.get('quiet', False) # Quiet run?
	
	include_reddening = include_reddening if not mdeg_val else False
	wave, galaxy, noise, FWHM_gal, velscale, normalising_factor = gdff.get_data_from_files(file1, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, fit_type='ppxf')
	if sky_spectrum is not None:
		sky_wave, sky_flux, sky_noise, FWHM_gal_sky, velscale_sky, normalising_factor_sky = gdff.get_data_from_files(sky_spectrum, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, norm_fact=normalising_factor, fit_type='ppxf')
		sky_val = np.array(np.transpose([np.arange(0, len(sky_wave), 1), sky_flux]))
	elif sky_flux_init is not None:
		wave_data_init, flux_data_init, flux_err_data_init = file1
		sky_flux_err_init = np.ones_like(wave_data_init)
		sky_spectrum = np.array([wave_data_init, sky_flux_init, sky_flux_err_init])
		sky_wave, sky_flux, sky_noise, FWHM_gal_sky, velscale_sky, normalising_factor_sky = gdff.get_data_from_files(sky_spectrum, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, norm_fact=normalising_factor, fit_type='ppxf')
		sky_val = np.array(np.transpose([np.arange(0, len(sky_wave), 1), sky_flux]))
	else:
		sky_flux = None
		sky_val = None
	
	miles = lib.miles(pathname, velscale, FWHM_gal)
	reg_dim = miles.templates.shape[1:]
	lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1. + z)
	dv = c_kms*(miles.ln_lam_temp[0] - np.log(wave[0]))  # eq.(8) of Cappellari (2017)
	vel = c_kms*np.log(1. + z)   # eq.(8) of Cappellari (2017)
	start_star_array = kwargs.get('start_star_array', [vel, 180., 0., 0.]) # initiating starting guess star for [V, sigma, h3, h4]
	start_star = start_star_array     # starting guess star for [V, sigma(--in km/s), h3, h4]
	fixed_array_stars = fixed_stellar_par
	stars_templates = miles.templates.reshape(miles.templates.shape[0], -1)
	V1_lo = c_kms*np.log(1. + redshift_lim_lower)   # eq.(8) of Cappellari (2017)
	V1_up = c_kms*np.log(1. + redshift_lim_upper)   # eq.(8) of Cappellari (2017)
	sigma1_lo = FWHM_gal/10.
	sigma1_up = 200.+(kin_lim_1/2.)
	h3_lo = 0.0-kin_lim_2
	h3_up = 0.0+kin_lim_2
	h4_lo = 0.0-kin_lim_2
	h4_up = 0.0+kin_lim_2
	bounds_stars = [[V1_lo, V1_up], [sigma1_lo, sigma1_up], [h3_lo, h3_up], [h4_lo, h4_up]]
	reddening_val = reddening_val if include_reddening else None
	bf.blockPrint()
	gas_templates, gas_names, line_wave = util.emission_lines(miles.ln_lam_temp, lam_range_gal, FWHM_gal, tie_balmer=tie_balmer_val, limit_doublets=limit_doublets)
	bf.enablePrint()
	n_temps = stars_templates.shape[1]
	start_gas_array = kwargs.get('start_gas_array', [vel, 180.]) # (km/s), initiating starting guess gas for [V, sigma]
	start_gas = start_gas_array     # starting guess gas for [V, sigma(--in km/s)]
	fixed_array_gas = fixed_gas_par
	V2_lo = c_kms*np.log(1. + redshift_lim_lower)   # eq.(8) of Cappellari (2017)
	V2_up = c_kms*np.log(1. + redshift_lim_upper)   # eq.(8) of Cappellari (2017)
	sigma2_lo = FWHM_gal/10.
	sigma2_up = 100.+(kin_lim_1/2.)
	bounds_gas = [[V2_lo, V2_up], [sigma2_lo, sigma2_up]]
	gas_reddening_rev = gas_reddening_val if tie_balmer_val else None
	templates_rev, component_rev, gas_component_rev, gas_names_rev, moments_rev, start_rev, bounds_val, fixed_array = get_params_for_ppxf(stars_templates, bounds_stars, fixed_array_stars, gas_templates, start_star, start_gas, bounds_gas, fixed_array_gas, gas_names, number_of_stellar_components=number_of_stellar_components, number_of_gas_components=number_of_gas_components, quiet=quiet_val, tie_balmer=tie_balmer_val, limit_doublets=limit_doublets)
	if not (fix_certain_component==None):
		for i in range(len(fix_certain_component)):
			moments_rev[fix_certain_component[i]-1] = -moments_rev[fix_certain_component[i]-1]
			
	noise = bf.clean_data(noise, type_of_data='err')
	#bf.print_cust(sky_val.shape)
	#bf.print_cust(galaxy.shape)
	if (number_of_stellar_components) and (number_of_gas_components):
		if not quiet_val: bf.print_cust(f"\nFitting {number_of_stellar_components} stellar absorption with {number_of_gas_components} gas emission components\n", quiet_val=quiet_val)
		#sol_pp = ppxf(templates_rev, galaxy, noise, velscale, start_rev, bias=bias_val, bounds=bounds_val, clean=clean_val, component=component_rev, constr_templ={}, constr_kinem={}, degree=adeg_val, fixed=fixed_array, fraction=fraction_val, ftol=ftol_val, gas_component=gas_component_rev, gas_names=gas_names_rev, gas_reddening=gas_reddening_rev, goodpixels=good_pixels_array, lam=wave, lam_temp=miles.lam_temp, linear=linear_fit_bool, linear_method=method_of_linear_fit, mask=mask_array, method=method_of_fit, mdegree=mdeg_val, moments=moments_rev, plot=plot_bool, quiet=quiet_val, reddening=reddening_val, reddening_func=reddening_func_val, sigma_diff=sigma_diff_val, sky=sky_val, templates_rfft=templates_rfft_val, tied=tied_val, trig=trig_bool, velscale_ratio=velscale_ratio_val, x0=x0_val, reg_dim=reg_dim)
		sol_pp = ppxf(templates_rev, galaxy, noise, velscale, start_rev, bias=bias_val, bounds=bounds_val, clean=clean_val, component=component_rev, constr_templ={}, constr_kinem={}, degree=adeg_val, fixed=fixed_array, fraction=fraction_val, ftol=ftol_val, gas_component=gas_component_rev, gas_names=gas_names_rev, gas_reddening=gas_reddening_rev, goodpixels=good_pixels_array, lam=wave, vsyst=dv, linear=linear_fit_bool, linear_method=method_of_linear_fit, mask=mask_array, method=method_of_fit, mdegree=mdeg_val, moments=moments_rev, plot=plot_bool, quiet=quiet_val, reddening=reddening_val, reddening_func=reddening_func_val, sigma_diff=sigma_diff_val, sky=sky_val, templates_rfft=templates_rfft_val, tied=tied_val, trig=trig_bool, velscale_ratio=velscale_ratio_val, x0=x0_val, reg_dim=reg_dim)
		
	elif not (number_of_stellar_components) and (number_of_gas_components):
		if not quiet_val: bf.print_cust(f"\nFitting {number_of_gas_components} gas emission components\n", quiet_val=quiet_val)
		sol_pp = ppxf(templates_rev, galaxy, noise, velscale, start_rev, bias=bias_val, bounds=bounds_val, clean=clean_val, component=component_rev, constr_templ={}, constr_kinem={}, degree=adeg_val+10, fixed=fixed_array, fraction=fraction_val, ftol=ftol_val, gas_component=gas_component_rev, gas_names=gas_names_rev, gas_reddening=gas_reddening_rev, goodpixels=good_pixels_array, lam=wave, vsyst=dv, linear=linear_fit_bool, linear_method=method_of_linear_fit, mask=mask_array, method=method_of_fit, mdegree=mdeg_val, moments=moments_rev, plot=plot_bool, quiet=quiet_val, reddening=reddening_val, reddening_func=reddening_func_val, sigma_diff=sigma_diff_val, sky=sky_val, templates_rfft=templates_rfft_val, tied=tied_val, trig=trig_bool, velscale_ratio=velscale_ratio_val, x0=x0_val, reg_dim=reg_dim)
	elif (number_of_stellar_components) and not (number_of_gas_components):
		if not quiet_val: bf.print_cust(f"\nFitting {number_of_stellar_components} stellar absorption components\n", quiet_val=quiet_val)
		sol_pp = ppxf(templates_rev, galaxy, noise, velscale, start_rev, bias=bias_val, bounds=bounds_val, clean=clean_val, component=component_rev, constr_templ={}, constr_kinem={}, degree=adeg_val, fixed=fixed_array, fraction=fraction_val, ftol=ftol_val, gas_component=gas_component_rev, gas_names=gas_names_rev, goodpixels=good_pixels_array, lam=wave, vsyst=dv, linear=linear_fit_bool, linear_method=method_of_linear_fit, mask=mask_array, method=method_of_fit, mdegree=mdeg_val, moments=moments_rev, plot=plot_bool, quiet=quiet_val, reddening=reddening_val, reddening_func=reddening_func_val, sigma_diff=sigma_diff_val, sky=sky_val, templates_rfft=templates_rfft_val, tied=tied_val, trig=trig_bool, velscale_ratio=velscale_ratio_val, x0=x0_val, reg_dim=reg_dim)
	else:
		bf.print_cust("\nCannot run ppxf. Select atleast one stellar or gas profile. Exiting.\n", quiet_val=quiet_val)
		quit()

	galaxy = galaxy*normalising_factor
	noise = noise*normalising_factor
	
	return (wave, galaxy, noise, sol_pp, miles, component_rev, gas_names, gas_names_rev, normalising_factor)
##################################PPXF_MAIN_FUNCTION##################################

'''
def get_data_from_file_ppxf(file1, file_type, require_air_to_vaccum, extra_redshift, quiet=True, **kwargs):
	if file_type=='SDSS':
		hdu = fits.open(file1)
		t = hdu[1].data
		wave_orig = 10**(t['loglam'])
		flux_orig = t['flux']
		flux_err_orig = np.full_like(flux_orig, 0.01635)
		wave_ref, flux_ref, flux_err_ref, fwhm_gal_init, velscale = refine_obs_data_using_scipy(wave_orig, flux_orig, flux_err_orig)
	elif file_type=='other':
		wave_orig, flux_orig, flux_err_orig = np.loadtxt(file1, unpack=True, comments='#')
		wave_ref, flux_ref, flux_err_ref, fwhm_gal_init, velscale = refine_obs_data_using_scipy(wave_orig, flux_orig, flux_err_orig)
	elif file_type=='direct':
		wave_orig, flux_orig, flux_err_orig = file1
		wave_ref, flux_ref, flux_err_ref, fwhm_gal_init, velscale = refine_obs_data_using_scipy(wave_orig, flux_orig, flux_err_orig)
		
	if (require_air_to_vaccum):
		#Air to Vacuum conversion
		wave_ref *= np.median(util.vac_to_air(wave_ref)/wave_ref)
		
	FWHM_gal = np.nanmean(fwhm_gal_init)
	bf.print_cust(f'Velscale: {velscale}, FWHM Gal: {FWHM_gal}', quiet_val=quiet)
	fit_type = kwargs.get('fit_type', 'ppxf') # fit_type
	if ('ppxf' in fit_type):
		mask = (wave_ref > 3540.) & (wave_ref < 7409.)
	else:
		mask = np.full([len(wave_ref)], fill_value=True, dtype=np.bool)

	flux = flux_ref[mask]
	normalising_factor = kwargs.get('norm_fact', np.nanmedian(flux)) # Required for setting additive polynomial
	galaxy = flux/normalising_factor   # Normalize spectrum to avoid numerical issues
	wave = wave_ref[mask]
	noise = flux_err_ref[mask]/normalising_factor
	wave = wave / (1.+float(extra_redshift))
	FWHM_gal = FWHM_gal / (1.+float(extra_redshift))
	return (wave, galaxy, noise, FWHM_gal, velscale, normalising_factor)
'''


##################################GET_DETAILS_PPXF_SOLUTION##################################
def get_detailed_ppxf_solution(pp, miles, component_rev, gas_names, gas_names_rev, number_of_stellar_components, number_of_gas_components, **kwargs):
	file1 = kwargs.get('file1', str(ppxf_dir + '/miles_models/Vazdekis2012_ssp_mass_Padova00_UN_baseFe_v10.0.txt')) # x0
	file2_sdss = kwargs.get('file2_sdss', str(ppxf_dir + '/miles_models/Vazdekis2012_ssp_sdss_miuscat_UN1.30_v9.txt')) # x0
	file2_vega = kwargs.get('file2_vega', str(ppxf_dir + '/miles_models/Vazdekis2012_ssp_phot_Padova00_UN_v10.0.txt')) # x0
	band_init = kwargs.get('band', "r") # optical band
	quiet_init = kwargs.get('quiet', False) # quiet run?
	age_limit = kwargs.get('age_limit', 0.001) # Limiting Factor on Age
	metallicity_limit = kwargs.get('metallicity_limit', 0.01) # Limiting Factor on Age
	slope_limit = kwargs.get('slope_limit', 0.01) # Limiting Factor on Age
	normalising_factor_val = kwargs.get('normalising_factor', 1.) # Limiting Factor on Age
	st_age_unique_start_init = [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]
	st_age_unique_start1 = kwargs.get('st_age_unique_start', st_age_unique_start_init) # st_age_unique_array_start
	st_age_unique_end_init = [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]
	st_age_unique_end1 = kwargs.get('st_age_unique_end', st_age_unique_end_init) # st_age_unique_array_end

	bestfit_solution_array = pp.bestfit
	stellar_spectrum_array = (pp.bestfit - pp.gas_bestfit) if (number_of_gas_components) else pp.bestfit

	if (pp.gas_flux) is not None:
		sol_pp_gas_flux_reshaped = pp.gas_flux
		sol_pp_gas_flux_err_reshaped = pp.gas_flux_error
	else:
		sol_pp_gas_flux_reshaped = None
		sol_pp_gas_flux_err_reshaped = None

	str_fit = r'fit ($\rm \chi^{2}$=%2.2f)' %(pp.chi2)
	if (pp.gas_flux is None):
		if (number_of_stellar_components==1):
			str_fit += '\n' + r'Vel(star)=%2.2f' %(pp.sol[0])
			str_fit += '\n' + r'$\rm \sigma$(star)=%2.2f' %(pp.sol[1])
		elif (number_of_stellar_components>1):
			str_fit += '\n' + r'Vel(star)=%2.2f' %(pp.sol[0][0])
			str_fit += '\n' + r'$\rm \sigma$(star)=%2.2f' %(pp.sol[0][1])
	else:
		if (number_of_stellar_components==0):
			if (number_of_gas_components>1):
				str_fit += '\n' + r'Vel(gas)=%2.2f' %(pp.sol[0][0])
				str_fit += '\n' + r'$\rm \sigma$(gas)=%2.2f' %(pp.sol[0][1])
			else:
				str_fit += '\n' + r'Vel(gas)=%2.2f' %(pp.sol[0])
				str_fit += '\n' + r'$\rm \sigma$(gas)=%2.2f' %(pp.sol[1])
		if (number_of_stellar_components>0):
			str_fit += '\n' + r'Vel(star)=%2.2f' %(pp.sol[0][0])
			str_fit += '\n' + r'$\rm \sigma$(star)=%2.2f' %(pp.sol[0][1])
			str_fit += '\n' + r'Vel(gas)=%2.2f' %(pp.sol[-1][0])
			str_fit += '\n' + r'$\rm \sigma$(gas)=%2.2f' %(pp.sol[-1][1])

	if (pp.reddening) is not None:
		str_fit += '\n' + r'E(B-V)=%2.2f' %(pp.reddening)
	if (pp.gas_reddening) is not None:
		str_fit += '\n' + r'E(B-V)(gas)=%2.2f' %(pp.gas_reddening)

	stellar_vel = []
	stellar_sigma = []
	stellar_h3 = []
	stellar_h4 = []
	stellar_vel_err_corrected = []
	stellar_sigma_err_corrected = []
	stellar_h3_err_corrected = []
	stellar_h4_err_corrected = []
	gas_balmer_vel = []
	gas_balmer_sigma = []
	gas_others_vel = []
	gas_others_sigma = []
	gas_balmer_vel_err_corrected = []
	gas_balmer_sigma_err_corrected = []
	gas_others_vel_err_corrected = []
	gas_others_sigma_err_corrected = []
	
	count = 0
	count_gas_balmer = 0
	count_gas_forbidden = 0

	if (number_of_stellar_components) and (number_of_gas_components):
		for i in range(number_of_stellar_components):
			stellar_vel.append(float(pp.sol[i][0]))
			stellar_sigma.append(float(pp.sol[i][1]))
			stellar_h3.append(float(pp.sol[i][2]))
			stellar_h4.append(float(pp.sol[i][3]))
			stellar_vel_err_corrected.append(float(pp.error[i][0]))
			stellar_sigma_err_corrected.append(float(pp.error[i][1]))
			stellar_h3_err_corrected.append(float(pp.error[i][2]))
			stellar_h4_err_corrected.append(float(pp.error[i][3]))
			count+=1
		for j in range(number_of_gas_components):
			gas_balmer_vel.append(float(pp.sol[count+j][0]))
			gas_balmer_sigma.append(float(pp.sol[count+j][1]))
			gas_balmer_vel_err_corrected.append(float(pp.error[count+j][0]))
			gas_balmer_sigma_err_corrected.append(float(pp.error[count+j][1]))
			count_gas_balmer+=1
		for k in range(number_of_gas_components):
			gas_others_vel.append(float(pp.sol[count+count_gas_balmer+k][0]))
			gas_others_sigma.append(float(pp.sol[count+count_gas_balmer+k][1]))
			gas_others_vel_err_corrected.append(float(pp.error[count+count_gas_balmer+k][0]))
			gas_others_sigma_err_corrected.append(float(pp.error[count+count_gas_balmer+k][1]))
			count_gas_forbidden+=1
		
	elif (number_of_stellar_components) and not (number_of_gas_components):
		for i in range(number_of_stellar_components):
			stellar_vel.append(float(pp.sol[i][0]))
			stellar_sigma.append(float(pp.sol[i][1]))
			stellar_h3.append(float(pp.sol[i][2]))
			stellar_h4.append(float(pp.sol[i][3]))
			stellar_vel_err_corrected.append(float(pp.error[i][0]))
			stellar_sigma_err_corrected.append(float(pp.error[i][1]))
			stellar_h3_err_corrected.append(float(pp.error[i][2]))
			stellar_h4_err_corrected.append(float(pp.error[i][3]))
			count+=1

	elif not (number_of_stellar_components) and (number_of_gas_components):
		for j in range(number_of_gas_components):
			gas_balmer_vel.append(float(pp.sol[count+j][0]))
			gas_balmer_sigma.append(float(pp.sol[count+j][1]))
			gas_balmer_vel_err_corrected.append(float(pp.error[count+j][0]))
			gas_balmer_sigma_err_corrected.append(float(pp.error[count+j][1]))
			count_gas_balmer+=1
		for k in range(number_of_gas_components):
			gas_others_vel.append(float(pp.sol[count+count_gas_balmer+k][0]))
			gas_others_sigma.append(float(pp.sol[count+count_gas_balmer+k][1]))
			gas_others_vel_err_corrected.append(float(pp.error[count+count_gas_balmer+k][0]))
			gas_others_sigma_err_corrected.append(float(pp.error[count+count_gas_balmer+k][1]))
			count_gas_forbidden+=1

	status_of_optimization = float(pp.status)
	reddening_fit_e_b_minus_v = pp.reddening
	gas_reddening_fit_e_b_minus_v = pp.gas_reddening
	goodpixels_idx = pp.goodpixels
	reduced_chi_sq = pp.chi2
	reg_dim = miles.templates.shape[1:]

	if (number_of_stellar_components):
		weights = pp.weights
		idx_cust = np.where(component_rev==0)[0]
		weights = weights[idx_cust].reshape(reg_dim)
		st_age_unique, st_mass_unique, st_lum_unique = mass_to_light_customized(miles, weights, file1, file2_sdss, file2_vega, band = band_init, quiet = quiet_init, age_limit = age_limit, metallicity_limit = metallicity_limit, slope_limit = slope_limit, normalising_factor = normalising_factor_val, st_age_unique_start=st_age_unique_start1, st_age_unique_end=st_age_unique_end1)

		weighted_logAge = miles.mean_age_metal(weights, quiet = quiet_init)[0]
		weighted_metallicity = miles.mean_age_metal(weights, quiet = quiet_init)[1]
		mlpop_sol = miles.mass_to_light(weights, band="r", quiet = quiet_init)
		
	else:
		st_age_unique = st_age_unique_start_init
		st_mass_unique = np.zeros_like(st_age_unique_start_init)
		st_lum_unique = np.zeros_like(st_age_unique_start_init)
		weighted_logAge = 0.0
		weighted_metallicity = 0.0
		mlpop_sol = 0.0

	bestfit_solution_array = bestfit_solution_array*normalising_factor_val
	stellar_spectrum_array = stellar_spectrum_array*normalising_factor_val

	return (bestfit_solution_array, stellar_spectrum_array, sol_pp_gas_flux_reshaped, sol_pp_gas_flux_err_reshaped, stellar_vel, stellar_vel_err_corrected, stellar_sigma, stellar_sigma_err_corrected, stellar_h3, stellar_h3_err_corrected, stellar_h4, stellar_h4_err_corrected, gas_balmer_vel, gas_balmer_vel_err_corrected, gas_balmer_sigma, gas_balmer_sigma_err_corrected, gas_others_vel, gas_others_vel_err_corrected, gas_others_sigma, gas_others_sigma_err_corrected, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, goodpixels_idx, reduced_chi_sq, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol, str_fit)
##################################GET_DETAILS_PPXF_SOLUTION##################################



##################################GET_PARAMETERS_FROM_PPXF##################################
def get_params_for_ppxf(stars_templates, bounds_stars, fixed_array_stars, gas_templates, start_star, start_gas, bounds_gas, fixed_array_gas, gas_names, **kwargs):
	number_of_stellar_components = kwargs.get('number_of_stellar_components', 1)  # Number of stellar absorption components
	number_of_gas_components = kwargs.get('number_of_gas_components', 0)  # Number of stellar absorption components
	quiet = kwargs.get('quiet', False)  # Number of stellar absorption components
	tie_balmer = kwargs.get('tie_balmer', False)  # tie_balmer
	limit_doublets = kwargs.get('limit_doublets', False)  # limit_doublets
	forbidden_line_names = []
	for i in range(len(gas_names)):
		if ("[" in gas_names[i]):
			forbidden_line_names.extend([gas_names[i]])
	balmer_line_names = np.setdiff1d(gas_names, forbidden_line_names)
	#n_forbidden = np.sum(["[" in a for a in gas_names])  # forbidden lines contain "[*]"
	#n_balmer = len(gas_names) - n_forbidden
	n_forbidden = len(forbidden_line_names)  # forbidden lines contain "[*]"
	n_balmer = len(balmer_line_names)
	stars_templates /= np.nanmedian(stars_templates)
	gas_templates /= np.nanmedian(gas_templates)
	n_temps = stars_templates.shape[1]
	n_gas_count = len(gas_names)
	stars_templates_rev = np.array([])
	gas_templates_rev = np.array([])
	templates_rev = np.array([])
	component_rev = np.array([], dtype=np.int)
	gas_component_rev = np.array([], dtype=np.bool)
	gas_names_rev = np.array([])
	moments_rev = np.array([], dtype=np.int)
	start_rev = []
	bounds_val = []
	fixed_array = []
	count = 0
	count_balmer = 0
	count_forbidden = 0
	if (number_of_stellar_components) and (number_of_gas_components):
		stars_templates_rev = np.column_stack([stars_templates] * number_of_stellar_components)
		gas_templates_rev = np.column_stack([gas_templates]* number_of_gas_components)
		templates_rev = np.column_stack([stars_templates_rev, gas_templates_rev])
		for i in range(number_of_stellar_components):
			component_rev = np.append(component_rev, [i]*n_temps)
			moments_rev = np.append(moments_rev, len(start_star))
			start_rev.append(start_star)
			bounds_val.append(bounds_stars)
			fixed_array.append(fixed_array_stars)
			count+=1
		for j in range(number_of_gas_components):
			component_rev = np.append(component_rev, [count+j]*n_balmer)
			gas_names_rev = np.append(gas_names_rev, balmer_line_names)
			moments_rev = np.append(moments_rev, len(start_gas))
			start_rev.append(start_gas)
			bounds_val.append(bounds_gas)
			fixed_array.append(fixed_array_gas)
			count_balmer+=1
		for k in range(number_of_gas_components):
			component_rev = np.append(component_rev, [count+count_balmer+k]*n_forbidden)
			gas_names_rev = np.append(gas_names_rev, forbidden_line_names)
			moments_rev = np.append(moments_rev, len(start_gas))
			start_rev.append(start_gas)
			bounds_val.append(bounds_gas)
			fixed_array.append(fixed_array_gas)
			count_forbidden+=1
		gas_component_rev = (count-1 < np.array(component_rev))
		
	elif (number_of_stellar_components) and not (number_of_gas_components):
		stars_templates_rev = np.column_stack([stars_templates] * number_of_stellar_components)
		templates_rev = stars_templates_rev
		for i in range(number_of_stellar_components):
			component_rev = np.append(component_rev, [i]*n_temps)
			moments_rev = np.append(moments_rev, len(start_star))
			start_rev.append(start_star)
			bounds_val.append(bounds_stars)
			fixed_array.append(fixed_array_stars)
			count+=1
		gas_component_rev = None
		gas_names_rev = None
		#start_rev = start_rev[0]
		
	elif not (number_of_stellar_components) and (number_of_gas_components):
		gas_templates_rev = np.column_stack([gas_templates]* number_of_gas_components)
		templates_rev = gas_templates_rev

		for j in range(number_of_gas_components):
			component_rev = np.append(component_rev, [count+j]*n_balmer)
			gas_names_rev = np.append(gas_names_rev, balmer_line_names)
			moments_rev = np.append(moments_rev, len(start_gas))
			start_rev.append(start_gas)
			bounds_val.append(bounds_gas)
			fixed_array.append(fixed_array_gas)
			count_balmer+=1
		for k in range(number_of_gas_components):
			component_rev = np.append(component_rev, [count+count_balmer+k]*n_forbidden)
			gas_names_rev = np.append(gas_names_rev, forbidden_line_names)
			moments_rev = np.append(moments_rev, len(start_gas))
			start_rev.append(start_gas)
			bounds_val.append(bounds_gas)
			fixed_array.append(fixed_array_gas)
			count_forbidden+=1
		gas_component_rev = np.full_like(component_rev, fill_value=True, dtype=np.bool)
		#start_rev = start_rev[0]

	if (len(start_rev)==1):
		start_rev = start_rev[0]
		bounds_val = bounds_val[0]
		fixed_array = fixed_array[0]
	
	if not quiet:
		bf.print_cust(f'\nShape of the template - {templates_rev.shape}', quiet_val=quiet)
		bf.print_cust(f'\nComponents - {component_rev}', quiet_val=quiet)
		bf.print_cust(f'\nGas Components - {gas_component_rev}', quiet_val=quiet)
		bf.print_cust(f'\nGas names - {gas_names_rev}', quiet_val=quiet)
		bf.print_cust(f'\nMoments - {moments_rev}', quiet_val=quiet)
		bf.print_cust(f'\nInitial guess - {start_rev}', quiet_val=quiet)
		bf.print_cust(f'\nBounds - {bounds_val}', quiet_val=quiet)
		bf.print_cust(f'\nFixed? - {fixed_array}', quiet_val=quiet)

	return (templates_rev, component_rev, gas_component_rev, gas_names_rev, moments_rev, start_rev, bounds_val, fixed_array)
##################################GET_PARAMETERS_FROM_PPXF##################################

##################################GET_MASS_TO_LIGHT_INFORMATION##################################
def mass_to_light_customized(miles_cust, weights, file1, file2_sdss, file2_vega, **kwargs):
	assert miles_cust.age_grid.shape == miles_cust.metal_grid.shape == weights.shape, "Input weight dimensions do not match"
	band = kwargs.get('band', "r") # optical band
	quiet = kwargs.get('quiet', False) # quiet run?
	age_limit = kwargs.get('age_limit', 0.001) # Limiting Factor on Age
	metallicity_limit = kwargs.get('metallicity_limit', 0.01) # Limiting Factor on Age
	slope_limit = kwargs.get('slope_limit', 0.01) # Limiting Factor on Age
	normalising_factor_val = kwargs.get('normalising_factor', 1.) # Limiting Factor on Age
	st_age_unique_start_init = [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]
	st_age_unique_start = kwargs.get('st_age_unique_start', st_age_unique_start_init) # st_age_unique_array_start
	st_age_unique_end_init = [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]
	st_age_unique_end = kwargs.get('st_age_unique_end', st_age_unique_end_init) # st_age_unique_array_end

	vega_bands = ["U", "B", "V", "R", "I", "J", "H", "K"]
	sdss_bands = ["u", "g", "r", "i"]
	vega_sun_mag = [5.600, 5.441, 4.820, 4.459, 4.148, 3.711, 3.392, 3.334]
	sdss_sun_mag = [6.55, 5.12, 4.68, 4.57]  # values provided by Elena Ricciardelli
	file_dir = path.dirname(path.realpath(__file__))  # path of this procedure
	if band in vega_bands:
		k = vega_bands.index(band)
		sun_mag = vega_sun_mag[k]
		file2 = file2_vega
	elif band in sdss_bands:
		k = sdss_bands.index(band)
		sun_mag = sdss_sun_mag[k]
		file2 = file2_sdss
	else:
		raise ValueError("Unsupported photometric band")

	slope1, MH1, Age1, m_no_gas = np.loadtxt(file1, usecols=[1, 2, 3, 5]).T
	slope2, MH2, Age2, mag = np.loadtxt(file2, usecols=[1, 2, 3, 4 + k]).T

	# The following loop is a brute force, but very safe and general,
	# way of matching the photometric quantities to the SSP spectra.
	# It makes no assumption on the sorting and dimensions of the files
	mass_no_gas_grid = np.empty_like(weights)
	lum_grid = np.empty_like(weights)
	stellar_age_grid = np.empty_like(weights)
	stellar_mass_grid = np.empty_like(weights)
	metallicity_grid = np.empty_like(weights)

	for j in range(miles_cust.n_ages):
		for k in range(miles_cust.n_metal):
			p1 = (np.abs(miles_cust.age_grid[j, k] - Age1) < age_limit) & (np.abs(miles_cust.metal_grid[j, k] - MH1) < metallicity_limit) & (np.abs(1.30 - slope1) < slope_limit)
			mass_no_gas_grid[j, k] = m_no_gas[p1]
			stellar_mass_grid[j, k] = m_no_gas[p1]

			p2 = (np.abs(miles_cust.age_grid[j, k] - Age2) < age_limit) & (np.abs(miles_cust.metal_grid[j, k] - MH2) < metallicity_limit) & (np.abs(1.30 - slope2) < slope_limit)
			lum_grid[j, k] = (10**(-0.4*(mag[p2] - sun_mag)))
			stellar_age_grid[j, k] = Age2[p2]
			metallicity_grid[j, k] = MH2[p2]

	st_mass = (weights*stellar_mass_grid).flatten()*normalising_factor_val
	st_age = ((stellar_age_grid)).flatten()
	st_lum = (weights*lum_grid).flatten()*normalising_factor_val

	#st_age_unique = np.arange(int(st_age.min()), int(st_age.max()), 1, dtype=np.float)
	#st_age_unique = np.array([0.0, 0.1, 0.5, 1.0, 5.0, 10.0])
	#st_age_unique2 = np.array([0.1, 0.5, 1.0, 5.0, 10.0, 20.0])
	st_age_unique = np.array(st_age_unique_start)
	st_age_unique2 = np.array(st_age_unique_end)
	st_mass_unique = np.zeros([len(st_age_unique)])
	st_lum_unique = np.zeros([len(st_age_unique)])

	for j in range(len(st_age)):
		for i in range(len(st_age_unique)):
			if (st_age_unique[i]<=st_age[j]<=st_age_unique2[i]):
				st_mass_unique[i]+=st_mass[j]
				st_lum_unique[i]+=st_lum[j]

	return (st_age_unique2, st_mass_unique, st_lum_unique)
##################################GET_MASS_TO_LIGHT_INFORMATION##################################

##################################GET_REVISED_RESULT_FROM_PPXF##################################
def ppxf_res_rev(stellar_vel, stellar_sigma, stellar_h3, stellar_h4, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, reduced_chi_sq, gas_balmer_vel, gas_balmer_sigma, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol):
	fit_array = []
	for i in range(len(stellar_vel)):
		fit_array.extend([stellar_vel[i]])
	for i in range(len(stellar_vel)):
		fit_array.extend([stellar_sigma[i]])
	for i in range(len(stellar_vel)):
		fit_array.extend([stellar_h3[i]])
	for i in range(len(stellar_vel)):
		fit_array.extend([stellar_h4[i]])
	fit_array.extend([status_of_optimization])
	fit_array.extend([reddening_fit_e_b_minus_v])
	fit_array.extend([gas_reddening_fit_e_b_minus_v])
	fit_array.extend([reduced_chi_sq])
	for i in range(len(gas_balmer_vel)):
		fit_array.extend([gas_balmer_vel[i]])
	for i in range(len(gas_balmer_vel)):
		fit_array.extend([gas_balmer_sigma[i]])
	for i in range(len(st_age_unique)):
		fit_array.extend([st_age_unique[i]])
	for i in range(len(st_age_unique)):
		fit_array.extend([st_mass_unique[i]])
	for i in range(len(st_age_unique)):
		fit_array.extend([st_lum_unique[i]])
	fit_array.extend([weighted_logAge])
	fit_array.extend([weighted_metallicity])
	fit_array.extend([mlpop_sol])
	return (fit_array)

def ppxf_res_retrieve(fit_array, number_of_stellar_components, number_of_gas_components, quiet=True):
	stellar_vel = []
	count = 0
	for i in range(number_of_stellar_components):
		stellar_vel.extend([fit_array[count]])
		count+=1

	stellar_sigma = []
	for i in range(number_of_stellar_components):
		stellar_sigma.extend([fit_array[count]])
		count+=1

	stellar_h3 = []
	for i in range(number_of_stellar_components):
		stellar_h3.extend([fit_array[count]])
		count+=1

	stellar_h4 = []
	for i in range(number_of_stellar_components):
		stellar_h4.extend([fit_array[count]])
		count+=1

	status_of_optimization = fit_array[count]
	count+=1
	reddening_fit_e_b_minus_v = fit_array[count]
	count+=1
	gas_reddening_fit_e_b_minus_v = fit_array[count]
	count+=1
	reduced_chi_sq = fit_array[count]
	count+=1

	gas_balmer_vel = []
	for i in range(number_of_gas_components):
		gas_balmer_vel.extend([fit_array[count]])
		count+=1

	gas_balmer_sigma = []
	for i in range(number_of_gas_components):
		gas_balmer_sigma.extend([fit_array[count]])
		count+=1

	weighted_logAge = fit_array[-3]
	weighted_metallicity = fit_array[-2]
	mlpop_sol = fit_array[-1]

	stellar_age_pars = fit_array[count:-3]
	index_val = int(len(fit_array[count:-3])/3)
	#bf.print_cust(f'{index_val}', quiet_val=quiet)
	st_age_unique = fit_array[count:count+index_val]
	st_mass_unique = fit_array[count+index_val:count+(2*index_val)]
	st_lum_unique = fit_array[count+(2*index_val):count+(3*index_val)]

	return (stellar_vel, stellar_sigma, stellar_h3, stellar_h4, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, reduced_chi_sq, gas_balmer_vel, gas_balmer_sigma, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol)
##################################GET_REVISED_RESULT_FROM_PPXF##################################

##################################PPXF_FITTING_FUNCTIONS##################################

