import numpy
import numpy as np
from os import path
from scipy.optimize import curve_fit
import pyneb as pn
from astropy.cosmology import Planck15 as LCDM
#from extinction_fitting_functions import *
import extinction_fitting_functions as eff
#from basic_functions import *
import basic_functions as bf
#from ppxf_fitting_functions import *
#import ppxf_fitting_functions as pff
import get_data_from_file as gdff
from astropy.constants import c as c_ms
c_kms = float(c_ms.value / 1e3)
quite_val = False


##################################CUSTOM_EMISSION_FITTING_FUNCTIONS##################################
class fitClass:
	def __init__(self):
		pass
	def gaus_group_with_cont_rev(self, wave_array, *popt, **kwargs):
		group_prof = np.zeros([len(wave_array)])
		quiet_val = kwargs.get('quiet_val', True)  # quiet_val
		#center_list = kwargs.get('center_list', center_list_init)  # Centre List
		#number_of_narrow_components = kwargs.get('number_of_narrow_components', number_of_narrow_components_init)  # Number of narrow components
		#number_of_wide_components = kwargs.get('number_of_wide_components', number_of_wide_components_init)  # Number of wide components
		flux_sky = self.sky_flux
		amp_length = self.length_of_amplitude
		coeff_init = self.coeff_init_val
		position_init_narrow_comp = self.position_init_narrow_comp
		position_init_wide_comp = self.position_init_wide_comp
		position_final_narrow_comp = self.position_final_narrow_comp
		position_final_wide_comp = self.position_final_wide_comp
		comments_on_tied = self.comments_on_tied
		comments_on_balmer = self.comments_on_balmer
		center_list = self.center_list
		number_of_narrow_components = self.number_of_narrow_components
		number_of_wide_components = self.number_of_wide_components
		redshift_val = self.redshift_val
		resolution = self.resolution

		amp_array, center_array, sigma_array, reddening_val, coeff = get_params(popt, number_of_narrow_components, number_of_wide_components, center_list, amp_length)
		if (coeff==None):
			coeff_val = coeff_init
		else:
			coeff_val = coeff
		amp_array_rev = retrieve_all_amplitude_list_rev2(list(amp_array), number_of_narrow_components, number_of_wide_components, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer)
		cont_array = bf.chebyshev_disp(wave_array, coeff_val)
		count_amp = 0
		#bf.print_cust(f'{len(center_list)}, {len(comments_on_balmer)}', quiet_val=quiet_val)
		for j in range(len(center_list)):
			centre = center_list[j]*(1.+redshift_val)
			vel_array = bf.vel_prof(wave_array, centre)
			count = 0
			if (comments_on_balmer[j]):
				for k in range(number_of_narrow_components):
					group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
					count+=1
					count_amp+=1
				for l in range(number_of_wide_components):
					group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
					count+=1
					count_amp+=1
			else:
				for m in range(number_of_narrow_components):
					#bf.print_cust(f'{count}, {count_amp}', quiet_val=quiet_val)
					group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
					count+=1
					count_amp+=1

            #print (count)
        #print (count_amp)
        #quit()
		group_prof+=cont_array
		reddening_array = [reddening_val, 1.0, 0.0, 0.0]
		group_prof_adv = eff.func_6(wave_array/(1.+redshift_val), group_prof, *reddening_array)
		group_prof_adv_rev = group_prof_adv + flux_sky
		group_prof_convolved = bf.convolved_prof5(wave_array, group_prof_adv_rev, resolution)
		return (group_prof_convolved)

def fitting_function_rev(popt_init, datax, datay, yerror, amp_length, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp, coeff_init_val, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, **kwargs):
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	FWHM_gal = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	#Revising FWHM_gal if it is less than 4.0 because it throws error in calculations
	if (FWHM_gal<4.0):
		FWHM_gal = 4.0
	resolution = c_kms/FWHM_gal
	number_of_narrow_components = kwargs.get('number_of_narrow_components', 1)  # number_of_narrow_components
	number_of_wide_components = kwargs.get('number_of_wide_components', 0)  # number_of_wide_components
	center_list = kwargs.get('center_list', [4861.333, 6562.819])  # center_list_init
	comments_on_balmer = kwargs.get('comments_on_balmer', [True, True])  # comments_on_balmer
	comments_on_tied = kwargs.get('comments_on_tied', ['tied_init1', 'tied_paired1'])  # comments_on_tied
	factor_for_tied = kwargs.get('factor_for_tied', [1., 2.86])  # factor_for_tied
	factor_fixed_for_tying = kwargs.get('factor_fixed_for_tying', [True, True])  # factor_fixed_for_tying
	minimal_tying_factor_variance = kwargs.get('minimal_tying_factor_variance', 1e-3)  # minimal_tying_factor_variance
	bound_lower_factor_for_tied = kwargs.get('bound_lower_factor_for_tied', np.array(factor_for_tied)-minimal_tying_factor_variance) #bound_lower_factor_for_tied
	bound_upper_factor_for_tied = kwargs.get('bound_upper_factor_for_tied', np.array(factor_for_tied)+minimal_tying_factor_variance) #bound_upper_factor_for_tied
	minimum_amplitude_val = kwargs.get('minimum_amplitude_val', -1)  # minimum_amplitude_val
	maximum_amplitude_val = kwargs.get('maximum_amplitude_val', 6)  # maximum_amplitude_val
	minimum_reddening_val = kwargs.get('minimum_reddening_val', 1e-4)  # minimum_reddening_val
	maximum_reddening_val = kwargs.get('maximum_reddening_val', 2.0)  # maximum_reddening_val
	fit_velocity_val = kwargs.get('fit_vel', True)  # fit_vel
	fit_vel_disp_val = kwargs.get('fit_sigma', True)  # fit_sigma
	sigma_bound_narrow_min = kwargs.get('sigma_bound_narrow_min', 10.)  # sigma_bound_narrow_min
	sigma_bound_narrow_max = kwargs.get('sigma_bound_narrow_max', 100.0)  # sigma_bound_narrow_max
	sigma_bound_wide_max = kwargs.get('sigma_bound_wide_max', 5000.0)  # sigma_bound_wide_max
	poly_bound_val = kwargs.get('poly_bound_val', 1e4)  # poly_bound_val
	maximum_accepted_reddening = kwargs.get('maximum_accepted_reddening', 2.0)  # maximum_accepted_reddening
	e_b_minus_v_init = kwargs.get('e_b_minus_v_init', 0.0)  # e_b_minus_v_init
	redshift_val = kwargs.get('redshift_val', 0.0)  # redshift_val
	window_for_fit = kwargs.get('window_for_fit', 2000.0)  # window_for_fit
	stopping_number_for_continuum = kwargs.get('stopping_number_for_continuum', 5e5)  # stopping_number_for_continuum
	center_init_array_tmp, sigma_init_array_tmp = bf.get_initial_kinematics_in_kms(int(number_of_narrow_components), int(number_of_wide_components))
	center_init_array = kwargs.get('velocity_init', center_init_array_tmp)  # init_vel
	sigma_init_array = kwargs.get('vel_disp_init', sigma_init_array_tmp)  # init_sigma
	fit_continuum = kwargs.get('fit_continuum', True)  # fit_continuum
	fit_reddening_val = kwargs.get('fit_dust', True)  # fit_dust
	plot_fit_refined = kwargs.get('plot_fit_refined', False)  # plot_fit_refined
	multiplicative_factor_for_plot = kwargs.get('multiplicative_factor_for_plot', 10)  # multiplicative_factor_for_plot
	init_sky_flux = kwargs.get('sky_flux', np.zeros_like(datax))  # sky_flux
	cont_init = kwargs.get('cont_init', np.ones_like(datax))  # cont_init
	method_str = kwargs.get('method_str', 'trf')  # Method of fit
	maxfev_val = kwargs.get('maxfev_val', 2000000)  # Maxfev
	non_coeff_elements = int(amp_length + ((number_of_narrow_components+number_of_wide_components)*2)+1)
	coeff_len = len(popt_init[non_coeff_elements:])
	if (fit_continuum):
		bf.print_cust('Fitting with continuum', quiet_val=quiet_val)
		params_new = popt_init
		bounds_lower, bounds_upper = get_bounds(amp_length, number_of_narrow_components, number_of_wide_components, coeff_len, center_list, factor_fixed_for_tying, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp, bound_lower_factor_for_tied, bound_upper_factor_for_tied, min_amp=minimum_amplitude_val, max_amp=maximum_amplitude_val, center_min_val=-window_for_fit, center_max_val=window_for_fit, fit_continuum=fit_continuum, fit_dust=fit_reddening_val, fit_velocity=fit_velocity_val, fit_vel_disp=fit_vel_disp_val, velocity_init=center_init_array, vel_disp_init=sigma_init_array, sigma_bound_narrow_min=sigma_bound_narrow_min, sigma_bound_narrow_max=sigma_bound_narrow_max, sigma_bound_wide_max=sigma_bound_wide_max, minimum_reddening_val=minimum_reddening_val, maximum_reddening_val=maximum_reddening_val, poly_bound_val_init=poly_bound_val, e_b_minus_v_init=e_b_minus_v_init)
	else:
		params_new = popt_init[0:non_coeff_elements]
		bounds_lower, bounds_upper = get_bounds(amp_length, number_of_narrow_components, number_of_wide_components, coeff_len, center_list, factor_fixed_for_tying, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp, bound_lower_factor_for_tied, bound_upper_factor_for_tied, min_amp=minimum_amplitude_val, max_amp=maximum_amplitude_val, center_min_val=-window_for_fit, center_max_val=window_for_fit, fit_continuum=False, fit_dust=fit_reddening_val, fit_velocity=fit_velocity_val, fit_vel_disp=fit_vel_disp_val, velocity_init=center_init_array, vel_disp_init=sigma_init_array, sigma_bound_narrow_min=sigma_bound_narrow_min, sigma_bound_narrow_max=sigma_bound_narrow_max, sigma_bound_wide_max=sigma_bound_wide_max, minimum_reddening_val=minimum_reddening_val, maximum_reddening_val=maximum_reddening_val, poly_bound_val_init=poly_bound_val, e_b_minus_v_init=e_b_minus_v_init)
		bf.print_cust('Fitting without continuum', quiet_val=quiet_val)

	params_new = params_new.astype(np.float64)
	for i in range(len(params_new)):
		if not (bounds_lower[i]<=params_new[i]<=bounds_upper[i]):
			bf.print_cust(f'Detected out of bound parameter.', quiet_val=quiet_val)
			bf.print_cust(f'{i}, {bounds_lower[i]}, {params_new[i]}, {bounds_upper[i]}', quiet_val=quiet_val)
			bf.print_cust('Fixing it....', quiet_val=quiet_val)
			bounds_lower[i] = params_new[i] - 1e-1
			bounds_upper[i] = params_new[i] + 1e-1
			bf.print_cust(f'{i}, {bounds_lower[i]}, {params_new[i]}, {bounds_upper[i]}', quiet_val=quiet_val)

	inst = fitClass()
	inst.sky_flux = init_sky_flux
	inst.redshift_val = redshift_val
	inst.resolution = resolution
	inst.length_of_amplitude = amp_length
	inst.coeff_init_val = coeff_init_val
	inst.position_init_narrow_comp = position_init_narrow_comp
	inst.position_init_wide_comp = position_init_wide_comp
	inst.position_final_narrow_comp = position_final_narrow_comp
	inst.position_final_wide_comp = position_final_wide_comp
	inst.comments_on_tied = comments_on_tied
	inst.comments_on_balmer = comments_on_balmer
	inst.center_list = center_list
	inst.number_of_narrow_components = number_of_narrow_components
	inst.number_of_wide_components = number_of_wide_components
    #print(params_new)
    #quit()

	try:
		pfit, pcov = curve_fit(inst.gaus_group_with_cont_rev, datax, datay, p0=params_new, bounds=((bounds_lower), (bounds_upper)), sigma=yerror, maxfev=maxfev_val, method=method_str)
	except RuntimeError:
		bf.print_cust('Error - curve_fit failed')
		pfit = np.ones_like(params_new)
		pcov = np.ones_like(params_new)

	error = []
	for i in range(len(pfit)):
		try:
			error.append(np.absolute(pcov[i][i])**0.5)
		except:
			error.append(0.00)
	pfit_curvefit = pfit
	perr_curvefit = np.array(error)
	return pfit_curvefit, perr_curvefit

def get_line_flux(wave, flux, err, cont, fit, sky, center_array, sigma_array, center_list, redshift):
	line_flux_meh = np.zeros([len(center_list)])
	line_flux_meh_err = np.zeros([len(center_list)])
	profile_meh = (fit/cont) - (sky/cont)
	profile_meh_err = err / cont
	wave_rest_meh = wave/(1.+redshift)
	for meh in range(len(center_list)):
		vel_array_meh = bf.vel_prof(wave_rest_meh, center_list[meh])
		idx1_meh = []
		idx2_meh = []
		for weh in range(len(center_array)):
			idx1_meh.extend([np.searchsorted(vel_array_meh, (center_array[weh]-2*sigma_array[weh]))])
			idx2_meh.extend([np.searchsorted(vel_array_meh, (center_array[weh]+2*sigma_array[weh]))])

		idx1_meh = np.nanmin(idx1_meh)
		idx2_meh = np.nanmax(idx2_meh)
		line_flux_meh[meh] = np.nansum(profile_meh[idx1_meh:idx2_meh]) * (wave_rest_meh[idx1_meh+1] - wave_rest_meh[idx1_meh])
		line_flux_meh_err[meh] = np.nansum(profile_meh_err[idx1_meh:idx2_meh]) * (wave_rest_meh[idx1_meh+1] - wave_rest_meh[idx1_meh])

	return(line_flux_meh, line_flux_meh_err)


def complex_fit_func(wave_init, flux_init, err_init, par_dict, **kwargs):
	par_dict_rev = par_dict
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_balmer'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_balmer: {len(par_dict_rev['comments_on_balmer'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_tied: {len(par_dict_rev['comments_on_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_for_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_for_tied: {len(par_dict_rev['factor_for_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_fixed_for_tying'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_fixed_for_tying: {len(par_dict_rev['factor_fixed_for_tying'])}"
	flux_sky_init = kwargs.get('flux_sky', np.zeros_like(wave_init))  # flux_sky
	flux_sky_init_err = kwargs.get('flux_sky_err', np.zeros_like(wave_init))  # flux_sky
	cont_init1 = kwargs.get('cont_init', np.ones_like(wave_init))  # cont_init
	cont_init_err = kwargs.get('cont_init_err', np.ones_like(wave_init))  # cont_init
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	FWHM_gal1 = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	file_type = kwargs.get('data_file_type', 'direct')  # Is this an SDSS data file?
	require_air_to_vaccum = kwargs.get('air_to_vac', True)  # Does the observed wavelength need air to vacuum correction?
	extra_redshift = kwargs.get('extra_redshift', 0.0) # Required for high redshift objects
	#quiet_val = kwargs.get('quiet', False) # Quiet run?

	data_file1 = [wave_init, flux_init, err_init]
	sky_file1 = [wave_init, flux_sky_init, flux_sky_init_err]
	cont_file1 = [wave_init, cont_init1, cont_init_err]
	wave, flux, err, FWHM_gal1, velscale1, normalising_factor1 = gdff.get_data_from_files(data_file1, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, fit_type='custom')
	if (FWHM_gal1<4.0):
		FWHM_gal1 = 4.0
	resolution = float(c_kms / FWHM_gal1)

	wave2, flux_sky, sky_flux_err, FWHM_gal2, velscale2, normalising_factor2 = gdff.get_data_from_files(sky_file1, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, norm_fact=normalising_factor1, fit_type='custom')
	sky_val = np.array(np.transpose([np.arange(0, len(wave2), 1), flux_sky]))
	wave3, cont_init, cont_init_err, FWHM_gal3, velscale3, normalising_factor3 = gdff.get_data_from_files(cont_file1, file_type, require_air_to_vaccum, extra_redshift, quiet=quiet_val, norm_fact=normalising_factor1, fit_type='custom')
	continuum_fitted = cont_init

	popt_init, masked_array, amp_init_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, coeff_init_val, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp = get_opt_param_array_rev2(wave, flux, cont_init, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], par_dict_rev['redshift_val'], par_dict_rev['comments_on_balmer'], par_dict_rev['comments_on_tied'], par_dict_rev['factor_for_tied'], par_dict_rev['e_b_minus_v_init'], par_dict_rev['window_for_fit'], par_dict_rev['stopping_number_for_continuum'], par_dict_rev['center_init_array'], par_dict_rev['sigma_init_array'])
	amp_length = len(amp_init_array)
	continuum = bf.chebyshev_disp(wave[masked_array], coeff_init_val)
	flux_sky_orig = flux_sky[masked_array]

	pfit_curvefit, perr_curvefit = fitting_function_rev(popt_init, wave[masked_array], flux[masked_array], err[masked_array], amp_length, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp, coeff_init_val, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied=par_dict_rev['comments_on_tied'], comments_on_balmer=par_dict_rev['comments_on_balmer'], number_of_narrow_components=par_dict_rev['number_of_narrow_components'], number_of_wide_components=par_dict_rev['number_of_wide_components'], fit_continuum=par_dict_rev['fit_continuum_val'], fit_dust=par_dict_rev['fit_reddening_val'], sky_flux = flux_sky_orig, center_list=par_dict_rev['center_list'], factor_for_tied = par_dict_rev['factor_for_tied'], factor_fixed_for_tying = par_dict_rev['factor_fixed_for_tying'], redshift_val=par_dict_rev['redshift_val'], FWHM_gal=FWHM_gal1, quiet_val=quiet_val)
	amp_array, center_array, sigma_array, reddening_val_fit, coeff_fit = get_params(pfit_curvefit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	if (reddening_val_fit>par_dict_rev['maximum_accepted_reddening']):
		bf.print_cust(f'Reddening val: {reddening_val_fit} exceeded maximum accepted value... Reverting to fixed continuum and reddening fit')
		pfit_curvefit, perr_curvefit = fitting_function_rev(popt_init, wave[masked_array], flux[masked_array], err[masked_array], amp_length, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp, coeff_init_val, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied=par_dict_rev['comments_on_tied'], comments_on_balmer=par_dict_rev['comments_on_balmer'], number_of_narrow_components=par_dict_rev['number_of_narrow_components'], number_of_wide_components=par_dict_rev['number_of_wide_components'], fit_continuum=False, fit_dust=False, center_list=par_dict_rev['center_list'], factor_for_tied = par_dict_rev['factor_for_tied'], factor_fixed_for_tying = par_dict_rev['factor_fixed_for_tying'], redshift_val=par_dict_rev['redshift_val'], FWHM_gal=FWHM_gal1, quiet_val=quiet_val)
	bf.print_cust('Fit complete... Saving results...', quiet_val=quiet_val)
	wave_fitted, res_fitted = plot_gaus_group_with_cont_rev(wave[masked_array], par_dict_rev['plot_fit_refined'], par_dict_rev['multiplicative_factor_for_plot'], par_dict_rev['center_list'], par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], amp_length, coeff_init_val, flux_sky_orig, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, par_dict_rev['comments_on_tied'], par_dict_rev['comments_on_balmer'], par_dict_rev['redshift_val'], resolution, *pfit_curvefit)
	dof = len(flux[masked_array]) - len(pfit_curvefit)
	red_chi_squared = np.sum( ( (flux[masked_array] - res_fitted)**2. / (err[masked_array])**2. )) / dof
	bf.print_cust(f'Reduced chi-sqaure: {red_chi_squared}', quiet_val=quiet_val)
	amp_array, center_array, sigma_array, reddening_val_fit, coeff_fit = get_params(pfit_curvefit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	amp_array_rev = retrieve_all_amplitude_list_rev2(list(amp_array), par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, par_dict_rev['comments_on_tied'], par_dict_rev['comments_on_balmer'])
	amp_err_array, center_err_array, sigma_err_array, reddening_err_val_fit, coeff_err_fit = get_params(perr_curvefit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	if (coeff_fit is None):
		coeff_val_fit = coeff_init_val
		coeff_err_val_fit = np.zeros([len(coeff_init_val)])
	else:
		coeff_val_fit = coeff_fit
		coeff_err_val_fit = coeff_err_fit

	continuum_fitted_partial = bf.chebyshev_disp(wave[masked_array], coeff_val_fit)
	continuum_fitted[masked_array] = continuum_fitted_partial
	res_fitted_full = np.zeros_like(wave)
	res_fitted_full[masked_array] = res_fitted
	pfit_curvefit = np.append(pfit_curvefit, coeff_val_fit)
	perr_curvefit = np.append(perr_curvefit, coeff_err_val_fit)

	amp_array, center_array, sigma_array, reddening_val_fit, coeff_fit = get_params(pfit_curvefit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	amp_err_array, center_err_array, sigma_err_array, reddening_err_val_fit, coeff_err_fit = get_params(perr_curvefit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	line_flux, line_flux_err = get_line_flux(wave, flux, err, continuum_fitted, res_fitted_full, flux_sky, center_array, sigma_array, par_dict_rev['center_list'], par_dict_rev['redshift_val'])

	continuum_fitted = continuum_fitted*normalising_factor1
	res_fitted_full = res_fitted_full*normalising_factor1
	flux = flux*normalising_factor1
	err = err*normalising_factor1

	return (wave, flux, err, continuum_fitted, res_fitted_full, masked_array, pfit_curvefit, perr_curvefit, line_flux, line_flux_err, red_chi_squared, amp_length)


##################################EXTRA_FUNCTION_FOR_CUSTOM_EMISSION_FIT##################################
##################################GET_PARAMS_AND_BOUNDS##################################

def get_params(popt, number_of_narrow_components, number_of_wide_components, center_list, amp_length):
	amplitude_len = amp_length
	center_len = int(amplitude_len+(number_of_narrow_components+number_of_wide_components))
	sigma_len = int(center_len+(number_of_narrow_components+number_of_wide_components))
	amp_array = popt[0:amplitude_len]
	#bf.print_cust(f'{type(amplitude_len)}, {type(center_len)}')
	center_array = popt[amplitude_len:center_len]
	sigma_array = popt[center_len:sigma_len]
	reddening_val = popt[sigma_len]
	if len(popt)==(sigma_len+1):
		coeff=None
	else:
		coeff = popt[sigma_len+1:]

	#print (amp_array, center_array, sigma_array, reddening_val, coeff)
	#quit()
	return(amp_array, center_array, sigma_array, reddening_val, coeff)

def get_bounds(amp_length, number_of_narrow_components, number_of_wide_components, coeff_len, center_list, amp_ratio_fixed, idx_for_fixing_bounds_factor_init_narrow_comp, idx_for_fixing_bounds_factor_init_wide_comp, bound_lower_factor_for_tied_amp, bound_upper_factor_for_tied_amp, **kwargs):
	min_amp = kwargs.get('min_amp', -1)  # Maximum Amplitude
	max_amp = kwargs.get('max_amp', 6)  # Maximum Amplitude
	center_min_val = kwargs.get('center_min_val', -1000.)  # Minimum center val
	center_max_val = kwargs.get('center_max_val', 1000.)  # Maximum center val
	minimum_reddening_val = kwargs.get('minimum_reddening_val', 1e-4)  # Maximum center val
	maximum_reddening_val = kwargs.get('maximum_reddening_val', 2.)  # Maximum center val
	fit_continuum = kwargs.get('fit_continuum', True)  # Fit continuum?
	fit_dust = kwargs.get('fit_dust', True)  # Fit dust?
	fit_vel = kwargs.get('fit_velocity', True)  # Fit Velocity?
	fit_sigma = kwargs.get('fit_vel_disp', True)  # Fit Velocity Dispersion?
	init_vel = kwargs.get('velocity_init', np.full([number_of_narrow_components+number_of_wide_components], fill_value=0.))  # Initial guess for Velocity
	init_sigma = kwargs.get('vel_disp_init', np.full([number_of_narrow_components+number_of_wide_components], fill_value=50.))  # Initial guess for Velocity Dispersion
	sigma_bound_narrow_min = kwargs.get('sigma_bound_narrow_min', 10.)  # Minimal sigma limit for Narrow line
	sigma_bound_narrow_max = kwargs.get('sigma_bound_narrow_max', 100.)  # Maximum sigma limit for Narrow line
	sigma_bound_wide_max = kwargs.get('sigma_bound_wide_max', 5000.)  # Maximum sigma limit for Wide line
	poly_bound_val = kwargs.get('poly_bound_val_init', 1e4)  # Maximum sigma limit for Wide line
	e_b_minus_v_init = kwargs.get('e_b_minus_v_init', 1e-3)  # Maximum sigma limit for Wide line
	
	if (fit_vel==True):
		center_bound_min = np.full([(number_of_narrow_components+number_of_wide_components)], center_min_val)
		center_bound_max = np.full([(number_of_narrow_components+number_of_wide_components)], center_max_val)
	else:
		center_bound_min = np.array(init_vel)-10.
		center_bound_max = np.array(init_vel)+10.
	if (fit_sigma==True):
		sigma_bound_min = np.full([(number_of_narrow_components)], sigma_bound_narrow_min)
		sigma_bound_min = np.append(sigma_bound_min, np.full([(number_of_wide_components)], sigma_bound_narrow_max))
		sigma_bound_max = np.full([(number_of_narrow_components)], sigma_bound_narrow_max)
		sigma_bound_max = np.append(sigma_bound_max, np.full([(number_of_wide_components)], sigma_bound_wide_max))
	else:
		sigma_bound_min = np.array(init_sigma)-10.
		sigma_bound_max = np.array(init_sigma)+10.
	#amp_ratio_fixed = kwargs.get('factor_fixed_for_tying', factor_fixed_for_tying)  # Fit dust?
	#idx_for_fixing_bounds_factor_init_narrow_comp = kwargs.get('index_for_fixing_bounds_factor_init_narrow_comp', index_for_fixing_bounds_factor_init_narrow_comp)  # Fit dust?
	#idx_for_fixing_bounds_factor_init_wide_comp = kwargs.get('index_for_fixing_bounds_factor_init_wide_comp', index_for_fixing_bounds_factor_init_wide_comp)  # Fit dust?
	#bound_lower_factor_for_tied_amp = kwargs.get('bound_lower_factor_for_tied', bound_lower_factor_for_tied)  # bound_lower_factor_for_tied
	#bound_upper_factor_for_tied_amp = kwargs.get('bound_upper_factor_for_tied', bound_upper_factor_for_tied)  # bound_upper_factor_for_tied
	amplitude_len = amp_length
	amp_bound_min = np.full([(amplitude_len)], min_amp)
	amp_bound_max = np.full([(amplitude_len)], max_amp)
	idx_for_fixing_bounds_factor_init_wide_comp = idx_for_fixing_bounds_factor_init_wide_comp[::-1]
	for i in range(len(idx_for_fixing_bounds_factor_init_wide_comp)):
		if (amp_ratio_fixed[int(idx_for_fixing_bounds_factor_init_wide_comp[i])]==True):
			amp_bound_min[-(i+1)] = bound_lower_factor_for_tied_amp[int(idx_for_fixing_bounds_factor_init_wide_comp[i])]
			amp_bound_max[-(i+1)] = bound_upper_factor_for_tied_amp[int(idx_for_fixing_bounds_factor_init_wide_comp[i])]
		else:
			amp_bound_min[-(i+1)] = 1e-5
			amp_bound_max[-(i+1)] = 1e5

	idx_for_fixing_bounds_factor_init_narrow_comp = idx_for_fixing_bounds_factor_init_narrow_comp[::-1]
	for i in range(len(idx_for_fixing_bounds_factor_init_narrow_comp)):
		if (amp_ratio_fixed[int(idx_for_fixing_bounds_factor_init_narrow_comp[i])]==True):
			amp_bound_min[-(i+len(idx_for_fixing_bounds_factor_init_wide_comp)+1)] = bound_lower_factor_for_tied_amp[int(idx_for_fixing_bounds_factor_init_narrow_comp[i])]
			amp_bound_max[-(i+len(idx_for_fixing_bounds_factor_init_wide_comp)+1)] = bound_upper_factor_for_tied_amp[int(idx_for_fixing_bounds_factor_init_narrow_comp[i])]
		else:
			amp_bound_min[-(i+len(idx_for_fixing_bounds_factor_init_wide_comp)+1)] = 1e-5
			amp_bound_max[-(i+len(idx_for_fixing_bounds_factor_init_wide_comp)+1)] = 1e5

	if (fit_dust):
		reddening_bound_min = np.full([(1)], minimum_reddening_val)
		reddening_bound_max = np.full([(1)], maximum_reddening_val)
	else:
		reddening_bound_min = np.full([(1)], (e_b_minus_v_init-(e_b_minus_v_init/10.)))
		reddening_bound_max = np.full([(1)], (e_b_minus_v_init+(e_b_minus_v_init/10.)))

	coeff_lower_bound = np.full([coeff_len], -poly_bound_val)
	coeff_upper_bound = np.full([coeff_len], poly_bound_val)
	bounds_lower = []
	bounds_lower.extend(amp_bound_min)
	bounds_lower.extend(center_bound_min)
	bounds_lower.extend(sigma_bound_min)
	bounds_lower.extend(reddening_bound_min)
	if (fit_continuum):
		bounds_lower.extend(coeff_lower_bound)

	bounds_upper = []
	bounds_upper.extend(amp_bound_max)
	bounds_upper.extend(center_bound_max)
	bounds_upper.extend(sigma_bound_max)
	bounds_upper.extend(reddening_bound_max)
	if (fit_continuum):
		bounds_upper.extend(coeff_upper_bound)

	return(bounds_lower, bounds_upper)

##################################GET_PARAMS_AND_BOUNDS##################################

##################################GET_FITTED_TABLE##################################

def get_table(new_fitted_popt_param_array_rev, new_fitted_perr_param_array, center_list_init_names, number_of_narrow_components_init, number_of_wide_components_init):
	final_array_to_save = np.chararray([(len(center_list_init_names)+2), (number_of_narrow_components_init+number_of_wide_components_init)], itemsize=100)
	final_array_to_save[:] = 'a'
	for i in range(final_array_to_save.shape[0]):
		for j in range(final_array_to_save.shape[1]):
			final_array_to_save[i,j] = str(np.round(new_fitted_popt_param_array_rev[i,j],3)) + str("+/-") + str(np.round(new_fitted_perr_param_array[i,j],3))
	column_text = np.append(np.array(['Comp']), center_list_init_names)
	column_text = np.append(column_text, np.array(['V', 'sigma']))
	row_text_1 = np.arange(1, number_of_narrow_components_init+1)
	row_text_1 = np.append(row_text_1, np.arange(number_of_narrow_components_init+1, number_of_narrow_components_init+1+number_of_wide_components_init))
	row_text_1 = list(row_text_1.astype(np.str))
	string = 'Comp'
	row_text = [string + x for x in row_text_1]
	final_array_to_save = numpy.vstack([row_text, final_array_to_save])
	final_array_to_save = np.append(np.array([column_text]).transpose(), final_array_to_save, axis=1)
	return (final_array_to_save)

def make_table_with_result(wave_array, new_fitted_popt_param_array, center_list, number_of_narrow_components, number_of_wide_components, amp_length, coeff_init, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer, redshift_val, *popt, **kwargs):
	check_data_type = kwargs.get('check_data_type', 'data')  # length of amplitude array to be fitted
	amp_array, center_array, sigma_array, reddening_val, coeff = get_params(popt, number_of_narrow_components, number_of_wide_components, center_list, amp_length)
	if (coeff==None):
		coeff_val = coeff_init
	else:
		coeff_val = coeff
	amp_array_rev = retrieve_all_amplitude_list_rev2(list(amp_array), number_of_narrow_components, number_of_wide_components, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer)
	if (check_data_type=='data'):
		amp_array_rev = np.power(10, amp_array_rev)
	count_amp = 0
	new_fitted_popt_param_array_rev = new_fitted_popt_param_array
	for j in range(len(center_list)):
		centre = center_list[j]*(1.+redshift_val)
		vel_array = bf.vel_prof(wave_array, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components):
				new_fitted_popt_param_array[j, count] = amp_array_rev[count_amp]
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components):
				new_fitted_popt_param_array[j, count] = amp_array_rev[count_amp]
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components):
				new_fitted_popt_param_array[j, count] = amp_array_rev[count_amp]
				count+=1
				count_amp+=1
	return (new_fitted_popt_param_array_rev)


#make_table_with_result(wave_array, new_fitted_popt_param_array, center_list, number_of_narrow_components, number_of_wide_components, amp_length, coeff_init, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer, *popt, check_data_type='data')

##################################GET_FITTED_TABLE##################################

##################################PLOT_RESULTS##################################

def plot_gaus_group_with_cont_rev(wave_array_rev, plot_fit_refined, multiplicative_factor_for_plot, center_list, number_of_narrow_components, number_of_wide_components, amp_length, coeff_init, flux_sky, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer, redshift_val, resolution, *popt):
	if (plot_fit_refined):
		wave_array = np.logspace(np.log10(wave_array_rev.min()), np.log10(wave_array_rev.max()), num=len(wave_array_rev)*multiplicative_factor_for_plot, endpoint=True, base=10)
	else:
		wave_array = wave_array_rev
	group_prof = np.zeros([len(wave_array)])
	amp_array, center_array, sigma_array, reddening_val, coeff = get_params(popt, number_of_narrow_components, number_of_wide_components, center_list, amp_length)
	if (coeff==None):
		coeff_val = coeff_init
	else:
		coeff_val = coeff
	amp_array_rev = retrieve_all_amplitude_list_rev2(list(amp_array), number_of_narrow_components, number_of_wide_components, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer)
	cont_array = bf.chebyshev_disp(wave_array, coeff_val)
	count_amp = 0
	for j in range(len(center_list)):
		centre = center_list[j]*(1.+redshift_val)
		vel_array = bf.vel_prof(wave_array, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components):
				group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components):
				group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components):
				group_prof += bf.gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
				count+=1
				count_amp+=1
	group_prof+=cont_array
	reddening_array = [reddening_val, 1.0, 0.0, 0.0]
	group_prof_adv = eff.func_6(wave_array/(1.+redshift_val), group_prof, *reddening_array)
	group_prof_adv_rev = group_prof_adv + flux_sky
	group_prof_convolved = bf.convolved_prof5(wave_array, group_prof_adv_rev, resolution)
	return (wave_array, group_prof_convolved)

#plot_gaus_group_with_cont_rev(wave_array_rev, plot_fit_refined, multiplicative_factor_for_plot, center_list, number_of_narrow_components, number_of_wide_components, amp_length, coeff_init, flux_sky, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer, redshift_val, resolution, *popt)

##################################PLOT_RESULTS##################################

##################################GET_OPTIMUM_PARAMETERS##################################

def get_opt_param_array_rev2(wave, flux, cont, number_of_narrow_components_init, number_of_wide_components_init, center_list_init, redshift_val, comments_on_balmer, comments_on_tied, factor_for_tied, e_b_minus_v_init, window_for_fit, stopping_number_for_continuum, center_init_array, sigma_init_array):
	popt_init = []
	amplitude_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp = get_fitting_amplitude_list_rev2(wave, flux, number_of_narrow_components_init, number_of_wide_components_init, center_list_init, redshift_val, comments_on_balmer, comments_on_tied, factor_for_tied)
	popt_init.extend(amplitude_array)
	popt_init.extend(center_init_array)
	popt_init.extend(sigma_init_array)
	popt_init.extend([e_b_minus_v_init])
	masked_array = numpy.zeros((len(wave)), dtype=bool)
	for i in range(len(center_list_init)):
		centre = center_list_init[i]*(1.+redshift_val)
		vel_array = bf.vel_prof(wave, centre)
		mask = (vel_array > -2*window_for_fit) & (vel_array < 2*window_for_fit)
		masked_array[mask] = True
	order_test_init, coeff_init_val = bf.chebyshev_order(wave[masked_array], cont[masked_array], stopping_number=stopping_number_for_continuum)
	popt_init.extend(coeff_init_val)
	popt_init = np.nan_to_num(popt_init, copy=True, nan=1e-6, posinf=1e-6, neginf=1e-6)
	return (popt_init, masked_array, amplitude_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, coeff_init_val, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp)


def get_fitting_amplitude_list_rev2(wave, flux, number_of_narrow_components_init, number_of_wide_components_init, center_list_init, redshift_val, comments_on_balmer, comments_on_tied, factor_for_tied):
	amplitude_array = []
	position_init_narrow_comp = []
	position_init_wide_comp = []
	factor_init_narrow_comp = []
	factor_init_wide_comp = []
	index_for_fixing_bounds_factor_init_narrow_comp = []
	index_for_fixing_bounds_factor_init_wide_comp = []
	position_final_narrow_comp = []
	position_final_wide_comp = []
	count_init = 0
	count_paired = 0
	for j in range(len(center_list_init)):
		#bf.print_cust(f'{type(wave)}, {type(center_list_init[j])}, {type(redshift_val)}')
		idx_gen = np.searchsorted(wave, (center_list_init[j]*(1.+redshift_val)))
		if (flux[idx_gen]<=0.0):
			flux[idx_gen] = 0.01
		if (comments_on_balmer[j]):
			if (comments_on_tied[j]):
				if ('tied_init' in comments_on_tied[j]):
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_narrow_components_init)
					position_init_narrow_comp.extend(np.arange(count_init, count_init+number_of_narrow_components_init))
					count_init+=number_of_narrow_components_init
					count_paired+=number_of_narrow_components_init
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_wide_components_init)
					position_init_wide_comp.extend(np.arange(count_init, count_init+number_of_wide_components_init))
					count_init+=number_of_wide_components_init
					count_paired+=number_of_wide_components_init
				elif ('tied_paired' in comments_on_tied[j]):
					position_final_narrow_comp.extend(np.arange(count_paired, count_paired+number_of_narrow_components_init))
					if (number_of_narrow_components_init):
						factor_init_narrow_comp.extend([factor_for_tied[j]])
						index_for_fixing_bounds_factor_init_narrow_comp.extend([j])
					position_final_wide_comp.extend(np.arange(count_paired+number_of_narrow_components_init, count_paired+number_of_narrow_components_init+number_of_wide_components_init))
					if (number_of_wide_components_init):
						factor_init_wide_comp.extend([factor_for_tied[j]])
						index_for_fixing_bounds_factor_init_wide_comp.extend([j])
					count_paired+=(number_of_narrow_components_init+number_of_wide_components_init)
					count_init+=(number_of_narrow_components_init+number_of_wide_components_init)
			else:
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_narrow_components_init)
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_wide_components_init)
		else:
			if (comments_on_tied[j]):
				if ('tied_init' in comments_on_tied[j]):
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_narrow_components_init)
					position_init_narrow_comp.extend(np.arange(count_init, count_init+number_of_narrow_components_init))
					count_init+=number_of_narrow_components_init
					count_paired+=number_of_narrow_components_init
				elif ('tied_paired' in comments_on_tied[j]):
					position_final_narrow_comp.extend(np.arange(count_paired, count_paired+number_of_narrow_components_init))
					if (number_of_narrow_components_init):
						factor_init_narrow_comp.extend([factor_for_tied[j]])
						index_for_fixing_bounds_factor_init_narrow_comp.extend([j])
					count_paired+=number_of_narrow_components_init
					count_init+=number_of_narrow_components_init
			else:
					amplitude_array.extend([np.log10(flux[idx_gen])]*number_of_narrow_components_init)
	amplitude_array.extend(factor_init_narrow_comp)
	amplitude_array.extend(factor_init_wide_comp)
	amplitude_array = np.nan_to_num(amplitude_array, copy=True, nan=1e-6, posinf=1e-6, neginf=1e-6)
	return (amplitude_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp)

def retrieve_all_amplitude_list_rev2(amplitude_array, number_of_narrow_components_init, number_of_wide_components_init, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, comments_on_tied, comments_on_balmer):
	mask1 = [index for index,value in enumerate(comments_on_tied) if 'tied_init' in value]
	mask2 = np.where(np.array(comments_on_balmer)==True)[0]
	if (number_of_wide_components_init):
		wide_component_index = len(list(set(mask1) - (set(mask1) - set(mask2))))
	else:
		wide_component_index = 0

	if (number_of_narrow_components_init):
		narrow_component_index = len(mask1)
	else:
		narrow_component_index = 0
	total_index = []
	if (wide_component_index):
		index_for_wide_factor = amplitude_array[-wide_component_index:]
		if (narrow_component_index):
			index_for_narrow_factor = amplitude_array[-(wide_component_index+narrow_component_index):-(wide_component_index)]
			if (type(index_for_narrow_factor) is not list):
				index_for_narrow_factor = [index_for_narrow_factor]
			index_for_narrow_factor_rev = np.repeat(index_for_narrow_factor,number_of_narrow_components_init)
			total_index.extend(index_for_narrow_factor_rev)
		if (type(index_for_wide_factor) is not list):
			index_for_wide_factor = [index_for_wide_factor]
		index_for_wide_factor_rev = np.repeat(index_for_wide_factor,number_of_wide_components_init)
		total_index.extend(index_for_wide_factor_rev)
	elif (narrow_component_index) and not (wide_component_index):
		index_for_narrow_factor = amplitude_array[-narrow_component_index:]
		if (type(index_for_narrow_factor) is not list):
			index_for_narrow_factor = [index_for_narrow_factor]
		index_for_narrow_factor_rev = np.repeat(index_for_narrow_factor,number_of_narrow_components_init)
		total_index.extend(index_for_narrow_factor_rev)
	rev_amp1 = amplitude_array[:-(wide_component_index+narrow_component_index)]
	rev_amp2 = rev_amp1
	position_init = list(np.append(position_init_narrow_comp, position_init_wide_comp))
	position_final = list(np.append(position_final_narrow_comp, position_final_wide_comp))
	position_final_sorted = [x for _, x in sorted(zip(position_init, position_final))]
	total_index_sorted = [x for _, x in sorted(zip(position_init, total_index))]
	position_init_sorted = list(np.sort(position_init))
	rev_amp3 = np.zeros([len(position_init_sorted)+len(position_final_sorted)])
	counter_test = 0
	for i in range(len(position_final_sorted)):
		rev_amp3[position_final_sorted[i].astype(np.int32)] = rev_amp2[counter_test]+np.log10(total_index_sorted[counter_test])
		rev_amp3[position_init_sorted[i].astype(np.int32)] = rev_amp2[counter_test]
		counter_test+=1
	rev_amp4 = np.nan_to_num(rev_amp3, copy=True, nan=1e-6, posinf=1e-6, neginf=1e-6)
	return(rev_amp4)

##################################GET_OPTIMUM_PARAMETERS##################################
##################################EXTRA_FUNCTION_FOR_CUSTOM_EMISSION_FIT##################################


##################################GET_PHYSICAL_CONDITIONS_EMISSION_LINES##################################
##################################GET_TEMPERATURE_AND_DENSITY_FROM_EMISSION_LINES##################################
def get_temp_and_density_from_emission_lines_pyneb(obs_data_nii_5755, obs_data_err_nii_5755, obs_data_nii_6584, obs_data_err_nii_6584, obs_data_sii_6731, obs_data_err_sii_6731, obs_data_sii_6716, obs_data_err_sii_6716, obs_data_halpha, obs_data_err_halpha, obs_data_hbeta, obs_data_err_hbeta, **kwargs):
	pn_log_level = kwargs.get('pn_log_level', 2)  # Centre List
	redenning_law = kwargs.get('redenning_law', 'CCM89')  # Centre List
	pn.log_.level = pn_log_level
	obs = pn.Observation()
	line5 = pn.EmissionLine(label ="H1r_6563A", obsIntens=obs_data_halpha, obsError=obs_data_err_halpha, corrected=False)
	line6 = pn.EmissionLine(label ="H1r_4861A", obsIntens=obs_data_hbeta, obsError=obs_data_err_hbeta, corrected=False)
	obs.addLine(line5)
	obs.addLine(line6)
	obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
	line1 = pn.EmissionLine('N', 2, 5755, obsIntens=obs_data_nii_5755, obsError=obs_data_err_nii_5755, corrected=False)
	line2 = pn.EmissionLine('N', 2, 6584, obsIntens=obs_data_nii_6584, obsError=obs_data_err_nii_6584, corrected=False)
	line3 = pn.EmissionLine('S', 2, 6731, obsIntens=obs_data_sii_6731, obsError=obs_data_err_sii_6731, corrected=False)
	line4 = pn.EmissionLine('S', 2, 6716, obsIntens=obs_data_sii_6716, obsError=obs_data_err_sii_6716, corrected=False)
	obs.addLine(line1)
	obs.addLine(line2)
	obs.addLine(line3)
	obs.addLine(line4)
	obs.correctData()
	diags = pn.Diagnostics()
	diags.addDiag(['[NII] 5755/6584', '[SII] 6731/6716'])
	Te, Ne = diags.getCrossTemDen('[NII] 5755/6584', '[SII] 6731/6716', obs=obs)
	dust = obs.extinction.E_BV
	return (np.array(Te), np.array(Ne), np.array(dust))
##################################GET_TEMPERATURE_AND_DENSITY_FROM_EMISSION_LINES##################################

##################################GET_METALLICITY_FROM_EMISSION_LINES##################################
def metallicity_from_emission_lines(flux_oiii_5007, flux_err_oiii_5007, flux_hbeta, flux_err_hbeta, flux_nii_6583, flux_err_nii_6583, flux_halpha, flux_err_halpha):
	A = flux_oiii_5007 / flux_hbeta
	A_err = np.abs(A)*(np.sqrt( (flux_err_oiii_5007/flux_oiii_5007)**2 + (flux_err_hbeta/flux_hbeta)**2 ))
	B = flux_nii_6583 / flux_halpha
	B_err = np.abs(B)*(np.sqrt( (flux_err_nii_6583/flux_nii_6583)**2 + (flux_err_halpha/flux_halpha)**2 ))
	C = A / B
	C_err = np.abs(C)*(np.sqrt( (A_err/A)**2 + (B_err/B)**2 ))
	D = np.log10(C)
	D_err = np.abs( C_err / (C*np.log(10)) )
	if (-1.1 < D < 1.7):
		metallicity = D
		metallicity_err = D_err
	else:
		metallicity = B
		metallicity_err = B_err
	return (metallicity, metallicity_err)
##################################GET_METALLICITY_FROM_EMISSION_LINES##################################

##################################GET_DUST_FROM_EMISSION_LINES##################################
def dust_from_emission_lines(flux_hbeta, flux_err_hbeta, flux_halpha, flux_err_halpha):
	A = flux_halpha / flux_hbeta
	A_err = np.abs(A)*(np.sqrt( (flux_err_halpha/flux_halpha)**2 + (flux_err_hbeta/flux_hbeta)**2 ))
	B = A / 2.86
	B_err = B / 2.86
	C = np.log10(B)
	C_err = np.abs( B_err / (B*np.log(10)) )
	D = 1.97*C
	D_err = 1.97*C_err
	return (D, D_err)
##################################GET_DUST_FROM_EMISSION_LINES##################################

##################################GET_SFR_FROM_EMISSION_LINES##################################
def sfr_from_emission_lines(flux_halpha, flux_err_halpha, vel_map, vel_map_err, sig_map, sig_map_err, emitter_redshift):
	flux_halpha = flux_halpha*1e-20
	flux_err_halpha = flux_err_halpha*1e-20
	fitted_area_emm = flux_halpha*1e5
	fitted_area_err_emm = flux_err_halpha*1e5
	dL = LCDM.luminosity_distance(emitter_redshift).to('cm').value
	luminosity = fitted_area_emm*4*np.pi*dL**2
	luminosity_err = fitted_area_err_emm*4*np.pi*dL**2
	SFR = luminosity*4.645e-42
	SFR_err = luminosity_err*4.645e-42
	return (SFR, SFR_err)
##################################GET_SFR_FROM_EMISSION_LINES##################################
##################################GET_PHYSICAL_CONDITIONS_EMISSION_LINES##################################

##################################CUSTOM_EMISSION_FITTING_FUNCTIONS##################################

