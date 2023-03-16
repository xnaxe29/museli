import numpy
import numpy as np
import scipy
import scipy as sp
from scipy import signal
from scipy import interpolate
import warnings
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from os import path
import sys
import time
from astropy.constants import c as c_ms
c_kms = float(c_ms.value / 1e3)
c = c_kms
from astropy.constants import k_B
k = float(k_B.value)
quite_val = False
#c = 299792.458 # Speed in Light in Km/s
#c_kms = 299792.458 # Speed in Light in Km/s
#k = 1.38064852e-23 #Boltzmann's constant in m^2 kg s^(-2) K^(-1)

if not quite_val: warnings.simplefilter(action='ignore', category=FutureWarning)
if not quite_val: np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
if not quite_val: np.seterr(divide = 'ignore')

# Disable Print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore Print
def enablePrint():
    sys.stdout = sys.__stdout__


def print_cust(print_str, **kwargs):
	sleep_time = kwargs.get('sleep_time', 0.1)  # Sleep Time
	quite_val = kwargs.get('quiet_val', False)  # Quite Val
	if not quite_val:
		print ("\n")
		print (print_str)
		print ("\n")
		sys.stdout.flush()
		time.sleep(sleep_time)


#######################CUSTOM_CREATION_OF_COLORBAR_IN_SUBPLOTS_OF_MATPLOTLIB#######################
#This function is called for creating custom colorbar for subplots
#See https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot for more details

def add_colorbar(mappable):
	last_axes = plt.gca()
	ax = mappable.axes
	fig = ax.figure
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.05)
	cbar = plt.colorbar(mappable, cax=cax1)
	cbar.set_ticks(ticker.LogLocator(), update_ticks=True)
	cbar.ax.tick_params(size=0)
	return cbar

def add_colorbar_lin(mappable):
	last_axes = plt.gca()
	ax = mappable.axes
	fig = ax.figure
	divider = make_axes_locatable(ax)
	cax1 = divider.append_axes("right", size="5%", pad=0.05)
	cbar = plt.colorbar(mappable, cax=cax1)
	cbar.set_ticks(ticker.LinearLocator(), update_ticks=True)
	cbar.ax.tick_params(size=0)
	return cbar

#######################CUSTOM_CREATION_OF_COLORBAR_IN_SUBPLOTS_OF_MATPLOTLIB#######################


######################FIND_THE_NEAREST_VALUE##################################
def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return array[idx]

def find_nearest_idx(array,value):
	if type(value) == list:
		idx = []
		for i in range(len(value)):
			idx.extend((np.abs(array-value[i])).argmin())
	elif (isinstance(value, np.ndarray)):
		idx = np.zeros_like(value, dtype=int)
		for i in range(len(value)):
			idx[i] = (np.abs(array-value[i])).argmin()
	else:
		idx = (np.abs(array-value)).argmin()
	return idx

def find_nearest_new(array,value):
	idx = (np.abs(array-value)).argmin()
	return int(idx), float(array[idx])
######################FIND_THE_NEAREST_VALUE##################################


######################ARG_MEDIAN##################################

def arg_median(a):
	if len(a) % 2 == 1:
		return np.where(a == np.nanmedian(a))[0][0]
	else:
		idx, value = find_nearest_new(a, np.nanmedian(a))
		return idx

######################ARG_MEDIAN##################################

########################################FLIP########################################

def flip(m, axis):
    if not hasattr(m, 'ndim'):
        m = asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]

########################################FLIP########################################

####################ARRAY_SMOOTHING_FUNCTION####################
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
####################ARRAY_SMOOTHING_FUNCTION####################

####################RETRIEVE_HEADER_INFORMATION_FROM_WAVE####################
def retrieve_wave_header_info(wave_array, wave_type='linear'):
	if ('log' in wave_type):
		wave_array = np.log10(wave_array)
	crval_cust = float(wave_array[0])
	naxis_cust = int(len(wave_array))
	cdn_n = float(wave_array[1]-wave_array[0])
	cr_pix = 1.
	return (crval_cust, naxis_cust, cdn_n, cr_pix)
####################RETRIEVE_HEADER_INFORMATION_FROM_WAVE####################

########################################STRING_TO_BOOL########################################
def str2bool(v):
	if not type(v)==bool:
		return v.lower() in ("yes", "true", "t", "1")
	else:
		return v
########################################STRING_TO_BOOL########################################


########################################CHANGE_STYLE_OF_RA_DEC########################################
def decdeg2dms(dd):
   is_positive = dd >= 0
   dd = abs(dd)
   minutes,seconds = divmod(dd*3600,60)
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   str_final = str(degrees)+str(r'$^{\deg}$')+str(minutes)+str(r"""$^{"}$""")+str(seconds)+str(r"""$^{'}$""")
   return (str_final)
########################################CHANGE_STYLE_OF_RA_DEC########################################

########################################CHECK_IF_A_DIRECTORY_EXISTS########################################
def check_directory(dir_name):
	if not os.path.exists(dir_name):
		print_cust(f'Making Directory - {dir_name} \n')
		sys.stdout.flush()
		time.sleep(0.1)
		os.makedirs(dir_name)
	else:
		print_cust(f'Directory - {dir_name} exists \n')
		sys.stdout.flush()
		time.sleep(0.1)
########################################CHECK_IF_A_DIRECTORY_EXISTS########################################

##################################CHEBYSHEV_FUNCTIONS##################################
def chebyshev_order(wave, cont, stopping_number):
    wave_new = np.linspace(-1, 1, len(wave))
    i=1
    while True:
        roots = numpy.polynomial.chebyshev.chebfit(wave_new, cont, i, rcond=None, full=False, w=None)
        poly = numpy.polynomial.chebyshev.chebval(wave_new, roots, tensor=True)
        chi_sq = ((poly - cont) ** 2)
        chi_sq_sum = (np.sum(chi_sq))/len(cont)
        i+=1
        if chi_sq_sum<(float(stopping_number)):
            break
    return (i, roots)

def chebyshev_fit(wave, cont, order_test):
    wave_new = np.linspace(-1, 1, len(wave))
    roots = numpy.polynomial.chebyshev.chebfit(wave_new, cont, int(order_test), rcond=None, full=False, w=None)
    poly_new = numpy.polynomial.chebyshev.chebval(wave_new, roots, tensor=True)
    return(roots, poly_new)

def chebyshev_disp(wave, coeff):
    wave_new = np.linspace(-1, 1, len(wave))
    poly_new = numpy.polynomial.chebyshev.chebval(wave_new, coeff, tensor=True)
    return(poly_new)
##################################CHEBYSHEV_FUNCTIONS##################################

##################################GAUSS_BASED_FUNCTIONS##################################
def vel_prof(x, centre):
    xnew = c_kms * ((x-centre)/x)
    return (xnew)

def wave_prof(vel_center, centre):
    xnew = (centre*c_kms) / (c_kms-vel_center)
    return (xnew)

def gaus_prof_wave(wave, wave_rest, log_amp, center, sigma):
	vel_array = vel_prof(wave, wave_rest)
	prof = (10**amp)*np.exp(-(vel_array-center)**2./(2.*sigma**2.))
	return (prof)

def gaus_prof_vel(vel_array, amp, center, sigma):
	prof = (10**amp)*np.exp(-(vel_array-center)**2./(2.*sigma**2.))
	return (prof)
##################################GAUSS_BASED_FUNCTIONS##################################


############################CONVOLUTION#####################################
def convolved_prof5(wave, profile, res):
	#res = 3026.62381895424
	wave_short = np.logspace(np.log10(wave.min()), np.log10(wave.max()), len(wave)*10)
	center = np.searchsorted(wave_short, np.log10(np.median(wave_short)))
	deltalam = wave_short[center + 1] - wave_short[center]
	sigma = (wave_short[center]) / (res*2. * (2 * np.sqrt(2 * np.log(2))) * deltalam)
	gauss = scipy.signal.gaussian(len(profile), sigma, sym=True)
	gauss = gauss / np.sum(gauss)
	prof_new = signal.fftconvolve(profile, gauss, mode='same')
	return (prof_new)
############################CONVOLUTION#####################################

############################AIR_TO_VACUUM_CONVERSION############################

def _wave_convert(lam):
    """
    Convert between vacuum and air wavelengths using
    equation (1) of Ciddor 1996, Applied Optics 35, 1566
        http://doi.org/10.1364/AO.35.001566

    :param lam - Wavelength in Angstroms
    :return: conversion factor

    """
    lam = np.asarray(lam)
    sigma2 = (1e4/lam)**2
    fact = 1 + 5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)

    return fact

def vac_to_air(lam_vac):
    """
    Convert vacuum to air wavelengths

    :param lam_vac - Wavelength in Angstroms
    :return: lam_air - Wavelength in Angstroms

    """
    return lam_vac/_wave_convert(lam_vac)

def air_to_vac(lam_air):
    """
    Convert air to vacuum wavelengths

    :param lam_air - Wavelength in Angstroms
    :return: lam_vac - Wavelength in Angstroms

    """
    return lam_air*_wave_convert(lam_air)

############################AIR_TO_VACUUM_CONVERSION############################


##################################LOG_TO_REAL_ERROR##################################
def log_to_real_err(log_array, log_error):
	linear_array = 10**(log_array)
	linear_err = np.abs(linear_array * (np.log(10) * log_error))
	return linear_err
##################################LOG_TO_REAL_ERROR##################################

##################################REAL_TO_LOG_ERROR##################################
def real_to_log_err(linear_array, linear_error):
	log_array = np.log10(linear_array)
	log_err = np.abs(linear_error/(linear_array*(np.log(10))))
	return log_err
##################################REAL_TO_LOG_ERROR##################################

##################################SMOOTHING_FUNCTION##################################

def smooth_with_error(y, y_err, box_pts):
	box = np.ones(box_pts)/box_pts
	y_err_rev = y_err / y
	y_smooth = np.convolve(y, box, mode='same')
	y_err_rev_smoothed_1 = np.convolve(y_err_rev, box, mode='same')
	y_err_smoothed = y_err_rev_smoothed_1*y_smooth
	return y_smooth, y_err_smoothed

##################################SMOOTHING_FUNCTION##################################


##################################CLEAN_DATA##################################

def clean_data(data, **kwargs):
	data_type = kwargs.get('type_of_data', 'data')  # Initial value for Velocity
	data_val = kwargs.get('val_data', 1e-6)  # Initial value for Velocity
	data = data.astype(np.float64)
	data_1 = np.nan_to_num(data, nan=data_val, posinf=data_val, neginf=-data_val)
	mask = np.where(data_1[:]=='')
	data_1[mask] = data_val
	mask2 = np.where(data_1[:]==' ')
	data_1[mask2] = data_val
	if 'err' in data_type:
		data_1[data_1<=0.]=data_val*10.
	return (data_1)

##################################CLEAN_DATA##################################


######################GET_EQUIVALENT_WIDTH##################################
def get_equivalendth_width(wave_full, data_full, err_full, cont_full, wave_start, wave_end, redshift):
	wave_rest = wave_full / (1.+redshift)
	idx_start = find_nearest_idx(wave_rest, wave_start)
	idx_end = find_nearest_idx(wave_rest, wave_end)
	wave = wave_rest[idx_start:idx_end]
	data = data_full[idx_start:idx_end] / cont_full[idx_start:idx_end]
	err = err_full[idx_start:idx_end] / cont_full[idx_start:idx_end]
	data_new = 1. - data
	err_new = err
	delta_lambda = wave[1] - wave[0]
	equivalent_width = (np.sum(data_new))*(delta_lambda)
	equivalent_width_err = (np.sqrt(np.sum(err_new**2)))*(delta_lambda)
	return (equivalent_width, equivalent_width_err)

def get_equivalendth_width_rev(wave_rest, data_norm, err_norm, center=6563., redshift=0.0, vel_window=2000.):
	vel_array = vel_prof(wave_rest, center)
	idx_start = find_nearest_idx(vel_array, -vel_window)
	idx_end = find_nearest_idx(vel_array, vel_window)
	wave = wave_rest[idx_start:idx_end]
	data = 1. - data_norm[idx_start:idx_end]
	err = err_norm[idx_start:idx_end]
	delta_lambda = float(wave[1] - wave[0])
	assert delta_lambda!=0.0, f"Delta Lambda: {delta_lambda} cannot be zero..."
	equivalent_width = (np.sum(data))*(delta_lambda)
	equivalent_width_err = (np.sqrt(np.nansum(err)))*(delta_lambda)
	return (equivalent_width, equivalent_width_err)

def get_equivalendth_width_rev_new(wave_rest, data_original, err_original, lick_index_array=np.array([6522., 6542., 6554., 6570., 6585., 6625.])):
	index_w11, index_w12, index_m11, index_m12, index_w21, index_w22 = find_nearest_idx(wave_rest, lick_index_array)
	wave_to_fit = np.append(wave_rest[index_w11:index_w12], wave_rest[index_w21:index_w22])
	flux_to_fit = np.append(data_original[index_w11:index_w12], data_original[index_w21:index_w22])
	flux_err_to_fit = np.append(err_original[index_w11:index_w12], err_original[index_w21:index_w22])
	f_fit_flux = interpolate.interp1d(wave_to_fit, flux_to_fit, axis=0, fill_value="extrapolate", kind='linear')
	norm_wave = wave_rest[index_w11:index_w22]
	norm_flux = data_original[index_w11:index_w22] / f_fit_flux(wave_rest[index_w11:index_w22])
	norm_flux_err = err_original[index_w11:index_w22] / f_fit_flux(wave_rest[index_w11:index_w22])
	norm_flux_rev = 1. - norm_flux
	delta_lambda = float(norm_wave[1] - norm_wave[0])
	assert delta_lambda!=0.0, f"Delta Lambda: {delta_lambda} cannot be zero..."
	equivalent_width = (np.nansum(norm_flux_rev))*(delta_lambda)
	equivalent_width_err = (np.sqrt(np.nansum(norm_flux_err)))*(delta_lambda)
	return (equivalent_width, equivalent_width_err)


######################GET_EQUIVALENT_WIDTH##################################


######################GET_INITIAL_KINEMATICS##################################
def get_initial_kinematics_in_kms(number_of_narrow_components_init, number_of_wide_components_init):
	velocity_init_narrow = np.linspace(-200., 200., number_of_narrow_components_init)
	velocity_init_wide = np.append(velocity_init_narrow, np.linspace(-200., 200., number_of_wide_components_init))
	#center_init_array = velocity_init_wide
	vel_disp_init_narrow = np.full([number_of_narrow_components_init], 100.)
	vel_disp_init_wide = np.append(vel_disp_init_narrow, np.full([number_of_wide_components_init], 300.))
	#sigma_init_array = vel_disp_init_wide
	return(velocity_init_wide, vel_disp_init_wide)
######################GET_INITIAL_KINEMATICS##################################
