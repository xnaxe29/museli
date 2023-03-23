import sys
import os
import numpy as np
import spectres
from scipy import interpolate
from astropy.io import fits
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import ast
import pandas as pd
import scipy as scp
from scipy.signal import argrelextrema
from scipy.ndimage.morphology import binary_erosion
from scipy.signal import find_peaks
from scipy.signal import savgol_filter
from scipy import interpolate
from scipy.interpolate import interp1d
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
from sklearn.cluster import DBSCAN
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.optimize import curve_fit

#from basic_functions import *
import basic_functions as bf
import handling_fits_files as hff
from astropy.constants import c as c_ms
c_kms = float(c_ms.value / 1e3)
quite_val = False



##################################GET_DATA_FROM_FILE##################################


def get_binned_data(file1, sky_file=None, min_snr=3.0, quiet_val=False):
	wave, galaxy, noise = file1
	snr_array = galaxy/noise
	wave_smooth_point = len(wave)
	smoothing_count = 1
	if sky_file is not None:
		sky_wave, sky_flux, sky_flux_err = sky_file
	while (np.nanmin(snr_array) < min_snr):
		bf.print_cust(f'{np.nanmin(snr_array)}', quiet_val=quiet_val)
		smoothing_count+=1
		galaxy_smoothed, noise_smoothed = bf.smooth_with_error(galaxy, noise, smoothing_count)
		new_wavs = np.logspace(np.log10(wave.min()), np.log10(wave.max()), num=wave_smooth_point, endpoint=True, base=10)
		galaxy_new, noise_new = spectres.spectres(new_wavs, wave, galaxy_smoothed, spec_errs=noise_smoothed, fill=np.nan, verbose=True)
		if sky_file is not None:
			sky_flux_new, sky_flux_err_new = spectres.spectres(new_wavs, sky_wave, sky_flux, spec_errs=sky_flux_err, fill=np.nan, verbose=True)

		snr_array = galaxy_new/noise_new
	galaxy_new = bf.clean_data(galaxy_new)
	noise_new = bf.clean_data(noise_new, type_of_data='err')
	file2 = [new_wavs, galaxy_new, noise_new]
	if sky_file is not None:
		sky_flux_new = bf.clean_data(sky_flux_new)
		sky_flux_err_new = bf.clean_data(sky_flux_err_new, type_of_data='err')
		file_sky = [new_wavs, sky_flux_new, sky_flux_err_new]
	else:
		file_sky = None
	return (file2, file_sky)



def combined_function_rebinning_continuum(wave_orig, data_real, err_real, spectral_smoothing_int=1, redshift_val=0.0, quiet_val=False, continuum_est_type='custom'):
	wave_rebinned, flux_rebinned, flux_err_rebinned, fwhm_gal_init, velscale = refine_obs_data_using_scipy(wave_orig, data_real, err_real, spectral_smoothing=spectral_smoothing_int)
	continuum_line_type = check_continuum_type(wave_rebinned, flux_rebinned, redshift=redshift_val)
	inst_cont = continuum_fitClass()
	inst_cont.flux = flux_rebinned
	inst_cont.continuum_finding_method = str(continuum_est_type)
	continuum_rebinned = inst_cont.continuum_finder_flux(wave_rebinned)
	#continuum_rebinned = get_initial_continuum_rev_2(wave_rebinned, flux_rebinned, line_type=str(continuum_line_type), printing=quiet_val)
	return (flux_rebinned, flux_err_rebinned, fwhm_gal_init, velscale, continuum_rebinned)
	
def check_continuum_type(wave, flux, redshift, vel_window=2000.):
	redshifted_halpha = (1.+redshift)*6563.
	if (np.nanmin(wave)<=redshifted_halpha<=np.nanmax(wave)):
		vel_array = bf.vel_prof(wave, redshifted_halpha)
	else:
		vel_array = bf.vel_prof(wave, wave[int(len(wave)/2)])

	idx1 = bf.find_nearest_idx(vel_array,-float(vel_window))
	idx2 = bf.find_nearest_idx(vel_array,float(vel_window))
	min, mean, max = np.nanmin(flux[idx1:idx2]), np.nanmean(flux[idx1:idx2]), np.nanmax(flux[idx1:idx2])
	if (np.abs(min-mean) < np.abs(max-mean)):
		cont_line_type = 'emission'
	else:
		cont_line_type = 'absorption'
	return (cont_line_type)

def refine_obs_data_using_scipy(wave_orig, data_real, err_real, spectral_smoothing=1):
	wave_rebinned = np.logspace(np.log10(wave_orig.min()), np.log10(wave_orig.max()), num=int(len(wave_orig)/int(spectral_smoothing)), endpoint=True, base=10)
	f_flux = interpolate.interp1d(wave_orig, data_real, axis=0, fill_value="extrapolate", kind='cubic')
	f_flux_err = interpolate.interp1d(wave_orig, err_real, axis=0, fill_value="extrapolate", kind='cubic')
	flux_rebinned = f_flux(wave_rebinned)
	flux_err_rebinned = f_flux_err(wave_rebinned)
	#flux_rebinned_smoothed, flux_err_rebinned_smoothed = bf.smooth_with_error(flux_rebinned, flux_err_rebinned, int(spectral_smoothing))
	frac = wave_rebinned[1]/wave_rebinned[0]    # Constant lambda fraction per pixel
	dlam_gal = (frac - 1)*wave_rebinned            # Size of every pixel in Angstrom
	wdisp = np.ones([len(wave_rebinned)])          # Intrinsic dispersion of every pixel, in pixels units
	fwhm_gal_init = 2.355*wdisp*dlam_gal              # Resolution FWHM of every pixel, in Angstroms.
	fwhm_gal = np.nanmean(fwhm_gal_init)            # Keeping it as mean of fwhm_gal as we need a specific number.
	velscale = c_kms*np.log(wave_rebinned[1]/wave_rebinned[0])
	return (wave_rebinned, flux_rebinned, flux_err_rebinned, fwhm_gal_init, velscale)


def get_data_from_files(file1, file_type, require_air_to_vaccum, extra_redshift, quiet=True, **kwargs):
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
		wave_ref *= np.median(bf.vac_to_air(wave_ref)/wave_ref)
	FWHM_gal = np.nanmean(fwhm_gal_init)
	bf.print_cust(f'Velscale: {velscale}, FWHM Gal: {FWHM_gal}', quiet_val=quiet)
	fit_type = kwargs.get('fit_type', 'custom') # fit_type
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
def get_data_from_file_custom(file1, file_type, require_air_to_vaccum, extra_redshift, quiet=True, **kwargs):
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
		wave_ref *= np.median(bf.vac_to_air(wave_ref)/wave_ref)
	FWHM_gal = np.nanmean(fwhm_gal_init)
	bf.print_cust(f'Velscale: {velscale}, FWHM Gal: {FWHM_gal}', quiet_val=quiet)
	fit_type = kwargs.get('fit_type', 'custom') # fit_type
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

##################################GET_CONTINUUM##################################


class continuum_fitClass:
	def __init__(self):
		self.plot=False
		self.print=False
		self.legend=False
		self.ax=False
		self.label=False
		self.color='gray'
		self.default_smooth=5
		self.default_order=8
		self.default_allowed_percentile=75
		self.default_filter_points_len=10
		self.data_height_upscale=2
		self.default_poly_order=3
		self.default_window_size_default=999
		self.default_fwhm_galaxy_min=10
		self.default_noise_level_sigma=10
		self.default_fwhm_ratio_upscale=10
		#self.spline_smoothing_factor = 1
		self.pca_component_number = 1
		self.lowess_cont_frac = 0.05
		self.gaussian_cont_fit = 200
		self.peak_prominence = 0.05
		self.peak_width = 10
		self.median_filter_window = 101
		self.continuum_finding_method = 'custom'
		pass
	def smooth(self, y, box_pts):
		box = np.ones(box_pts)/box_pts
		y_smooth = np.convolve(y, box, mode='same')
		return y_smooth
	def gaussian(self, x, amp, mu, sig):
		return amp*(1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2))
	#This function has been created taking inference from Martin+2021 (2021MNRAS.500.4937M)
	def continuum_finder(self, wave, *pars, **kwargs):
		flux = self.flux
		#n_smooth=5, order=8, allowed_percentile=75, poly_order=3, window_size_default=None
		if pars:
			#n_smooth, allowed_percentile, poly_order, window_size_default, fwhm_galaxy_min = pars
			#print (pars)
			n_smooth, allowed_percentile, filter_points_len, data_height_upscale, poly_order, window_size_default, fwhm_galaxy_min, noise_level_sigma, fwhm_ratio_upscale = pars
		else:
			#n_smooth, allowed_percentile, poly_order, window_size_default, fwhm_galaxy_min = np.array([self.default_smooth, self.default_allowed_percentile, self.default_poly_order, self.default_window_size_default, self.default_fwhm_galaxy_min])
			n_smooth, allowed_percentile, filter_points_len, data_height_upscale, poly_order, window_size_default, fwhm_galaxy_min, noise_level_sigma, fwhm_ratio_upscale = np.array([self.default_smooth, self.default_allowed_percentile, self.default_filter_points_len, self.data_height_upscale, self.default_poly_order, self.default_window_size_default, self.default_fwhm_galaxy_min, self.default_noise_level_sigma, self.default_fwhm_ratio_upscale])

		n_smooth = int(n_smooth)
		#order = int(order)
		allowed_percentile = int(allowed_percentile)
		filter_points_len = int(filter_points_len)
		poly_order = int(poly_order)
		window_size_default = int(window_size_default)
		fwhm_galaxy_min = int(fwhm_galaxy_min)
		noise_level_sigma = int(noise_level_sigma)
		fwhm_ratio_upscale = int(fwhm_ratio_upscale)

		pick = np.isfinite(flux) #remove NaNs
		flux = flux[pick] #remove NaNs
		#smoothed_data = scp.ndimage.convolve1d(flux, np.asarray([1.]*n_smooth)/n_smooth) #smooth data
		smoothed_data = self.smooth(flux, int(n_smooth))
		local_std = np.median([ np.std(s) for s in np.array_split(flux, int(n_smooth)) ]) #find local standard deviation
		mask_less = argrelextrema(smoothed_data, np.less)[0] #find relative extreme points in absorption
		mask_greater = argrelextrema(smoothed_data, np.greater)[0] #find relative extreme points in emission
		mask_less_interpolate_func = interpolate.interp1d(wave[mask_less], flux[mask_less], kind='cubic', fill_value="extrapolate") #interpolate wavelength array like function from relative extreme points in absorption
		mask_greater_interpolate_func = interpolate.interp1d(wave[mask_greater], flux[mask_greater], kind='cubic', fill_value="extrapolate") #interpolate wavelength array like function from relative extreme points in emission
		absolute_array = mask_greater_interpolate_func(wave)-mask_less_interpolate_func(wave) #obtain the absolute array for find_peaks algorithm
		filter_points = np.array([int(i*len(absolute_array)/filter_points_len) for i in range(1,filter_points_len)])
		noise_height_max_default = noise_level_sigma*np.nanmin(np.array([np.nanstd(absolute_array[filter_points[i]-10:filter_points[i]+10]) for i in range(len(filter_points))]))
		data_height_max_default = data_height_upscale*np.nanmax(np.array([np.abs(np.nanmax(absolute_array)), np.abs(np.nanmin(absolute_array))]))
		noise_height_max = kwargs.get('noise_height_max', noise_height_max_default)  # Maximal height for noise
		data_height_max = kwargs.get('data_height_max', data_height_max_default)  # Maximal height for data
		peaks = find_peaks(absolute_array, height=[noise_height_max, data_height_max], prominence=(local_std*3.), width = [fwhm_galaxy_min, int(fwhm_ratio_upscale*fwhm_galaxy_min)]) #run scipy.signal find_peaks algorithm to find peak points
		edges = np.int32([np.round(peaks[1]['left_ips']), np.round(peaks[1]['right_ips'])]) #find edges of peaks
		d = (np.diff(flux, n=1))
		w = 1./np.concatenate((np.asarray([np.median(d)]*1),d))
		w[0] = np.max(w)
		w[-1] = np.max(w)
		for edge in edges.T:
			#print (wave[pick][edge[0]], wave[pick][edge[1]])
			diff_tmp = int((edge[1] - edge[0])/2)
			w[edge[0]-diff_tmp:edge[1]+diff_tmp] = 1./10000.
		w = np.abs(w)
		pick_2 = np.where(w > np.percentile(w, allowed_percentile * (float(len(flux)) / float(len(wave)))))[0]
		#fit = np.poly1d(np.polyfit(tsteps[pick][pick_2], a[pick_2], order))
		#print (wave[pick][pick_2])
	
		if len(wave[pick][pick_2])>3:
			xx = np.linspace(np.min(wave[pick][pick_2]), np.max(wave[pick][pick_2]), 1000)
			itp = interpolate.interp1d(wave[pick][pick_2], flux[pick_2], kind='linear')

		else:
			mask = np.ones_like(wave, dtype=np.bool8)
			ynew = np.abs(np.diff(flux[mask], prepend=1e-10))
			ynew2 = np.percentile(ynew, allowed_percentile)
			xx = wave[mask][ynew < ynew2]
			y_rev = flux[mask][ynew < ynew2]
			itp = interpolate.interp1d(xx, y_rev, axis=0, fill_value="extrapolate", kind='linear')
			#y_rev2 = f_flux(wave)

		#window_size = int(((1.0 / (step_to_t[pick][pick_2][-1] - step_to_t[pick][pick_2][0])) * 1000.))
		#window_size = window_size_default
		window_size = int(fwhm_ratio_upscale*fwhm_galaxy_min)
		if window_size % 2 == 0:
			window_size = window_size + 1
		fit_savgol = savgol_filter(itp(xx), window_size, poly_order)
		fit = interpolate.interp1d(xx, fit_savgol, kind='cubic', fill_value="extrapolate")
		#std_cont = np.std((a[pick_2] - fit(tsteps[pick][pick_2])) / fit(tsteps[pick][pick_2]))
		std_cont = np.std(flux[pick_2] - fit(wave[pick][pick_2]))
		#r = (a[pick_2] - fit(tsteps[pick][pick_2])) / fit(tsteps[pick][pick_2])
		#r = flux[pick_2] - fit(wave[pick][pick_2])
		return fit, std_cont, flux, pick, pick_2, peaks, std_cont


	# Function for fitting a polynomial continuum
	def poly_fit_continuum(self, wave):
		degree = self.default_poly_order
		flux = self.flux
		wavelength = wave
		coefficients = np.polyfit(wavelength, flux, degree)
		continuum = np.polyval(coefficients, wavelength)
		return continuum



	# Function for fitting a spline continuum
	def spline_fit_continuum(self, wave):
		flux = self.flux
		wavelength = wave
		#smoothing_factor = self.spline_smoothing_factor
		continuum = interp1d(wavelength, flux, kind='cubic')
		return continuum(wavelength)



	# Function for estimating the continuum using PCA
	def pca_continuum(self, wave):
		flux = self.flux
		wavelength = wave
		pca_component_number = self.pca_component_number
		pca = PCA(n_components = pca_component_number)
		X = flux.reshape(-1, 1)
		pca.fit(X)
		continuum = pca.inverse_transform(pca.transform(X)).flatten()
		return continuum



	# Function for estimating the continuum using a lowess smoother
	def lowess_continuum(self, wave):
		wavelength = wave
		flux = self.flux
		lowess_continuum_fraction = self.lowess_cont_frac
		continuum = lowess(flux, wavelength, frac=lowess_continuum_fraction)[:, 1]
		return continuum


	# Function for estimating the continuum using a Gaussian fit
	def gaussian_func(self, x, a, b, c, d):
		return a * np.exp(-((x-b)/c)**2) + d

	def gaussian_fit_continuum(self, wave):
		wavelength = wave
		flux = self.flux
		window = self.gaussian_cont_fit
		continuum = np.zeros_like(flux)
		for i in range(len(flux)):
			low = max(0, i-window//2)
			high = min(len(flux), i+window//2)
			x = wavelength[low:high]
			y = flux[low:high]
			try:
				popt, _ = curve_fit(self.gaussian_func, x, y, p0=[1, wavelength[i], 5, 0])
				continuum[i] = self.gaussian_func(wavelength[i], *popt)
			except RuntimeError:
				continuum[i] = np.nan
		mask = np.isnan(continuum)
		continuum[mask] = np.interp(wavelength[mask], wavelength[~mask], continuum[~mask])
		return continuum


	# Function for estimating the continuum using peak finding
	def peak_find_continuum(self, wave):
		wavelength = wave
		flux = self.flux
		prominence = self.peak_prominence
		width = self.peak_width
		peaks, _ = find_peaks(flux, prominence=prominence, width=width)
		troughs, _ = find_peaks(-flux, prominence=prominence, width=width)
		indices = np.concatenate([peaks, troughs, [0, len(flux)-1]])
		continuum = np.interp(wavelength, wavelength[indices], flux[indices])
		return continuum


	# Function for estimating the continuum using median filtering
	def median_filter_continuum(self, wave):
		wavelength = wave
		flux = self.flux
		window = self.median_filter_window
		continuum = np.zeros_like(flux)
		for i in range(len(flux)):
			low = max(0, i-window//2)
			high = min(len(flux), i+window//2)
			continuum[i] = np.nanmedian(flux[low:high])
		mask = np.isnan(continuum)
		continuum[mask] = np.interp(wavelength[mask], wavelength[~mask], continuum[~mask])
		return continuum

	# Define distance metric function
	def dist_metric(self, x, y):
		return np.abs(x - y)

	def continuum_using_fof(self, wave):
		wavelength = wave
		flux = self.flux
		# Calculate distance matrix
		X = wavelength.reshape(-1, 1)
		D = cdist(X, X, self.dist_metric)
		# Find clusters using DBSCAN
		db = DBSCAN(eps=3, min_samples=3, metric='precomputed').fit(D)
		labels = db.labels_
		# Calculate continuum
		mask = (labels == 0)
		continuum = np.median(flux[mask])
		return continuum


	def continuum_finder_flux(self, wave, *pars, **kwargs):
		cont_find_method = self.continuum_finding_method
		if (cont_find_method=='custom'):
			fit, std_cont, flux, pick, pick_2, peaks, std_cont = self.continuum_finder(wave, *pars)
			cont_flux = fit(wave[pick])
		elif (cont_find_method=='poly'):
			cont_flux = self.poly_fit_continuum(wave)
		elif (cont_find_method=='spline'):
			cont_flux = self.spline_fit_continuum(wave)
		elif (cont_find_method=='pca'):
			cont_flux = self.pca_continuum(wave)
		elif (cont_find_method=='lowess'):
			cont_flux = self.lowess_continuum(wave)
		elif (cont_find_method=='gauss'):
			cont_flux = self.gaussian_fit_continuum(wave)
		elif (cont_find_method=='peak_find'):
			cont_flux = self.peak_find_continuum(wave)
		elif (cont_find_method=='median_filtering'):
			cont_flux = self.median_filter_continuum(wave)
		elif (cont_find_method=='fof'):
			continuum_fof = self.continuum_using_fof(wave)
			cont_flux = continuum_fof*np.ones_like(wave)
		else:
			bf.print_cust(f'Continuum finding method: {cont_find_method} not found. Reverting back to custom method.')
			fit, std_cont, flux, pick, pick_2, peaks, std_cont = self.continuum_finder(wave, *pars)
			cont_flux = fit(wave[pick])
		return (cont_flux)







#This function has been created from literature (Martin+2021, 2021MNRAS.500.4937M)
def get_initial_continuum_rev_2(wave, flux, smoothing_par=5, allowed_percentile=50, savgol_filter_window_length=500, savgol_filter_polyorder=5, fwhm_galaxy=10, line_type='emission', plot=False, printing=False):
	if (plot):
		plt.plot(wave, flux, zorder=1)
	local_std = np.median([ np.std(s) for s in np.array_split(flux, smoothing_par) ])
	#fwhm_index = int(np.ceil((3.*fwhm_galaxy / np.mean(np.diff(wave)))))
	fwhm_index = int(fwhm_galaxy)
	flux_rev = bf.smooth(flux, smoothing_par)
	if line_type=='absorption':
		bf.print_cust('Estimating continuum for absorption lines.', quiet_val=printing)
		#peaks, dict_cust = find_peaks(flux_rev, width = [fwhm_index, None], prominence=(None, local_std), height=[None,0])
		peaks, dict_cust = find_peaks(flux_rev, width = [fwhm_index, None], prominence=(None, local_std), height=[0,None])
	else:
		bf.print_cust('Estimating continuum for emission lines.', quiet_val=printing)
		peaks, dict_cust = find_peaks(flux_rev, width = [fwhm_index, None], prominence=(None, local_std), height=[0,None])

	left_edge = dict_cust['left_bases']
	right_edge = dict_cust['right_bases']
	left_width = dict_cust['left_ips']
	right_width = dict_cust['right_ips']
	mask = np.ones_like(wave, dtype=np.bool8)
	for i in range(len(left_width)):
		dw = (int(right_width[i]) - int(left_width[i]))
		mask[int(left_width[i])-dw:int(right_width[i])+dw] = False
	if not (len(mask[mask==True])):
		mask[int(len(mask)/2.)-int(len(mask)/10.):int(len(mask)/2.)+int(len(mask)/10.)] = True
	ynew = np.abs(np.diff(flux[mask], prepend=1e-10))
	ynew2 = np.percentile(ynew, allowed_percentile)
	x_rev = wave[mask][ynew < ynew2]
	y_rev = flux[mask][ynew < ynew2]
	if (plot):
		plt.plot(x_rev, y_rev, 'ro', zorder=4)

	if (savgol_filter_window_length <= savgol_filter_polyorder):
		savgol_filter_window_length = savgol_filter_polyorder+1
	if (savgol_filter_window_length % 2) == 0:
		savgol_filter_window_length+=1

	if len(x_rev)<3:
		y_rev3 = np.full([len(flux)], fill_value=1e-6)
	else:
		f_flux = interpolate.interp1d(x_rev, y_rev, axis=0, fill_value="extrapolate", kind='linear')
		y_rev2 = f_flux(wave)
		y_rev3 = savgol_filter(y_rev2, window_length = savgol_filter_window_length, polyorder=savgol_filter_polyorder)
	if (plot):
		#plt.plot(wave, flux, zorder=1)
		#plt.plot(x_rev, y_rev, 'g--', zorder=4)
		plt.plot(wave, y_rev3, 'r--', zorder=5)
		plt.show()
	return (y_rev3)

##################################GET_CONTINUUM##################################

##################################GET_DATA_FROM_FILE##################################




########################################CHECK_FOR_SKY########################################
def check_for_sky(par_dict, wave_rebinned_central, quiet_val=False):
	if ('sky_spectra_file' in par_dict):
		assert os.path.exists(str(par_dict['sky_spectra_file'])), f"File: {str(par_dict['sky_spectra_file'])} not found..."
		bf.print_cust(f"Sky spectra found: {str(par_dict['sky_spectra_file'])}", quiet_val=quiet_val)
		if ('fits' in str(par_dict['sky_spectra_file'])):
			header_muse_sky, header_muse_sky_err, data_muse_sky_init, err_muse_sky_init = hff.open_ifu_fits(str(par_dict['sky_spectra_file']))
			data_muse_sky = bf.clean_data(data_muse_sky_init, val_data=np.nanmin(data_muse_sky_init))
			err_muse_sky = bf.clean_data(err_muse_sky_init, type_of_data='err', val_data=np.nanmin(err_muse_sky_init))
			sky_ra_muse, sky_dec_muse, sky_wave_muse = hff.obtain_physical_axis(header_muse_sky)
			data_muse_sky_flattened = data_muse_sky.reshape(data_muse_sky.shape[0], -1)
			err_muse_sky_flattened = err_muse_sky.reshape(err_muse_sky.shape[0], -1)
			bar = IncrementalBar('Countdown', max = int(data_muse_sky_flattened.shape[0]))
			flux_sky_muse = np.zeros([data_muse_sky_flattened.shape[0]])
			flux_err_sky_muse = np.zeros([data_muse_sky_flattened.shape[0]])
			for i in range(data_muse_sky_flattened.shape[0]):
				bar.next()
				idx_flux_sky_muse = bf.arg_median(data_muse_sky_flattened[i,:])
				flux_sky_muse[i] = data_muse_sky_flattened[i,idx_flux_sky_muse]
				flux_err_sky_muse[i] = np.sqrt(err_muse_sky_flattened[i,idx_flux_sky_muse])
		
			wave_sky_full_rebinned, flux_sky_full_rebinned, flux_err_sky_full_rebinned, sky_fwhm_gal_init, sky_velscale = refine_obs_data_using_scipy(sky_wave_muse, flux_sky_muse, flux_err_sky_muse, spectral_smoothing=int(par_dict['spectral_smoothing']))
			if not ((len(flux_sky_full_rebinned)==len(wave_rebinned_central)) and (len(flux_err_sky_full_rebinned)==len(wave_rebinned_central))):
				bf.print_cust(f"Mismatch between data and sky rebinned length: Data: {len(wave_rebinned_central)}, Sky: {len(flux_sky_full_rebinned)}. Setting sky to zeros with shape of the wavelength array.", quiet_val=quiet_val)
				wave_sky_full_rebinned = wave_rebinned_central
				flux_sky_full_rebinned = np.zeros([len(wave_rebinned_central)])
				flux_err_sky_full_rebinned = np.ones([len(wave_rebinned_central)])

		else:
			#wave_sky_full, flux_sky_full, flux_err_sky_full = np.loadtxt(str(par_dict['sky_spectra_file']), unpack=True)
			data_sky = pd.read_csv(str(par_dict['sky_spectra_file']), sep="\t", comment='#', encoding=None)
			wave_sky_full = data_sky['wave']
			flux_sky_full = data_sky['sky_flux']
			flux_err_sky_full = data_sky['sky_flux_err']
			wave_sky_full_rebinned, flux_sky_full_rebinned, flux_err_sky_full_rebinned, sky_fwhm_gal_init, sky_velscale = refine_obs_data_using_scipy(wave_sky_full, flux_sky_full, flux_err_sky_full, spectral_smoothing=int(par_dict['spectral_smoothing']))
			if not ((len(flux_sky_full_rebinned)==len(wave_rebinned_central)) and (len(flux_err_sky_full_rebinned)==len(wave_rebinned_central))):
				bf.print_cust(f"Mismatch between data and sky rebinned length: Data: {len(wave_rebinned_central)}, Sky: {len(flux_sky_full_rebinned)}. Setting sky to zeros with shape of the wavelength array.", quiet_val=quiet_val)
				wave_sky_full_rebinned = wave_rebinned_central
				flux_sky_full_rebinned = np.zeros([len(wave_rebinned_central)])
				flux_err_sky_full_rebinned = np.ones([len(wave_rebinned_central)])
	else:
		bf.print_cust(f"Sky spectra not found. Setting sky to zeros with shape of the wavelength array.", quiet_val=quiet_val)
		wave_sky_full_rebinned = wave_rebinned_central
		flux_sky_full_rebinned = np.zeros([len(wave_rebinned_central)])
		flux_err_sky_full_rebinned = np.ones([len(wave_rebinned_central)])
	return (wave_sky_full_rebinned, flux_sky_full_rebinned, flux_err_sky_full_rebinned)
########################################CHECK_FOR_SKY########################################





########################################INITIAL_GUESS_FROM_PARAMETER_FILE########################################

def initial_guess_from_parfile(parameter_file_string_current, quiet_val=False):
	assert os.path.exists(parameter_file_string_current), f"File: {parameter_file_string_current} not found..."
	initial_guesses = {}
	if (len(sys.argv)!=4):
		bf.print_cust("No parameter file given along command line. Searching current directory for parameter file.", quiet_val=quiet_val)
		if (os.path.isfile(parameter_file_string_current)):
			bf.print_cust(f"Parameter file found in the current directory. {str(parameter_file_string_current)}", quiet_val=quiet_val)
			parameter_file_name_final = parameter_file_string_current
		else:
			bf.print_cust("No parameter file (with default name - initial_parameters.dat) found in the current directory.", quiet_val=quiet_val)
			save_prompt= input("Would you like to provide the name of the parameter file? (y/n) : ")
			if (save_prompt=='y'):
				save_prompt2 = input("Please enter the name of the parameter file?: ")
				parameter_file_name_final = str(sys.argv[1]) + str("/") + str(save_prompt2)
			else:
				bf.print_cust(f"No parameter file given. Quitting...", quiet_val=quiet_val)
				quit()
	else:
		parameter_file_name_final = str(sys.argv[1]) + str("/") + str(sys.argv[3])
	bf.print_cust(f"Executing script with data from: {str(parameter_file_name_final)}", quiet_val=quiet_val)

	initial_guesses_rev = initial_guesses
	with open(str(parameter_file_name_final)) as f:
		for line in f:
			if '#' not in line:
				if (len(line.split())>2):
					(key, val) = line.split(':')
					key = key.replace(':', '').replace('-', '').lower()
					initial_guesses_rev[str(key)] = ast.literal_eval(val.replace(' ', ''))
				else:
					(key, val) = line.split()
					key = key.replace(':', '').replace('-', '').lower()
					initial_guesses_rev[str(key)] = val

	bf.print_cust(initial_guesses_rev, quiet_val=quiet_val)
	return (initial_guesses_rev)

########################################INITIAL_GUESS_FROM_PARAMETER_FILE########################################

######################REVISING_DICTIONARY##################################

def modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated):
	par_dict_updated = par_dict
	par_dict_updated['binning_quant'] = int(binning_quant_updated)
	par_dict_updated['binning_type'] = str(radio_binning_type_updated)
	par_dict_updated['voronoi_snr_type'] = str(radio_binning_significance_type_updated)
	par_dict_updated['execution_type'] = str(radio_work_type_updated)
	par_dict_updated['execution_fit_type'] = str(radio_fit_type_updated)
	par_dict_updated['emission_fit_type'] = str(radio_emission_fit_type_updated)
	return (par_dict_updated)

def revise_dictionary(par_dict_init, dir_name_6):
	default_sky = dir_name_6 + str('/default_sky.dat')
	default_lick_index_file = dir_name_6 + str('/default_lick_indexes.dat')
	default_emission_file = dir_name_6 + str('/default_emission_lines.dat')
	assert os.path.exists(str(default_sky)), f"File: {str(default_sky)} not found..."
	assert os.path.exists(str(default_lick_index_file)), f"File: {str(default_lick_index_file)} not found..."
	assert os.path.exists(str(default_emission_file)), f"File: {str(default_emission_file)} not found..."

	par_dict_emission = par_dict_init
	keys = ['number_of_narrow_components', 'number_of_wide_components', 'stopping_number_for_continuum', 'minimal_tying_factor_variance', 'minimum_amplitude_val', 'maximum_amplitude_val', 'minimum_reddening_val', 'maximum_reddening_val', 'sigma_bound_narrow_min', 'sigma_bound_narrow_max', 'sigma_bound_wide_max', 'poly_bound_val', 'maximum_accepted_reddening', 'e_b_minus_v_init', 'redshift_val', 'window_for_fit', 'multiplicative_factor_for_plot', 'spectral_smoothing', 'ppxf_stars_comp', 'ppxf_gas_comp', 'lick_idx_for_halpha', 'binning_quant', 'ellipse_angle_for_advanced_binning', 'window_for_choosing_snr_in_binning', 'max_ew_width_abs', 'allowed_ew_significance', 'decider_ew', 'num_processes', 'contour_count', 'image_allowed_percentile']
	values = [1, 0, 5000000.0, 1e-3, -1, 6, 1e-4, 2.0, 10.0, 200.0, 5000.0, 1000.0, 2.0, 0.0, 0.0, 2000.0, 10, 1, 1, 1, 6, 10, 0, 100, 100, 0.1, 0.0, 4, 15, 95]
	for key, value in zip(keys, values):
		if (key not in par_dict_emission):
			par_dict_emission[key] = float(value)
		else:
			par_dict_emission[key] = float(par_dict_emission[key])
	par_dict_emission['number_of_narrow_components'] = int(par_dict_emission['number_of_narrow_components'])
	par_dict_emission['number_of_wide_components'] = int(par_dict_emission['number_of_wide_components'])
	par_dict_emission['multiplicative_factor_for_plot'] = int(par_dict_emission['multiplicative_factor_for_plot'])
	par_dict_emission['spectral_smoothing'] = int(par_dict_emission['spectral_smoothing'])
	par_dict_emission['ppxf_stars_comp'] = int(par_dict_emission['ppxf_stars_comp'])
	par_dict_emission['ppxf_gas_comp'] = int(par_dict_emission['ppxf_gas_comp'])
	par_dict_emission['num_processes'] = int(par_dict_emission['num_processes'])
	keys2 = ['fit_velocity_val', 'fit_vel_disp_val', 'fit_continuum_val', 'fit_reddening_val', 'plot_fit_refined', 'ppxf_emission_tie_balmer', 'ppxf_emission_limit_doublets', 'quiet']
	values2 = [True, True, True, True, False, True, False, False]
	for key2, value2 in zip(keys2, values2):
		if (key2 not in par_dict_emission):
			par_dict_emission[key2] = bf.str2bool(value2)
		else:
			par_dict_emission[key2] = bf.str2bool(par_dict_emission[key2])
	center_init_array_tmp, sigma_init_array_tmp = bf.get_initial_kinematics_in_kms(int(par_dict_emission['number_of_narrow_components']), int(par_dict_emission['number_of_wide_components']))
	if 'center_init_array' not in par_dict_emission:
		par_dict_emission['center_init_array'] = center_init_array_tmp
	if 'sigma_init_array' not in par_dict_emission:
		par_dict_emission['sigma_init_array'] = sigma_init_array_tmp

	keys3 = ['observed_instrument', 'sky_spectra_file', 'lick_index_file', 'emission_file', 'binning_type', 'region_of_binning_interest_type', 'voronoi_snr_type', 'execution_type', 'execution_fit_type', 'emission_fit_type', 'continuum_finding_method', 'two_dimensional_map_name', ]
	values3 = ['VLT-MUSE', default_sky, default_lick_index_file, default_emission_file, 'None', 'logical', 'snr', 'snr_map', 'auto', 'custom', 'custom', 'default', ]
	for key3, value3 in zip(keys3, values3):
		if (key3 not in par_dict_emission):
			par_dict_emission[key3] = str(value3)
		else:
			par_dict_emission[key3] = str(par_dict_emission[key3])
	default_region_of_binning_interest = [0, 0, -1, -1]
	if ('region_of_binning_interest' not in par_dict_emission):
		par_dict_emission['region_of_binning_interest'] = list(default_region_of_binning_interest)
	default_binning_radius_list = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 500]
	if ('binning_radius_list' not in par_dict_emission):
		par_dict_emission['binning_radius_list'] = list(default_binning_radius_list)
	default_center_coordinate_for_advanced_binning = [0, 0]
	if ('center_coordinate_for_advanced_binning' not in par_dict_emission):
		par_dict_emission['center_coordinate_for_advanced_binning'] = list(default_center_coordinate_for_advanced_binning)
	default_st_age_list_start = [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]
	default_st_age_list_end = [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]
	if 'st_age_list_start' not in par_dict_emission:
		par_dict_emission['st_age_list_start'] = default_st_age_list_start
		par_dict_emission['st_age_unique_end'] = default_st_age_list_end
	else:
		if not ('st_age_unique_end' in par_dict_emission):
			par_dict_emission['st_age_unique_end'] = list(np.array(par_dict_emission['st_age_list_start'])*2)
		else:
			par_dict_emission['st_age_list_start'] = par_dict_emission['st_age_list_start']
			par_dict_emission['st_age_unique_end'] = par_dict_emission['st_age_unique_end']

	default_lick_index_species_list = [66, 70, 73, 74, 82, 84, 90, 100]
	if ('lick_index_species_list' not in par_dict_emission):
		par_dict_emission['lick_index_species_list'] = list(default_lick_index_species_list)

	default_emission_index_species_list = [15, 23, 24, 25, 26, 27]
	if ('emission_index_species_list' not in par_dict_emission):
		par_dict_emission['emission_index_species_list'] = list(default_emission_index_species_list)

	return (par_dict_emission)

def get_emission_keys(par_dict, wave_array):
	par_dict_rev = par_dict
	file_path = par_dict_rev['emission_file']
	#idx_new = [15, 23, 24, 25, 26, 27]
	idx_new = par_dict_rev['emission_index_species_list']
	data = pd.read_csv(file_path, sep="\t", comment='#', encoding=None)
	A = (data['Name'][idx_new].to_numpy().astype(np.str_))
	B = (data['Lambda'][idx_new].to_numpy().astype(np.int32))
	C = np.vstack([A,B])
	name_new = np.chararray([len(A)], itemsize=20)
	for i in range(len(A)):
		name_new[i] = str(A[i]) + str(B[i])
	name_new = name_new.decode("utf-8")
	wave_redshifted = (data['Lambda'][idx_new].to_numpy().astype(np.float32))*(1.+par_dict['redshift_val'])
	mask = (np.nanmin(wave_array)<wave_redshifted) & (wave_redshifted<np.nanmax(wave_array))
	par_dict_rev['center_list_names'] = name_new[mask]
	par_dict_rev['center_list'] = data['Lambda'][idx_new].to_numpy().astype(np.float32)[mask]
	par_dict_rev['comments_on_balmer'] = data['comments_on_balmer'][idx_new].to_numpy().astype(np.bool_)[mask]
	par_dict_rev['comments_on_tied'] = data['comments_on_tied'][idx_new].to_numpy().astype(np.str_)[mask]
	par_dict_rev['factor_for_tied'] = data['factor_for_tied'][idx_new].to_numpy().astype(np.float32)[mask]
	par_dict_rev['factor_fixed_for_tying'] = data['factor_fixed_for_tying'][idx_new].to_numpy().astype(np.bool_)[mask]
	return(par_dict_rev)


######################REVISING_DICTIONARY##################################



