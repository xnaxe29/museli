import os
import sys
import h5py
from progress.bar import IncrementalBar
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
import time
import ppxf as ppxf_package
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
import ppxf.miles_util as lib
from os import path
ppxf_dir = path.dirname(path.realpath(ppxf_package.__file__))
from astropy.constants import c as c_ms
c_kms = float(c_ms.value / 1e3)


quiet_val = False
#from basic_functions import *
import basic_functions as bf
bf.quite_val = quiet_val

#from plotting_functions import *
import plotting_functions as pf
pf.quite_val = quiet_val

#from get_data_from_file import *
import get_data_from_file as gdff
gdff.quite_val = quiet_val

#from handling_fits_files import *
import handling_fits_files as hff
hff.quite_val = quiet_val

#from spatial_binning_functions import *
import spatial_binning_functions as sbf
sbf.quite_val = quiet_val

#from emission_fit_function import *
import emission_fit_function as emff
emff.quite_val = quiet_val

#from extinction_fitting_functions import *
import extinction_fitting_functions as eff
eff.quite_val = quiet_val

#from ppxf_fitting_functions import *
import ppxf_fitting_functions as pff
pff.quite_val = quiet_val


########################################STAGE_ONE_ANALYSIS########################################

def stage_one_analysis(par_dict, new_filename, quiet_val=False):
	start_time1 = time.time()
	bf.print_cust('Performing STAGE-I analysis...', quiet_val=quiet_val)
	#LOAD ORIGINAL MUSE CUBE
	bf.print_cust('Loading Original Data...', quiet_val=quiet_val)
	header_original, header_original_err, data_original_file, err_original_file = hff.open_ifu_fits(str(par_dict['muse_data_filename']), str(par_dict['observed_instrument']))
	bf.print_cust('Original Data Loaded...', quiet_val=quiet_val)
	#OBTAINING PHYSICAL AXES
	bf.print_cust('Obtaining Physical Axes...', quiet_val=quiet_val)
	ra_array, dec_array, wave_array = hff.obtain_physical_axis(header_original)
	bf.print_cust('Physical Axes Retrieved...', quiet_val=quiet_val)
	#CLEANING DATA, REMOVING NANs, EMPTY ELEMENTS, etc.
	bf.print_cust('Cleaning Data...', quiet_val=quiet_val)
	data_original_cleaned_file = bf.clean_data(data_original_file)
	err_original_cleaned_file = bf.clean_data(err_original_file, type_of_data='err')
	bf.print_cust('Data Cleaned...', quiet_val=quiet_val)
	#Flattening data for convinience
	bf.print_cust('Flattening Original Data...', quiet_val=quiet_val)
	data_flattened = data_original_cleaned_file.reshape(data_original_cleaned_file.shape[0], -1)
	data_err_flattened = err_original_cleaned_file.reshape(err_original_cleaned_file.shape[0], -1)
	bf.print_cust('Original Data Flattened...', quiet_val=quiet_val)
	#Rebinning to logscale, spectral smoothing, sky and continuum estimation
	bf.print_cust('Rebinning to logscale, spectral smoothing, sky and continuum estimation...', quiet_val=quiet_val)
	wave_rebinned_central, data_rebinned_central, data_err_rebinned_central, fwhm_gal_central, velscale_central = gdff.refine_obs_data_using_scipy(wave_array, data_flattened[:, int(data_flattened.shape[1]/2)], data_err_flattened[:, int(data_flattened.shape[1]/2)], spectral_smoothing=int(par_dict['spectral_smoothing']))
	cont_init_central = gdff.get_initial_continuum_rev_2(wave_rebinned_central, data_rebinned_central)
	wave_sky_full_rebinned, flux_sky_full_rebinned, flux_err_sky_full_rebinned = gdff.check_for_sky(par_dict, wave_rebinned_central)
	data_rebinned = np.zeros([len(wave_rebinned_central), data_flattened.shape[1]])
	data_err_rebinned = np.ones([len(wave_rebinned_central), data_flattened.shape[1]])
	fwhm_gal_datacube = np.zeros([len(wave_rebinned_central), data_flattened.shape[1]])
	continuum_datacube = np.zeros([len(wave_rebinned_central), data_flattened.shape[1]])
	velscale_datacube = np.zeros([data_flattened.shape[1]])
	bar = IncrementalBar('Countdown', max = int(data_flattened.shape[1]))
	#bar = IncrementalBar('Countdown', max = int(data_flattened.shape[1]/100))
	for i in range(data_flattened.shape[1]):
	#for i in range(0, data_flattened.shape[1], 100):
		bar.next()
		data_rebinned[:, i], data_err_rebinned[:,i], fwhm_gal_datacube[:, i], velscale_datacube[i], continuum_datacube[:, i] = gdff.combined_function_rebinning_continuum(wave_array, data_flattened[:, i], data_err_flattened[:, i], redshift_val=float(par_dict['redshift_val']), quiet_val=quiet_val)
	bf.print_cust('Rebinning to logscale, spectral smoothing, sky and continuum estimation is now complete...', quiet_val=quiet_val)
	#Recovering new data back to the 3D shape
	bf.print_cust('Recovering new data back to the 3D shape...', quiet_val=quiet_val)
	data_rebinned_reshaped = data_rebinned.reshape((len(wave_rebinned_central), data_original_cleaned_file.shape[1], data_original_cleaned_file.shape[2]))
	data_err_rebinned_reshaped = data_err_flattened.reshape((len(wave_rebinned_central), data_original_cleaned_file.shape[1], data_original_cleaned_file.shape[2]))
	fwhm_gal_datacube_reshaped = fwhm_gal_datacube.reshape((len(wave_rebinned_central), data_original_cleaned_file.shape[1], data_original_cleaned_file.shape[2]))
	continuum_datacube_reshaped = continuum_datacube.reshape((len(wave_rebinned_central), data_original_cleaned_file.shape[1], data_original_cleaned_file.shape[2]))
	velscale_datacube_reshaped = velscale_datacube.reshape((data_original_cleaned_file.shape[1], data_original_cleaned_file.shape[2]))
	bf.print_cust('3D shape for new data recovered...', quiet_val=quiet_val)
	#Saving revised data
	bf.print_cust('Saving revised data...', quiet_val=quiet_val)
	hff.modify_ifu_fits(str(par_dict['muse_data_filename']), new_filename, wave_rebinned_central, data_rebinned_reshaped, data_err_cube=data_err_rebinned_reshaped, data_cont_cube=continuum_datacube_reshaped, data_fwhm_cube=fwhm_gal_datacube_reshaped, data_velscale_image=velscale_datacube_reshaped, sky_array=flux_sky_full_rebinned)
	bf.print_cust('Revised data saved...', quiet_val=quiet_val)
	bf.print_cust('STAGE-I analysis complete...', quiet_val=quiet_val)
	bf.print_cust(f"---STAGE-I took {float(time.time() - start_time1)} seconds ---", quiet_val=quiet_val)
	
########################################STAGE_ONE_ANALYSIS########################################








########################################STAGE_TWO_ANALYSIS########################################

def stage_two_analysis(par_dict, new_filename, file_name_eq_width, quiet_val=False):
	start_time2 = time.time()
	bf.print_cust('Performing STAGE-II analysis...', quiet_val=quiet_val)
	#To obtain information from revised data
	bf.print_cust('Loading information from revised data...', quiet_val=quiet_val)
	header_rev, header_rev_err, header_rev_cont, header_rev_fwhm, header_rev_velscale, header_rev_sky, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, velscale_rev_file, sky_rev_file = hff.open_ifu_fits_custom_file(new_filename)
	ra_rev, dec_rev, wave_rev = hff.obtain_physical_axis(header_rev)
	bf.print_cust('information from revised data loaded...', quiet_val=quiet_val)
	#bf.print_cust(f'{wave_rev.min()}, {wave_rev.max()}', quiet_val=quiet_val)
	bf.print_cust('Making EW Maps...', quiet_val=quiet_val)
	bf.print_cust('Flattening and Normalising Data...', quiet_val=quiet_val)
	data_flattened = data_rev_file.reshape(data_rev_file.shape[0], -1)
	data_err_flattened = data_rev_file.reshape(data_rev_file.shape[0], -1)
	data_cont_flattened = cont_rev_file.reshape(data_rev_file.shape[0], -1)
	data_flattened_normalised = data_flattened / data_cont_flattened
	data_err_flattened_normalised = data_err_flattened / data_cont_flattened
	bf.print_cust('Data Flattened and Normalised...', quiet_val=quiet_val)
	wave_rest = wave_rev / (1.+float(par_dict['redshift_val']))
	assert os.path.exists(str(par_dict['lick_index_file'])), f"File: {str(par_dict['lick_index_file'])} not found..."
	if ("lick_index_species_list" in par_dict):
		lick_index_species_list = par_dict['lick_index_species_list']
	else:
		lick_index_species_list = [0, 1]

	lick_index_wave11, lick_index_wave12, lick_index_wave21, lick_index_wave22, lick_index_wave31, lick_index_wave32, lick_index_sign, lick_index_species, lick_index_reference = np.genfromtxt(str(par_dict['lick_index_file']), unpack=True, names=True, encoding=None, dtype=None)
	mean_wave_list_full = (lick_index_wave21 + lick_index_wave22) / 2.
	lick_index_species_rev = np.array(lick_index_species[lick_index_species_list[:]])
	mean_wave_list_rev = np.array(mean_wave_list_full[lick_index_species_list[:]])
	bf.print_cust(f"Estimating EW for: {lick_index_species_rev}", quiet_val=quiet_val)
	ew_map_array = np.zeros([len(mean_wave_list_rev), data_flattened_normalised.shape[1], 2])
	bar = IncrementalBar('Countdown', max = int(data_flattened_normalised.shape[1]))
	#bar = IncrementalBar('Countdown', max = int(data_flattened_normalised.shape[1]/100))
	for i in range(data_flattened_normalised.shape[1]):
	#for i in range(0, data_flattened_normalised.shape[1], 100):
		bar.next()
		for j in range(len(mean_wave_list_rev)):
			ew_map_array[j, i, 0], ew_map_array[j, i, 1] = bf.get_equivalendth_width_rev(wave_rest, data_flattened_normalised[:, i], data_err_flattened_normalised[:, i], center=float(mean_wave_list_rev[j]), redshift=float(par_dict['redshift_val']))

	ew_map_array_cube = ew_map_array.reshape((len(mean_wave_list_rev), data_rev_file.shape[1], data_rev_file.shape[2], 2))
	#file_name_eq_width = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_eq_width.fits')
	header_ew, header_ew_err, data_ew_file, err_ew_file = hff.copy_ifu_fits(str(par_dict['muse_data_filename']), file_name_eq_width)
	with fits.open(file_name_eq_width, mode='update', output_verify='fix') as hdu_list:
		hdu_list[1].data = ew_map_array_cube[:,:,:,0]
		hdu_list[2].data = ew_map_array_cube[:,:,:,1]
		hdu_list[1].header['NAXIS1'] = (len(mean_wave_list_rev), 'length of data axis 1 ')
		hdu_list[2].header['NAXIS1'] = (len(mean_wave_list_rev), 'length of data axis 1 ')
		hdu_list[1].header['BUNIT'] = 'Ang'
		hdu_list[2].header['BUNIT'] = 'Ang'
		#del hdu_list[1].header['CTYPE3']; del hdu_list[1].header['CUNIT3']; del hdu_list[1].header['CD3_3']; del hdu_list[1].header['CRPIX3']; del hdu_list[1].header['CRVAL3']; del hdu_list[1].header['CD1_3']; del hdu_list[1].header['CD2_3']; del hdu_list[1].header['CD3_1']; del hdu_list[1].header['CD3_2']; del hdu_list[2].header['CTYPE3']; del hdu_list[2].header['CUNIT3']; del hdu_list[2].header['CD3_3']; del hdu_list[2].header['CRPIX3']; del hdu_list[2].header['CRVAL3']; del hdu_list[2].header['CD1_3']; del hdu_list[2].header['CD2_3']; del hdu_list[2].header['CD3_1']; del hdu_list[2].header['CD3_2']
	bf.print_cust('EW Maps generated...', quiet_val=quiet_val)
	bf.print_cust('STAGE-II analysis complete...', quiet_val=quiet_val)
	bf.print_cust(f"---STAGE-II took {float(time.time() - start_time2)} seconds ---", quiet_val=quiet_val)
	
########################################STAGE_TWO_ANALYSIS########################################







########################################STAGE_THREE_ANALYSIS########################################

def stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict, display_binning=False, quiet_val=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	start_time3 = time.time()
	if ('region_of_binning_interest' in par_dict):
		region_of_binning_interest = list(par_dict['region_of_binning_interest'])
	if ('region_of_binning_interest_type' in par_dict):
		region_of_binning_interest_type = str(par_dict['region_of_binning_interest_type'])

	if 'voronoi' in str(par_dict['binning_type'].lower()):
		assert ('binning_quant' in par_dict), f"Binning quant required..."
		assert ('voronoi_snr_type' in par_dict), f"voronoi_snr_type required..."
		file_name_rev_binned = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_snr_type_') + str(par_dict['voronoi_snr_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('.dat')
		bf.print_cust(file_name_rev_binned)
		sbf.voronoi_binning_func(file_name_rev_linear, file_name_rev_binned, int(par_dict['binning_quant']), str(par_dict['voronoi_snr_type']), quiet_val=False, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
	elif 'advanced' in str(par_dict['binning_type'].lower()):
		assert ('binning_quant' in par_dict), f"Binning quant required..."
		assert ('ellipse_angle_for_advanced_binning' in par_dict), f"ellipse_angle_for_advanced_binning required..."
		file_name_rev_binned = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_ellipse_type_') + str(par_dict['ellipse_angle_for_advanced_binning']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('.dat')
		sbf.advanced_binning_func(file_name_rev_linear, file_name_rev_binned, par_dict, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
	elif 'square' in str(par_dict['binning_type'].lower()):
		assert ('binning_quant' in par_dict), f"Binning quant required..."
		file_name_rev_binned = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('.dat')
		sbf.square_binning_func(file_name_rev_linear, file_name_rev_binned, par_dict, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
	elif 'none' in str(par_dict['binning_type'].lower()):
		bf.print_cust('Spatial Binning not required by user', quiet_val=quiet_val)
		file_name_rev_binned = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_no_spatial_binning.dat')
		sbf.no_spatial_binning_func(file_name_rev_linear, file_name_rev_binned, par_dict)
	else:
		bf.print_cust('Spatial Binning not specified by user. Setting to no spatial binning.', quiet_val=quiet_val)
		file_name_rev_binned = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_no_spatial_binning.dat')
		sbf.no_spatial_binning_func(file_name_rev_linear, file_name_rev_binned, par_dict)
	if (display_binning):
		sbf.binning_result_display(file_name_rev_linear, file_name_rev_binned, par_dict)
	bf.print_cust(f"---STAGE-III took {float(time.time() - start_time3)} seconds ---", quiet_val=quiet_val)
	return(file_name_rev_binned)
	
########################################STAGE_THREE_ANALYSIS########################################






########################################STAGE_FOUR_ANALYSIS########################################

########################################STAGE_FOUR_EXECUTION########################################
def stage_four_execution(par_dict, file_name_rev_linear, file_name_rev_binned, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, quiet_val=False, plot_val=True):
	start_time4 = time.time()
	bf.print_cust("Executing stage four...", quiet_val=quiet_val)
	if 'snr_map' in str(par_dict['execution_type'].lower()):
		if (plot_val):
			bf.print_cust('Getting SNR Map...', quiet_val=quiet_val)
			pf.get_snr_map_revised(file_name_rev_linear, file_name_rev_binned, par_dict)
	if 'fit' in str(par_dict['execution_type'].lower()):
		bf.print_cust('Fitting...', quiet_val=quiet_val)
		x_axis_binned_rev, y_axis_binned_rev, bin_num = np.loadtxt(file_name_rev_binned).T
		bin_num_unique = np.unique(bin_num)
		main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = get_output_filenames_and_arrays(dir_name_4, bin_num_unique, wave_rev, par_dict)
		#bf.print_cust(main_result_cube_filename)
		if (os.path.exists(main_result_cube_filename)):
			prompt_4=input(f"File {main_result_cube_filename} already exists. Overwrite?(y/n) : ")
			if ('y' in prompt_4.lower()):
				bf.print_cust(f"deleting...{main_result_cube_filename}", quiet_val=quiet_val)
				string_for_deleting = str('rm ') + main_result_cube_filename
				bf.print_cust(string_for_deleting)
				os.system(string_for_deleting)
				bf.print_cust(f"deleting...{main_absorption_result_cube_filename}", quiet_val=quiet_val)
				string_for_deleting = str('rm ') + main_absorption_result_cube_filename
				bf.print_cust(string_for_deleting)
				os.system(string_for_deleting)
				bf.print_cust(f"deleting...{main_emission_result_cube_filename}", quiet_val=quiet_val)
				string_for_deleting = str('rm ') + main_emission_result_cube_filename
				bf.print_cust(string_for_deleting)
				os.system(string_for_deleting)
				fitting_regime(file_name_rev_binned, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, main_result_cube, main_absorption_result_cube, main_emission_result_cube, main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict)
		else:
			fitting_regime(file_name_rev_binned, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, main_result_cube, main_absorption_result_cube, main_emission_result_cube, main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict)
	bf.print_cust("Stage four execution complete...", quiet_val=quiet_val)
	bf.print_cust(f"---STAGE-IV took {float(time.time() - start_time4)} seconds ---", quiet_val=quiet_val)
	#return None
	#continue
########################################STAGE_FOUR_EXECUTION########################################


########################################FITTING_REGIME_FUNCTION########################################
def fitting_regime(file_name_rev_binned, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, main_result_cube, main_absorption_result_cube, main_emission_result_cube, main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, quiet_val=False):
	bf.print_cust("Starting fitting...", quiet_val=quiet_val)
	x_axis_binned_rev, y_axis_binned_rev, bin_num = np.loadtxt(file_name_rev_binned).T
	bin_num_unique = np.unique(bin_num)
	x_axis_binned_rev = x_axis_binned_rev-1
	y_axis_binned_rev = y_axis_binned_rev-1
	#wave_rev = wave_rev / (1.+float(par_dict['redshift_val']))
	halpha_eq_width_map[np.abs(halpha_eq_width_map)>float(par_dict['max_ew_width_abs'])]=np.nan
	#halpha_eq_width_err_map[np.abs(halpha_eq_width_map)>100.]=np.nan
	flux_sky1 = sky_rev_file
	temporary_data_cube = np.zeros([len(bin_num_unique), len(wave_rev), 3])
	temporary_information_cube = np.zeros([len(bin_num_unique), 2])
	bar = IncrementalBar('Countdown', max = len(bin_num_unique))
	for i in range(len(bin_num_unique)):
		bar.next()
		mask = np.in1d(bin_num, [bin_num_unique[i]])
		flux_tmp = data_rev_file[:, y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)]
		err_tmp = err_rev_file[:, y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)]
		cont_tmp = cont_rev_file[:, y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)]
		flux1 = np.nansum(flux_tmp, axis=(1))/len(x_axis_binned_rev[mask])
		err1 = np.sqrt(np.nansum(err_tmp, axis=(1)))/len(x_axis_binned_rev[mask])
		#halpha_ew_val1 = np.nanmedian(halpha_eq_width_map[y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)])
		halpha_ew_val1 = np.nanmedian(halpha_eq_width_map[y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)])
		FWHM_gal1 = np.nanmedian(fwhm_rev_file[:, y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)])
		init_cont1 = np.nansum(cont_tmp, axis=(1))/len(x_axis_binned_rev[mask])
		wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = main_fitting_func(wave_rev, flux1, err1, par_dict, halpha_ew_val1, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1)
		if ('emission' in par_dict['execution_fit_type']):
			main_emission_result_cube[i,:,0] = solution_fit
			main_emission_result_cube[i,:,1] = solution_err_fit
		elif ('absorption' in par_dict['execution_fit_type']):
			main_absorption_result_cube[i,:,0] = solution_fit
			main_absorption_result_cube[i,:,1] = solution_err_fit
		else:
			if (halpha_ew_val1>=float(par_dict['decider_ew'])):
				main_absorption_result_cube[i,:,0] = solution_fit
				main_absorption_result_cube[i,:,1] = solution_err_fit
			else:
				main_emission_result_cube[i,:,0] = solution_fit
				main_emission_result_cube[i,:,1] = solution_err_fit
		main_result_cube[i,0:len(wave_fit),0] = wave_fit
		main_result_cube[i,0:len(wave_fit),1] = flux_fit
		main_result_cube[i,0:len(wave_fit),2] = flux_err_fit
		main_result_cube[i,0:len(wave_fit),3] = continuum_fit
		main_result_cube[i,0:len(wave_fit),4] = result_fit
		main_result_cube[i,0:len(wave_fit),5] = mask_fit
	bf.print_cust("Saving files...", quiet_val=quiet_val)
	np.save(main_result_cube_filename, main_result_cube)
	bf.print_cust(f"File: {main_result_cube_filename} saved...", quiet_val=quiet_val)
	np.save(main_absorption_result_cube_filename, main_absorption_result_cube)
	bf.print_cust(f"File: {main_absorption_result_cube_filename} saved...", quiet_val=quiet_val)
	np.save(main_emission_result_cube_filename, main_emission_result_cube)
	bf.print_cust(f"File: {main_emission_result_cube_filename} saved...", quiet_val=quiet_val)
########################################FITTING_REGIME_FUNCTION########################################

##################################STELLAR_MAIN_FITTING_FUNCTION##################################
def main_fitting_func(wave_rev, flux, err, par_dict, halpha_ew_val, **kwargs):
	FWHM_gal1 = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	flux_sky1 = kwargs.get('flux_sky', np.zeros_like(wave_rev))  # Sky
	init_cont1 = kwargs.get('init_cont', np.ones_like(wave_rev))  # continuum
	quiet_run = kwargs.get('quiet_run', True)  # quiet_run
	if ('emission' in par_dict['execution_fit_type']):
		wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = emission_fit(wave_rev, flux, err, par_dict, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1, quiet_val=quiet_run)
	elif ('absorption' in par_dict['execution_fit_type']):
		wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = absorption_fit(wave_rev, flux, err, par_dict, flux_sky = flux_sky1, quiet_val=quiet_run)
	else:
		if (halpha_ew_val>=float(par_dict['decider_ew'])):
			wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = absorption_fit(wave_rev, flux, err, par_dict, flux_sky = flux_sky1, quiet_val=quiet_run)
		else:
			wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = emission_fit(wave_rev, flux, err, par_dict, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1, quiet_val=quiet_run)
	return (wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit)

#wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = main_fitting_func(wave_rev, flux, err, par_dict, halpha_ew_val, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1)
##################################STELLAR_MAIN_FITTING_FUNCTION##################################

##################################STELLAR_ABSORPTION_LINE_FITTING##################################
def absorption_fit(wave_rev, flux, err, par_dict, **kwargs):
	flux_sky1 = kwargs.get('flux_sky', np.zeros_like(wave_rev))  # Sky
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	if (par_dict['ppxf_stars_comp']<1):
		par_dict['ppxf_stars_comp']=1
	file1 = np.array([wave_rev, flux, err])
	lam_gal_rev, fit_galaxy, fit_noise, fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, fit_normalising_factor = pff.ppxf_function(file1, redshift=par_dict['redshift_val'], number_of_stellar_components=par_dict['ppxf_stars_comp'], number_of_gas_components=par_dict['ppxf_gas_comp'], sky_flux_init=flux_sky1, first_order_kinematic_limits=1000., tie_balmer=par_dict['ppxf_emission_tie_balmer'], limit_doublets=par_dict['ppxf_emission_limit_doublets'], quiet=quiet_val)
	bestfit_solution_array, stellar_spectrum_array, sol_pp_gas_flux_reshaped, sol_pp_gas_flux_err_reshaped, stellar_vel, stellar_vel_err, stellar_sigma, stellar_sigma_err, stellar_h3, stellar_h3_err, stellar_h4, stellar_h4_err, gas_balmer_vel, gas_balmer_vel_err, gas_balmer_sigma, gas_balmer_sigma_err, gas_others_vel, gas_others_vel_err, gas_others_sigma, gas_others_sigma_err, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, goodpixels_idx, reduced_chi_sq, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol, str_fit = pff.get_detailed_ppxf_solution(fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, par_dict['ppxf_stars_comp'], par_dict['ppxf_gas_comp'], normalising_factor = fit_normalising_factor, st_age_unique_start=par_dict['st_age_list_start'], st_age_unique_end=par_dict['st_age_unique_end'], quiet=quiet_val)
	masked_array = np.full([len(bestfit_solution_array)], True)
	for i in range(len(goodpixels_idx)):
		masked_array[goodpixels_idx[i]] = True
	solution_array = np.append(stellar_vel, stellar_sigma)
	solution_array = np.append(solution_array, stellar_h3)
	solution_array = np.append(solution_array, stellar_h4)
	solution_array = np.append(solution_array, reddening_fit_e_b_minus_v)
	solution_array = np.append(solution_array, gas_reddening_fit_e_b_minus_v)
	solution_array = np.append(solution_array, gas_balmer_vel)
	solution_array = np.append(solution_array, gas_balmer_sigma)
	solution_array = np.append(solution_array, gas_others_vel)
	solution_array = np.append(solution_array, gas_others_sigma)
	solution_array = np.append(solution_array, st_age_unique)
	solution_array = np.append(solution_array, st_mass_unique)
	solution_array = np.append(solution_array, st_lum_unique)
	solution_array = np.append(solution_array, weighted_logAge)
	solution_array = np.append(solution_array, weighted_metallicity)
	solution_err_array = np.append(stellar_vel_err, stellar_sigma_err)
	solution_err_array = np.append(solution_err_array, stellar_h3_err)
	solution_err_array = np.append(solution_err_array, stellar_h4_err)
	solution_err_array = np.append(solution_err_array, status_of_optimization)
	solution_err_array = np.append(solution_err_array, reduced_chi_sq)
	solution_err_array = np.append(solution_err_array, gas_balmer_vel_err)
	solution_err_array = np.append(solution_err_array, gas_balmer_sigma_err)
	solution_err_array = np.append(solution_err_array, gas_others_vel_err)
	solution_err_array = np.append(solution_err_array, gas_others_sigma_err)
	solution_err_array = np.append(solution_err_array, st_age_unique)
	solution_err_array = np.append(solution_err_array, st_mass_unique)
	solution_err_array = np.append(solution_err_array, st_lum_unique)
	solution_err_array = np.append(solution_err_array, mlpop_sol)
	solution_err_array = np.append(solution_err_array, mlpop_sol)
	return (lam_gal_rev, fit_galaxy, fit_noise, stellar_spectrum_array, bestfit_solution_array, masked_array, solution_array, solution_err_array)
##################################STELLAR_ABSORPTION_LINE_FITTING##################################


##################################GAS_EMISSION_LINE_FITTING##################################
def emission_fit(wave_rev, flux, err, par_dict, **kwargs):
	FWHM_gal1 = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	flux_sky1 = kwargs.get('flux_sky', np.zeros_like(wave_rev))  # Sky
	init_cont1 = kwargs.get('init_cont', np.ones_like(wave_rev))  # continuum
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	if ('custom'in par_dict['emission_fit_type']):
		wave_fit, flux_fit, err_fit, continuum_fitted, res_fitted_full_fitted, masked_array, pfit_curvefit, perr_curvefit, line_flux, line_flux_err, red_chi_squared, amp_length = emff.complex_fit_func(wave_rev, flux, err, par_dict, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, cont_init = init_cont1, quiet_val=quiet_val)
		solution_array = np.append(pfit_curvefit, line_flux)
		solution_array = np.append(solution_array, red_chi_squared)
		solution_err_array = np.append(perr_curvefit, line_flux_err)
		solution_err_array = np.append(solution_err_array, amp_length)
		return (wave_fit, flux_fit, err_fit, continuum_fitted, res_fitted_full_fitted, masked_array, solution_array, solution_err_array)
	elif ('ppxf'in par_dict['emission_fit_type']):
		if (par_dict['ppxf_gas_comp']<1):
			par_dict['ppxf_gas_comp']=1
		file1 = np.array([wave_rev, flux, err])
		lam_gal_rev, fit_galaxy, fit_noise, fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, fit_normalising_factor = pff.ppxf_function(file1, redshift=par_dict['redshift_val'], number_of_stellar_components=par_dict['ppxf_stars_comp'], number_of_gas_components=par_dict['ppxf_gas_comp'], sky_flux_init=flux_sky1, first_order_kinematic_limits=1000., tie_balmer=par_dict['ppxf_emission_tie_balmer'], limit_doublets=par_dict['ppxf_emission_limit_doublets'], quiet=quiet_val)
		bestfit_solution_array, stellar_spectrum_array, sol_pp_gas_flux_reshaped, sol_pp_gas_flux_err_reshaped, stellar_vel, stellar_vel_err, stellar_sigma, stellar_sigma_err, stellar_h3, stellar_h3_err, stellar_h4, stellar_h4_err, gas_balmer_vel, gas_balmer_vel_err, gas_balmer_sigma, gas_balmer_sigma_err, gas_others_vel, gas_others_vel_err, gas_others_sigma, gas_others_sigma_err, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, goodpixels_idx, reduced_chi_sq, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol, str_fit = pff.get_detailed_ppxf_solution(fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, par_dict['ppxf_stars_comp'], par_dict['ppxf_gas_comp'], normalising_factor = fit_normalising_factor, st_age_unique_start=par_dict['st_age_list_start'], st_age_unique_end=par_dict['st_age_unique_end'], quiet=quiet_val)
		masked_array = np.full([len(bestfit_solution_array)], True)
		for i in range(len(goodpixels_idx)):
			masked_array[goodpixels_idx[i]] = True
		
		'''
		revised_line_flux = np.zeros([sol_pp_gas_flux_reshaped.shape[1]])
		revised_line_flux_err = np.zeros([sol_pp_gas_flux_reshaped.shape[1]])
		for i in range(sol_pp_gas_flux_reshaped.shape[1]):
			for j in range(sol_pp_gas_flux_reshaped.shape[0]):
				revised_line_flux[i] += sol_pp_gas_flux_reshaped[j][i]
				revised_line_flux_err[i] += sol_pp_gas_flux_err_reshaped[j][i]**2
		revised_line_flux_err = np.sqrt(revised_line_flux_err)
		'''
		revised_line_flux = sol_pp_gas_flux_reshaped
		revised_line_flux_err = sol_pp_gas_flux_err_reshaped
		solution_array = np.append(gas_balmer_vel, gas_balmer_sigma)
		solution_array = np.append(solution_array, gas_others_vel)
		solution_array = np.append(solution_array, gas_others_sigma)
		solution_array = np.append(solution_array, revised_line_flux)
		solution_array = np.append(solution_array, reduced_chi_sq)
		solution_err_array = np.append(gas_balmer_vel_err, gas_balmer_sigma_err)
		solution_err_array = np.append(solution_err_array, gas_others_vel_err)
		solution_err_array = np.append(solution_err_array, gas_others_sigma_err)
		solution_err_array = np.append(solution_err_array, revised_line_flux_err)
		solution_err_array = np.append(solution_err_array, gas_reddening_fit_e_b_minus_v)
		return (lam_gal_rev, fit_galaxy, fit_noise, stellar_spectrum_array, bestfit_solution_array, masked_array, solution_array, solution_err_array)

#wave_em_fit, flux_em_fit, flux_err_em_fit, continuum_em_fit, result_em_fit, mask_em_fit, solution_em_fit, solution_err_em_fit = emission_fit(wave_rev, flux, err, par_dict, **kwargs)
##################################GAS_EMISSION_LINE_FITTING##################################

##################################GET_OUTPUT_FILENAMES_AND_ARRAYS##################################
def get_output_filenames_and_arrays(dir_name_4, bin_num_unique, wave_rev, par_dict, quiet_val=True):
	par_dict = gdff.get_emission_keys(par_dict, wave_rev)
	main_result_cube_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_snr_type_') + str(par_dict['voronoi_snr_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('_results_cube') + str('.npy')
	main_result_cube = np.zeros([len(bin_num_unique), len(wave_rev), 6])
	main_absorption_result_cube_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_snr_type_') + str(par_dict['voronoi_snr_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('_custom_absorption_results') + str('.npy')
	total_count_of_abs_pars = (par_dict['ppxf_stars_comp']*4) + 2 + (par_dict['ppxf_gas_comp']*4) + (len(par_dict['st_age_list_start'])*3) + 2
	main_absorption_result_cube = np.zeros([len(bin_num_unique), total_count_of_abs_pars, 2])
	if ('custom'in par_dict['emission_fit_type']):
		flux = np.ones_like(wave_rev)
		init_cont = np.ones_like(wave_rev)
		popt_init, masked_array, amp_init_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, coeff_init_val, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp = emff.get_opt_param_array_rev2(wave_rev, flux, init_cont, par_dict['number_of_narrow_components'], par_dict['number_of_wide_components'], par_dict['center_list'], par_dict['redshift_val'], par_dict['comments_on_balmer'], par_dict['comments_on_tied'], par_dict['factor_for_tied'], par_dict['e_b_minus_v_init'], par_dict['window_for_fit'], par_dict['stopping_number_for_continuum'], par_dict['center_init_array'], par_dict['sigma_init_array'])
		center_list_len = len(par_dict['center_list'])
		bf.print_cust(f'{len(popt_init)}, {center_list_len}', quiet_val=quiet_val)
		main_emission_result_cube_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_snr_type_') + str(par_dict['voronoi_snr_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('_emission_type_') + str(par_dict['emission_fit_type'])  + str('_custom_emission_results') + str('.npy')
		main_emission_result_cube = np.zeros([len(bin_num_unique), len(popt_init)+len(coeff_init_val)+len(par_dict['center_list'])+1, 2])
	elif ('ppxf'in par_dict['emission_fit_type']):
		if (par_dict['ppxf_gas_comp']<1):
			par_dict['ppxf_gas_comp']=1
		pathname_ppxf = ppxf_dir + '/miles_models/*.fits'  # Initialize a default stellar template
		miles = lib.miles(pathname_ppxf, 3.0, 55.0)
		lam_range_gal = np.array([np.min(wave_rev), np.max(wave_rev)])/(1. + par_dict['redshift_val'])
		bf.blockPrint()
		gas_templates, gas_names, line_wave = util.emission_lines(miles.ln_lam_temp, lam_range_gal, 3.0, tie_balmer=par_dict['ppxf_emission_tie_balmer'], limit_doublets=par_dict['ppxf_emission_limit_doublets'])
		bf.enablePrint()
		main_emission_result_cube_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_binned_') + str(par_dict['binning_type']) + str('_snr_type_') + str(par_dict['voronoi_snr_type']) + str('_quant_') + str(int(par_dict['binning_quant'])) + str('_emission_type_') + str(par_dict['emission_fit_type']) + str('_custom_emission_results') + str('.npy')
		main_emission_result_cube = np.zeros([len(bin_num_unique), int(((4+len(gas_names))*par_dict['ppxf_gas_comp'])+1), 2])
	return (main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube)

#main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = get_output_filenames_and_arrays(dir_name_4, bin_num_unique, wave_rev, par_dict)
##################################GET_OUTPUT_FILENAMES_AND_ARRAYS##################################

########################################STAGE_FOUR_ANALYSIS########################################





########################################STAGE_FIVE_ANALYSIS########################################

########################################MAIN_FUNCTION########################################

def stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_linear, file_name_rev_binned, expanded_hdr_filename, header_rev, expanded_filename, quiet_val=False):
	start_time7 = time.time()
	bf.print_cust("Executing stage five...", quiet_val=quiet_val)
	if (os.path.exists(expanded_hdr_filename) or os.path.exists(expanded_filename)):
		prompt_5=input(f"File {expanded_hdr_filename} or {expanded_filename} already exists. Repeat expansion?(y/n) : ")
		if ('y' in prompt_5.lower()):
			bf.print_cust(f"deleting...{expanded_hdr_filename} and {expanded_filename}")
			string_for_deleting = str('rm ') + expanded_hdr_filename
			bf.print_cust(string_for_deleting)
			os.system(string_for_deleting)
			string_for_deleting2 = str('rm ') + expanded_filename
			bf.print_cust(string_for_deleting2)
			os.system(string_for_deleting2)
			expanded_hdr_filename1, expanded_filename1 = expand_results(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_linear, file_name_rev_binned, header_rev, expanded_hdr_filename, expanded_filename)
		else:
			expanded_hdr_filename1 = expanded_hdr_filename
			expanded_filename1 = expanded_filename
	else:
		expanded_hdr_filename1, expanded_filename1 = expand_results(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_linear, file_name_rev_binned, header_rev, expanded_hdr_filename, expanded_filename)

	bf.print_cust("Stage five execution complete...", quiet_val=quiet_val)
	bf.print_cust(f"---STAGE-V took {float(time.time() - start_time7)} seconds ---", quiet_val=quiet_val)
	return(expanded_hdr_filename1, expanded_filename1)

def expand_results(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_linear, file_name_rev_binned, header_rev, expanded_hdr_filename, expanded_filename):
	assert os.path.exists(str(main_result_cube_filename)), f"File: {str(main_result_cube_filename)} not found..."
	assert os.path.exists(str(main_absorption_result_cube_filename)), f"File: {str(main_absorption_result_cube_filename)} not found..."
	assert os.path.exists(str(main_emission_result_cube_filename)), f"File: {str(main_emission_result_cube_filename)} not found..."
	bf.print_cust("Expanding Results...")
	start_time5 = time.time()
	main_result_cube_init = np.load(main_result_cube_filename)
	main_absorption_result_init = np.load(main_absorption_result_cube_filename)
	main_emission_result_init = np.load(main_emission_result_cube_filename)
	y_axis_unbinned, x_axis_unbinned, ra_unbinned, dec_unbinned, signal11, noise11, signal12, noise12 = np.loadtxt(file_name_rev_linear).T
	#ra_unbinned_unique = np.unique(ra_unbinned)
	#dec_unbinned_unique = np.unique(dec_unbinned)
	x_axis_binned_rev, y_axis_binned_rev, bin_num = np.loadtxt(file_name_rev_binned).T
	bin_num_unique = np.unique(bin_num)
	x_axis_binned_rev = x_axis_binned_rev-1
	y_axis_binned_rev = y_axis_binned_rev-1
	x_axis_unique = np.unique(y_axis_binned_rev)
	y_axis_unique = np.unique(x_axis_binned_rev)
	expanded_data_cube = np.zeros([len(wave_rev), len(x_axis_unique), len(y_axis_unique), 5])
	expanded_information_cube_absorption = np.zeros([main_absorption_result_init.shape[1], len(x_axis_unique), len(y_axis_unique), 2])
	header_information_absorption, header_information_err_absorption = get_expanded_information_abs(par_dict)
	header_information_emission, header_information_err_emission, amp_length = get_expanded_header_information_emission(wave_rev, par_dict)
	#print (header_information_emission)
	expanded_information_cube_emission = np.zeros([main_emission_result_init.shape[1], len(x_axis_unique), len(y_axis_unique), 2])
	bar = IncrementalBar('Countdown', max = len(bin_num_unique))
	wave_new = main_result_cube_init[0, :, 0]
	for i in range(len(bin_num_unique)):
		bar.next()
		mask = np.in1d(bin_num, [bin_num_unique[i]])
		#wave_new = main_result_cube_init[i, :, 0]
		flux_new = main_result_cube_init[i, :, 1]
		flux_err_new = main_result_cube_init[i, :, 2]
		cont_new = main_result_cube_init[i, :, 3]
		result_new = main_result_cube_init[i, :, 4]
		mask_new = main_result_cube_init[i, :, 5]
		x_indexes = x_axis_binned_rev[mask.astype(bool)].astype(int)
		y_indexes = y_axis_binned_rev[mask.astype(bool)].astype(int)
		for j in range(len(y_indexes)):
			#expanded_data_cube[:, y_indexes[j], x_indexes[j], 0] = wave_new
			expanded_data_cube[:, y_indexes[j], x_indexes[j], 0] = flux_new
			expanded_data_cube[:, y_indexes[j], x_indexes[j], 1] = flux_err_new
			expanded_data_cube[:, y_indexes[j], x_indexes[j], 2] = cont_new
			expanded_data_cube[:, y_indexes[j], x_indexes[j], 3] = result_new
			expanded_data_cube[:, y_indexes[j], x_indexes[j], 4] = mask_new
			expanded_information_cube_absorption[:, y_indexes[j], x_indexes[j], 0] = main_absorption_result_init[i,:,0]
			expanded_information_cube_absorption[:, y_indexes[j], x_indexes[j], 1] = main_absorption_result_init[i,:,1]
			expanded_information_cube_emission[:, y_indexes[j], x_indexes[j], 0] = main_emission_result_init[i,:,0]
			expanded_information_cube_emission[:, y_indexes[j], x_indexes[j], 1] = main_emission_result_init[i,:,1]

	if (os.path.exists(expanded_filename)):
		rem_str = str('rm ') + str(expanded_filename)
		os.system(rem_str)
	if (os.path.exists(expanded_hdr_filename)):
		rem_str = str('rm ') + str(expanded_hdr_filename)
		os.system(rem_str)

	'''
	expanded_hdr_data = []
	expanded_hdr_data.extend([header_information_absorption])
	expanded_hdr_data.extend([header_information_err_absorption])
	expanded_hdr_data.extend([header_information_emission])
	expanded_hdr_data.extend([header_information_err_emission])
	np.save(expanded_hdr_filename, expanded_hdr_data)
	'''

	primary_hdu = fits.PrimaryHDU(header=header_rev)
	c1 = fits.Column(name='header_information_absorption', array=np.array(header_information_absorption), format='A20V')
	c2 = fits.Column(name='header_information_err_absorption', array=np.array(header_information_err_absorption), format='A20V')
	c3 = fits.Column(name='header_information_emission', array=np.array(header_information_emission), format='A20V')
	c4 = fits.Column(name='header_information_err_emission', array=np.array(header_information_err_emission), format='A20V')
	table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
	hdul = fits.HDUList([primary_hdu, table_hdu])
	hdul.writeto(expanded_hdr_filename, overwrite=True)
	ra_unbinned_unique, dec_unbinned_unique, wave_unbinned_unique = hff.obtain_physical_axis(header_rev)
	with h5py.File(expanded_filename, 'w') as hf:
		hf.create_dataset("amp_length",  data=amp_length)
		hf.create_dataset("wave_new",  data=wave_new)
		hf.create_dataset("ra",  data=ra_unbinned_unique)
		hf.create_dataset("dec",  data=dec_unbinned_unique)
		hf.create_dataset("expanded_data_cube",  data=expanded_data_cube)
		hf.create_dataset("expanded_information_cube_absorption",  data=expanded_information_cube_absorption)
		hf.create_dataset("expanded_information_cube_emission",  data=expanded_information_cube_emission)
	bf.print_cust("Expansion successful")
	bf.print_cust(f"Expanding data took {float(time.time() - start_time5)} seconds ---")
	return(expanded_hdr_filename, expanded_filename)

########################################MAIN_FUNCTION########################################

###############################GET_EXPANDED_DATA_INFORMATION_EMISSION###############################

def get_expanded_data_information_emission(popt_fit, perr_fit, wave, par_dict, **kwargs):
	par_dict_rev = par_dict
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_balmer'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_balmer: {len(par_dict_rev['comments_on_balmer'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_tied: {len(par_dict_rev['comments_on_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_for_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_for_tied: {len(par_dict_rev['factor_for_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_fixed_for_tying'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_fixed_for_tying: {len(par_dict_rev['factor_fixed_for_tying'])}"
	flux = kwargs.get('flux', np.zeros_like(wave))  # flux_sky
	flux_sky = kwargs.get('flux_sky', np.zeros_like(wave))  # flux_sky
	cont_init = kwargs.get('cont_init', np.ones_like(wave))  # cont_init
	continuum_fitted = cont_init
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	FWHM_gal1 = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	if (FWHM_gal1<4.0):
		FWHM_gal1 = 4.0
	resolution = float(c_kms / FWHM_gal1)
	popt_init, masked_array, amp_init_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, coeff_init_val, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp = emff.get_opt_param_array_rev2(wave, flux, cont_init, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], par_dict_rev['redshift_val'], par_dict_rev['comments_on_balmer'], par_dict_rev['comments_on_tied'], par_dict_rev['factor_for_tied'], par_dict_rev['e_b_minus_v_init'], par_dict_rev['window_for_fit'], par_dict_rev['stopping_number_for_continuum'], par_dict_rev['center_init_array'], par_dict_rev['sigma_init_array'])
	amp_length = len(amp_init_array)
	amp_array, center_array, sigma_array, reddening_val_fit, coeff_fit = emff.get_params(popt_fit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	amp_array_rev = emff.retrieve_all_amplitude_list_rev2(list(amp_array), par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, par_dict_rev['comments_on_tied'], par_dict_rev['comments_on_balmer'])
	amp_err_array, center_err_array, sigma_err_array, reddening_err_val_fit, coeff_err_fit = emff.get_params(perr_fit, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], amp_length)
	amp_err_array_rev = emff.retrieve_all_amplitude_list_rev2(list(amp_err_array), par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, par_dict_rev['comments_on_tied'], par_dict_rev['comments_on_balmer'])

	new_data = []
	new_data.extend([amp_array_rev])
	new_data.extend([center_array])
	new_data.extend([sigma_array])
	new_data.extend([reddening_val_fit])
	new_data_err = []
	new_data_err.extend([amp_err_array_rev])
	new_data_err.extend([center_err_array])
	new_data_err.extend([sigma_err_array])
	new_data_err.extend([reddening_err_val_fit])
	return (new_data, new_data_err)

###############################GET_EXPANDED_DATA_INFORMATION_EMISSION###############################

###############################GET_EXPANDED_HEADER_INFORMATION_EMISSION###############################

def get_expanded_header_information_emission(wave, par_dict, **kwargs):
	par_dict_rev = par_dict
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_balmer'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_balmer: {len(par_dict_rev['comments_on_balmer'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['comments_on_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and comments_on_tied: {len(par_dict_rev['comments_on_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_for_tied'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_for_tied: {len(par_dict_rev['factor_for_tied'])}"
	assert (len(par_dict_rev['center_list'])==len(par_dict_rev['factor_fixed_for_tying'])), f"Length mismatch between center_list: {len(par_dict_rev['center_list'])} and factor_fixed_for_tying: {len(par_dict_rev['factor_fixed_for_tying'])}"
	flux = kwargs.get('flux', np.zeros_like(wave))  # flux_sky
	flux_sky = kwargs.get('flux_sky', np.zeros_like(wave))  # flux_sky
	cont_init = kwargs.get('cont_init', np.ones_like(wave))  # cont_init
	continuum_fitted = cont_init
	quiet_val = kwargs.get('quiet_val', True)  # quiet_val
	FWHM_gal1 = kwargs.get('FWHM_gal', 5.0)  # FWHM_gal
	if (FWHM_gal1<4.0):
		FWHM_gal1 = 4.0
	resolution = float(c_kms / FWHM_gal1)
	popt_init, masked_array, amp_init_array, position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, coeff_init_val, index_for_fixing_bounds_factor_init_narrow_comp, index_for_fixing_bounds_factor_init_wide_comp = emff.get_opt_param_array_rev2(wave, flux, cont_init, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], par_dict_rev['redshift_val'], par_dict_rev['comments_on_balmer'], par_dict_rev['comments_on_tied'], par_dict_rev['factor_for_tied'], par_dict_rev['e_b_minus_v_init'], par_dict_rev['window_for_fit'], par_dict_rev['stopping_number_for_continuum'], par_dict_rev['center_init_array'], par_dict_rev['sigma_init_array'])
	amp_length = len(amp_init_array)
	if ('custom'in par_dict_rev['emission_fit_type']):
		amp_array, center_array, sigma_array, reddening_val_fit, coeff_fit = emff.get_params(popt_init, par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], par_dict_rev['center_list'], len(amp_init_array))
		amp_array_rev = emff.retrieve_all_amplitude_list_rev2(list(amp_array), par_dict_rev['number_of_narrow_components'], par_dict_rev['number_of_wide_components'], position_init_narrow_comp, position_init_wide_comp, position_final_narrow_comp, position_final_wide_comp, par_dict_rev['comments_on_tied'], par_dict_rev['comments_on_balmer'])
		header_information_emission = []
		extra_label = []
		header_information_emission.extend(['Amp-']*len(amp_array_rev))
		index_cust = 0
		for i in range(len(amp_array_rev)):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['V-']*len(center_array))
		index_cust = 0
		for i in range(len(center_array)):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['sigma-']*len(sigma_array))
		index_cust = 0
		for i in range(len(sigma_array)):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['E(B-V)(gas)'])
		extra_label.extend(['-']*1)
		header_information_emission2 = [str(m)+str(n) for m,n in zip(header_information_emission,extra_label)]
		header_information_err_emission2 = [s + "_err" for s in header_information_emission2]
		header_information_emission3 = np.array(header_information_emission2)
		header_information_err_emission3 = np.array(header_information_err_emission2)
	elif ('ppxf'in par_dict_rev['emission_fit_type']):
		if (par_dict['ppxf_gas_comp']<1):
			par_dict['ppxf_gas_comp']=1
		pathname_ppxf = ppxf_dir + '/miles_models/*.fits'  # Initialize a default stellar template
		miles = lib.miles(pathname_ppxf, 3.0, 55.0)
		lam_range_gal = np.array([np.min(wave), np.max(wave)])/(1. + par_dict['redshift_val'])
		bf.blockPrint()
		gas_templates, gas_names, line_wave = util.emission_lines(miles.ln_lam_temp, lam_range_gal, 3.0, tie_balmer=par_dict['ppxf_emission_tie_balmer'], limit_doublets=par_dict['ppxf_emission_limit_doublets'])
		bf.enablePrint()
		header_information_emission = []
		extra_label = []
		header_information_emission.extend(['V-']*int(par_dict['ppxf_gas_comp']))
		index_cust = 0
		for i in range(int(par_dict['ppxf_gas_comp'])):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['sigma-']*int(par_dict['ppxf_gas_comp']))
		index_cust = 0
		for i in range(int(par_dict['ppxf_gas_comp'])):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['V_other-']*int(par_dict['ppxf_gas_comp']))
		index_cust = 0
		for i in range(int(par_dict['ppxf_gas_comp'])):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['sigma_other-']*int(par_dict['ppxf_gas_comp']))
		index_cust = 0
		for i in range(int(par_dict['ppxf_gas_comp'])):
			index_cust+=1
			extra_label.extend([index_cust]*1)
		header_information_emission.extend(['Amp-']*len(gas_names))
		extra_label.extend(gas_names)
		#header_information_emission.extend(['E(B-V)(gas)'])
		#extra_label.extend(['-']*1)
		header_information_emission.extend(['chi-sq'])
		extra_label.extend(['-']*1)
		header_information_emission2 = [str(m)+str(n) for m,n in zip(header_information_emission,extra_label)]
		header_information_err_emission2 = [s + "_err" for s in header_information_emission2]
		header_information_emission3 = np.array(header_information_emission2)
		header_information_err_emission3 = np.array(header_information_err_emission2)
		header_information_err_emission3[-1] = "E(B-V)gas"
	return (header_information_emission3, header_information_err_emission3, amp_length)

###############################GET_EXPANDED_HEADER_INFORMATION_EMISSION###############################

###############################GET_EXPANDED_INFORMATION_ABSORPTION###############################

def get_expanded_information_abs(par_dict):
	header_information_absorption = []
	extra_label = []
	header_information_absorption.extend(['V-', 'sigma', 'h3-', 'h4-']*par_dict['ppxf_stars_comp'])
	index_stars = 0
	for i in range(par_dict['ppxf_stars_comp']):
		index_stars+=1
		extra_label.extend([index_stars]*4)
	header_information_absorption.extend(['E(B-V)', 'E(B-V)(gas)'])
	extra_label.extend(['-']*2)
	header_information_absorption.extend(['V(gas-balmer)-', 'sigma(gas-balmer)-', 'V(gas-other)-', 'sigma(gas-other)-']*par_dict['ppxf_gas_comp'])
	index_gas = 0
	for i in range(par_dict['ppxf_gas_comp']):
		index_gas+=1
		extra_label.extend([index_gas]*4)
	header_information_absorption.extend(['age']*len(par_dict['st_age_unique_end']))
	header_information_absorption.extend(['mass']*len(par_dict['st_age_unique_end']))
	header_information_absorption.extend(['lum']*len(par_dict['st_age_unique_end']))
	extra_label.extend(par_dict['st_age_unique_end']*3)
	header_information_absorption.extend(['StAge', 'StMetallicity'])
	extra_label.extend(['-']*2)
	header_information_absorption2 = [str(m)+str(n) for m,n in zip(header_information_absorption,extra_label)]
	header_information_err_absorption = [s + "_err" for s in header_information_absorption2]
	header_information_absorption3 = np.array(header_information_absorption2)
	header_information_err_absorption3 = np.array(header_information_err_absorption)
	return (header_information_absorption3, header_information_err_absorption)

###############################GET_EXPANDED_INFORMATION_ABSORPTION###############################

########################################STAGE_FIVE_ANALYSIS########################################


















