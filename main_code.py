import matplotlib.pyplot as plt
#import wx
#import matplotlib.backends.backend_wxagg
import time
import warnings
import os
import sys
import pyfiglet
from astropy.table import Table
from astropy.io import fits
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astropy.wcs import WCS

#import pyfits
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
sys.path.append(str(sys.argv[2]) + '/custom_functions')
#from custom_functions_by_axe import *
import custom_functions_by_axe as cfba
#from basic_functions import *
import basic_functions as bf
import get_data_from_file as gdff
import handling_fits_files as hff
import spatial_binning_functions as sbf
import plotting_functions as pf
#
ascii_banner = pyfiglet.figlet_format("MUSELI", font='isometric1', justify='center')
bf.print_cust(ascii_banner)
bf.print_cust('Authors: Axe, Saulder, Martin, and Bag')
# 크리스, Martin & Axe corp.')
#
bf.print_cust(str("Successfully imported all modules...." + sys.argv[2]))
#
bf.print_cust(str("Running script: " + sys.argv[0]))
bf.print_cust(str("Current Working Directory: " + sys.argv[1]))
bf.print_cust(str("Base Working Directory: " + sys.argv[2]))
#
#######################OBTAINING_PARAMETER_INFORMATION#######################
#Obtains parameter information from a given parameter file
#If that is not available, it obtains default parameters from code,
#except for the main filename.
#
parameter_file_string_current = str(sys.argv[1]) + str('/parameter_file.dat')
par_dict_init = gdff.initial_guess_from_parfile(parameter_file_string_current)
dir_name_6 = str(sys.argv[2]) + str('/custom_functions/basic_files')
par_dict = gdff.revise_dictionary(par_dict_init, dir_name_6)
par_dict['quiet'] = True
#
#Telling the code to not print anything on screen (It will still be there in the logfile regardless)
quiet_val = bf.str2bool(par_dict['quiet'])
if not quiet_val: warnings.simplefilter(action='ignore', category=FutureWarning)
if not quiet_val: np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
if not quiet_val: np.seterr(divide = 'ignore')
cfba.quiet_val = quiet_val

#######################MAKING_DIRECTORIES#######################

dir_name_1 = str(sys.argv[1]) + str("/results")
dir_name_2 = str(sys.argv[1]) + str("/images")
dir_name_3 = str(sys.argv[1]) + str("/post_process_images")
dir_name_4 = str(sys.argv[1]) + str("/data")
bf.check_directory(dir_name_1)
bf.check_directory(dir_name_2)
bf.check_directory(dir_name_3)
bf.check_directory(dir_name_4)


#######################MAKING_DIRECTORIES#######################

'''
with fits.open(str(str(par_dict['muse_data_filename']))) as hdul:
	header_original = hdul[1].header
	header_original_err = hdul[2].header
	data_original_file = hdul[1].data
	err_original_file = hdul[2].data

print (repr(header_original))
primary_hdu = fits.PrimaryHDU(header=header_original)
c1 = fits.Column(name='a', array=np.array(['a', 'b', 'c', 'd']), format='A20V')
c2 = fits.Column(name='b', array=np.array(['e', 'f']), format='A20V')
c3 = fits.Column(name='c', array=np.array(['g', 'h', 'i']), format='A20V')
c4 = fits.Column(name='d', array=np.array(['j', 'kkk', 'lk', 'mj', 'we', '12e', 'asd3', '4gsd4']), format='A20V')
table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4])
hdul = fits.HDUList([primary_hdu, table_hdu])
hdul.writeto('table.fits', overwrite=True)

with fits.open('table.fits', 'append') as hdul:
	hdr_test = hdul[0].header
	d1 = hdul[1].data['a']
	d2 = hdul[1].data['b']
	d3 = hdul[1].data['c']
	d4 = hdul[1].data['d']

print (repr(hdr_test))
print (d1)
print (d2)
print (d3)
print (d4)
quit()
with fits.open(str(str(par_dict['muse_data_filename']))) as hdul:
header_original = hdul[1].header
muse_wcs = WCS(header_original).celestial
data2d = np.nanmean(data_original_file, axis=(0))
fig = plt.figure()
ax1 = plt.subplot(111, projection=muse_wcs)
ax1.imshow(data2d, origin='lower')
plt.show()
quit()
'''

#######################EXECUTING_STAGE_ONE#######################

new_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_revised.fits')
if (os.path.exists(new_filename)):
	prompt_1=input(f"File {new_filename} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_1.lower()):
		bf.print_cust(f"deleting...{new_filename}")
		string_for_deleting = str('rm ') + new_filename
		bf.print_cust(string_for_deleting)
		os.system(string_for_deleting)
		cfba.stage_one_analysis(par_dict, new_filename, quiet_val=quiet_val)
else:
	cfba.stage_one_analysis(par_dict, new_filename, quiet_val=quiet_val)

#######################EXECUTING_STAGE_ONE#######################


'''
#STAGE ONE LOOKS GOOD
wave_min = 7082.0
wave_max = 7107.0
fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

header_original, header_original_err, data_original_file, err_original_file = hff.open_ifu_fits(str(par_dict['muse_data_filename']), str(par_dict['observed_instrument']))
ra_array, dec_array, wave_array = hff.obtain_physical_axis(header_original)
idx1 = bf.find_nearest_idx(wave_array, wave_min)
idx2 = bf.find_nearest_idx(wave_array, wave_max)
imdata1 = np.nanmedian(data_original_file[idx1:idx2, :, :,], axis=(0))
im1 = ax1.imshow(imdata1, origin='lower', vmin=0, vmax=500)
bf.add_colorbar(im1)

header_rev, header_rev_err, header_rev_cont, header_rev_fwhm, header_rev_velscale, header_rev_sky, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, velscale_rev_file, sky_rev_file = hff.open_ifu_fits_custom_file(new_filename)
ra_rev, dec_rev, wave_rev = hff.obtain_physical_axis(header_rev)
idx1 = bf.find_nearest_idx(wave_rev, wave_min)
idx2 = bf.find_nearest_idx(wave_rev, wave_max)
imdata2 = np.nanmedian(data_rev_file[idx1:idx2, :, :,], axis=(0))
im2 = ax2.imshow(imdata2, origin='lower', vmin=0, vmax=500)
bf.add_colorbar(im2)
#plt.show()
plt.savefig('first_stage.pdf')
plt.close()
#quit()
'''

#######################EXECUTING_STAGE_TWO#######################

file_name_eq_width = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_eq_width.fits')
if (os.path.exists(file_name_eq_width)):
	prompt_2=input(f"File {file_name_eq_width} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_2.lower()):
		bf.print_cust(f"deleting...{file_name_eq_width}")
		string_for_deleting = str('rm ') + file_name_eq_width
		bf.print_cust(string_for_deleting)
		os.system(string_for_deleting)
		cfba.stage_two_analysis(par_dict, new_filename, file_name_eq_width)
else:
		cfba.stage_two_analysis(par_dict, new_filename, file_name_eq_width)

#######################EXECUTING_STAGE_TWO#######################

quit()
'''
#To obtain EW information from EW data cube
header_ew_original, header_ew_original_err, data_ew_original_file, err_ew_original_file = hff.open_ifu_fits(file_name_eq_width)
lick_index_species_list = par_dict['lick_index_species_list']
lick_index_wave11, lick_index_wave12, lick_index_wave21, lick_index_wave22, lick_index_wave31, lick_index_wave32, lick_index_sign,lick_index_species, lick_index_reference = np.genfromtxt(str(par_dict['lick_index_file']), unpack=True, names=True, encoding=None, dtype=None)
mean_wave_list_full = (lick_index_wave21 + lick_index_wave22) / 2.
lick_index_species_rev = np.array(lick_index_species[lick_index_species_list[:]])
with fits.open(str(str(par_dict['muse_data_filename']))) as hdul:
    header_original = hdul[1].header
muse_wcs = WCS(header_original).celestial
print (lick_index_species_rev)
#names_of_species = str(lick_index_species_rev[:]) + str(mean_wave_list_full.astype(np.int32)[:])
#print (names_of_species)
#quit()
#print (data_ew_original_file.shape[0])
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(15, 15),subplot_kw=dict(projection=muse_wcs))
#ax1 = fig.add_subplot(121)
counter_test = 0
data_ew_original_file[data_ew_original_file<-1000.] = np.nan
data_ew_original_file[data_ew_original_file>1000.0] = np.nan
data_ew_original_file = np.abs(data_ew_original_file)
count = 0
counter_test = [0, 2, 3, 6]
for i in range(2):
    for j in range(2):
        if ('5007' in  lick_index_species_rev[counter_test[count]]) or ('alpha' in  lick_index_species_rev[counter_test[count]]):
            #im1 = axs[i, j].imshow(data_ew_original_file[counter_test[count], :, :], vmin=5, vmax=100, cmap='viridis')
            im1 = axs[i, j].imshow(np.log10(data_ew_original_file[counter_test[count], :, :]), vmin=-1, vmax=2, cmap='viridis')
        else:
            #im1 = axs[i, j].imshow(data_ew_original_file[counter_test[count], :, :], vmin=0.1, vmax=100, cmap='viridis')
            im1 = axs[i, j].imshow(np.log10(data_ew_original_file[counter_test[count], :, :]), vmin=-1, vmax=2, cmap='viridis')

        divider = make_axes_locatable(axs[i, j])
        cax = divider.append_axes('right', size='5%', pad=0.55)
        fig.colorbar(im1, cax=cax, orientation='vertical')
        axs[i, j].set_title(lick_index_species_rev[counter_test[count]])
        #axs[i, j].invert_xaxis()
        #counter_test+=1
        count+=1

#ax2 = fig.add_subplot(122)
#im2 = ax2.imshow(err_ew_original_file[-1, :, :], origin='lower', vmin=-10, vmax=1)
#divider = make_axes_locatable(ax2)
#cax = divider.append_axes('right', size='5%', pad=0.05)
#fig.colorbar(im2, cax=cax, orientation='vertical')
#fig.suptitle(r'H$\alpha$ Equivalent Width Map (Data/Error)')
#plt.show()
plt.savefig('test_fig.pdf', dpi=100)
quit()
'''

#To obtain information from revised data
bf.print_cust('Loading information from revised data...')
header_rev, header_rev_err, header_rev_cont, header_rev_fwhm, header_rev_velscale, header_rev_sky, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, velscale_rev_file, sky_rev_file = hff.open_ifu_fits_custom_file(new_filename)
ra_rev, dec_rev, wave_rev = hff.obtain_physical_axis(header_rev)
#print (header_rev)
#print (ra_rev, dec_rev)
#quit()
header_ew_original, header_ew_original_err, data_ew_original_file, err_ew_original_file = hff.open_ifu_fits(file_name_eq_width)
bf.print_cust('information from revised data loaded...')
#bf.print_cust(f'{data_ew_original_file.shape}')
#bf.print_cust(header_ew_original)
#print (par_dict['lick_idx_for_halpha'])
halpha_eq_width_map = data_ew_original_file[int(par_dict['lick_idx_for_halpha']), :, :]
#halpha_eq_width_map = data_ew_original_file[int(par_dict['lick_idx_for_halpha']), :, :]
halpha_eq_width_err_map = err_ew_original_file[int(par_dict['lick_idx_for_halpha']), :, :]
#halpha_eq_width_err_map = err_ew_original_file[-1, :, :]
halpha_eq_width_map[np.abs(halpha_eq_width_map)>200.] = 1e-6
halpha_eq_width_err_map[np.abs(halpha_eq_width_map)>200.] = 1e-5

'''
#par_dict = gdff.revise_dictionary(par_dict)
#pos_x = 162
#pos_y = 168
#pos_x = 150
#pos_y = 150
pos_x = 105
pos_y = 250
flux1 = data_rev_file[:, pos_x, pos_y]
err1 = np.sqrt(err_rev_file[:, pos_x, pos_y])
halpha_ew_val1 = halpha_eq_width_map[pos_x, pos_y]
FWHM_gal1 = np.nanmedian(fwhm_rev_file[:, pos_x, pos_y])
flux_sky1 = sky_rev_file
init_cont1 = cont_rev_file[:, pos_x, pos_y]
par_dict['ppxf_stars_comp'] = 1
par_dict['ppxf_gas_comp'] = 1
par_dict['ppxf_emission_tie_balmer'] = True
par_dict['ppxf_emission_limit_doublets'] = False
par_dict = get_emission_keys(par_dict, wave_rev)
#wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = main_fitting_func(wave_rev, flux1, err1, par_dict, halpha_ew_val1, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1, quiet_run=False)
file1 = np.array([wave_rev, flux1, err1])
lam_gal_rev, fit_galaxy, fit_noise, fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, fit_normalising_factor = ppxf_function(file1, redshift=par_dict['redshift_val'], number_of_stellar_components=par_dict['ppxf_stars_comp'], number_of_gas_components=par_dict['ppxf_gas_comp'], sky_flux_init=flux_sky1, first_order_kinematic_limits=1000., tie_balmer=par_dict['ppxf_emission_tie_balmer'], limit_doublets=par_dict['ppxf_emission_limit_doublets'], quiet=False)

bestfit_solution_array, stellar_spectrum_array, sol_pp_gas_flux_reshaped, sol_pp_gas_flux_err_reshaped, stellar_vel, stellar_vel_err, stellar_sigma, stellar_sigma_err, stellar_h3, stellar_h3_err, stellar_h4, stellar_h4_err, gas_balmer_vel, gas_balmer_vel_err, gas_balmer_sigma, gas_balmer_sigma_err, gas_others_vel, gas_others_vel_err, gas_others_sigma, gas_others_sigma_err, status_of_optimization, reddening_fit_e_b_minus_v, gas_reddening_fit_e_b_minus_v, goodpixels_idx, reduced_chi_sq, st_age_unique, st_mass_unique, st_lum_unique, weighted_logAge, weighted_metallicity, mlpop_sol, str_fit = get_detailed_ppxf_solution(fit_sol_pp, fit_miles, fit_component_rev, fit_gas_names, fit_gas_names_rev, par_dict['ppxf_stars_comp'], par_dict['ppxf_gas_comp'], normalising_factor = fit_normalising_factor, st_age_unique_start=par_dict['st_age_list_start'], st_age_unique_end=par_dict['st_age_unique_end'], quiet=False)

plt.errorbar(lam_gal_rev, fit_galaxy, yerr=fit_noise, ds='steps-mid', color='tab:blue', zorder=1, label='data')
plt.plot(lam_gal_rev, bestfit_solution_array, 'r--', zorder=2, label='fit')
plt.plot(lam_gal_rev, stellar_spectrum_array, 'g--', zorder=2, label='star')
plt.legend()
plt.show()
quit()
'''

#######################EXECUTING_STAGE_THREE#######################
#par_dict = gdff.revise_dictionary(par_dict)
#par_dict['region_of_binning_interest'] = [50, 50, 250, 250]
file_name_rev_linear = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_linear_snr.dat')
#x_array, y_array, ra_array, dec_array, signal1, noise1, signal2, noise2 = np.loadtxt(file_name_rev_linear).T
if (os.path.exists(file_name_rev_linear)):
	prompt_3=input(f"File {file_name_rev_linear} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_3.lower()):
		bf.print_cust(f"deleting...{file_name_rev_linear}")
		string_for_deleting = str('rm ') + file_name_rev_linear
		bf.print_cust(string_for_deleting)
		os.system(string_for_deleting)
		sbf.create_linear_file_for_binning(file_name_rev_linear, par_dict, header_rev, cont_rev_file, data_rev_file, err_rev_file, data_ew_original_file, err_ew_original_file)
else:
	sbf.create_linear_file_for_binning(file_name_rev_linear, par_dict, header_rev, cont_rev_file, data_rev_file, err_rev_file, data_ew_original_file, err_ew_original_file)

file_name_rev_binned = cfba.stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict)

#######################EXECUTING_STAGE_THREE#######################

'''
y1, x1, ra1, dec1, signal11, noise11, signal12, noise12 = np.loadtxt(file_name_rev_linear).T
fig, axs = plt.subplots(2, figsize=(6,10), sharex=True, sharey=True)
im1 = axs[0].scatter(x1, y1, c = signal11/noise11, cmap='viridis', vmin=0, vmax=20)
axs[0].set_title("SNR")
add_colorbar_lin(im1)
x_new, y_new, signal_new = np.loadtxt(file_name_rev_binned).T
im2 = axs[1].scatter(x_new, y_new, c = signal_new, cmap='rainbow')
axs[1].set_title("Binning")
add_colorbar_lin(im2)
plt.show()
quit()
'''

#######################EXECUTING_STAGE_FOUR#######################
#file_name_rev_binned = cfba.stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict)
#par_dict = gdff.revise_dictionary(par_dict)
x_axis_binned_rev, y_axis_binned_rev, bin_num = np.loadtxt(file_name_rev_binned).T
bin_num_unique = np.unique(bin_num)
main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = cfba.get_output_filenames_and_arrays(dir_name_4, bin_num_unique, wave_rev, par_dict)

cfba.stage_four_execution(par_dict, file_name_rev_linear, file_name_rev_binned, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file)
#######################EXECUTING_STAGE_FOUR#######################

'''
x_axis_binned_rev = x_axis_binned_rev-1
y_axis_binned_rev = y_axis_binned_rev-1
x_axis_unique = np.unique(y_axis_binned_rev)
y_axis_unique = np.unique(x_axis_binned_rev)
ew_sum_array = np.zeros([len(x_axis_unique), len(y_axis_unique)])
#print (halpha_eq_width_map.shape)
for i in range(len(bin_num_unique)):
	halpha_sum = []
	#halpha_sum = 0.0
	mask = np.in1d(bin_num, [bin_num_unique[i]])
	halpha_sum = np.nanmean(halpha_eq_width_map[y_axis_binned_rev[mask].astype(int), x_axis_binned_rev[mask].astype(int)])
	x_indexes = x_axis_binned_rev[mask.astype(bool)].astype(int)
	y_indexes = y_axis_binned_rev[mask.astype(bool)].astype(int)
	#for k in range(len(y_indexes)):
		#halpha_sum.extend([(halpha_eq_width_map[y_indexes[k]-1, x_indexes[k]-1])])
		#halpha_sum+=np.nan_to_num(halpha_eq_width_map[y_indexes[k]-1, x_indexes[k]-1])
	for j in range(len(y_indexes)):
		ew_sum_array[y_indexes[j], x_indexes[j]] = halpha_sum

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
#ew_sum_array = np.transpose(ew_sum_array)
ax11 = ax1.imshow(halpha_eq_width_map, origin='lower')
fig.colorbar(ax11, ax=ax1)

#v_expanded_array = np.transpose(v_expanded_array)
ax22 = ax2.imshow(ew_sum_array, origin='lower')
fig.colorbar(ax22, ax=ax2)
plt.show()
#quit()

#x_axis_unique = np.unique(y_axis_binned_rev)
#y_axis_unique = np.unique(x_axis_binned_rev)
test = np.load(main_absorption_result_cube_filename)
v_array = test[:, 0, 0]
#plt.hist(v_array)
#plt.show()
v_expanded_array = np.zeros([len(x_axis_unique), len(y_axis_unique)])
for i in range(len(bin_num_unique)):
	mask = np.in1d(bin_num, [bin_num_unique[i]])
	x_indexes = x_axis_binned_rev[mask.astype(bool)].astype(int)
	y_indexes = y_axis_binned_rev[mask.astype(bool)].astype(int)
	for j in range(len(y_indexes)):
		v_expanded_array[y_indexes[j]-1, x_indexes[j]-1] = v_array[i]


fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
#ew_sum_array = np.transpose(ew_sum_array)
ax11 = ax1.imshow(ew_sum_array, origin='lower', vmin=-10, vmax=10)
fig.colorbar(ax11, ax=ax1)

#v_expanded_array = np.transpose(v_expanded_array)
ax22 = ax2.imshow(v_expanded_array, origin='lower')
fig.colorbar(ax22, ax=ax2)
plt.show()
quit()
'''

#######################EXECUTING_STAGE_FIVE#######################

expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.hdf5')
expanded_hdr_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_hdr.fits')

expanded_hdr_filename1, expanded_filename1 = cfba.stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_linear, file_name_rev_binned, expanded_hdr_filename, header_rev, expanded_filename)

#######################EXECUTING_STAGE_FIVE#######################

snr_vmin_val1 = 0.1
snr_vmax_val1 = 20
physical_axes1 = False
x_array_binned, y_array_binned, snr_revised = pf.get_snr_map_revised(file_name_rev_linear, file_name_rev_binned, par_dict, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)

fig, ax = plt.subplots()
ax.cla()
im = ax.scatter(x_array_binned, y_array_binned, c = snr_revised, cmap='viridis', vmin=snr_vmin_val1, vmax=snr_vmax_val1)
im_cl = bf.add_colorbar_lin(im)

'''
def update(val):
	bf.print_cust("Current Status: ")
	binning_quant_updated = int(sbinning_quant_slider.val)
	bf.print_cust(binning_quant_updated)
	sspectral_smoothing_updated = int(sspectral_smoothing_slider.val)
	bf.print_cust(binning_quant_updated)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	bf.print_cust(binning_quant_updated)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	bf.print_cust(binning_quant_updated)
	radio_work_type_updated = str(radio_work_type.value_selected)
	bf.print_cust(binning_quant_updated)
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	bf.print_cust(binning_quant_updated)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	bf.print_cust(binning_quant_updated)
'''

binning_quant_slider_axes = plt.axes([0.1, 0.01, 0.3, 0.02])
sbinning_quant_slider = Slider(binning_quant_slider_axes, 'binning_quant', 2, 100, valinit=int(par_dict['binning_quant']), valfmt="%i")
#sbinning_quant_slider.on_changed(update)

#spectral_smoothing_slider_axes = plt.axes([0.4, 0.01, 0.1, 0.02])
#sspectral_smoothing_slider = Slider(spectral_smoothing_slider_axes, 'spectral_smoothing', 1, 100, valinit=int(par_dict['spectral_smoothing']), valfmt="%i")
#sspectral_smoothing_slider.on_changed(update)

binning_type_buttons_axes = plt.axes([0.01, 0.2, 0.1, 0.2])
dict_binning_type_buttons = {'None', 'square', 'voronoi', 'advanced'}
active_binning_type = list(dict_binning_type_buttons).index(str(par_dict['binning_type']))
#print (active_binning_type)
#print (list(dict_binning_type_buttons)[active_binning_type])
radio_binning_type_buttons = RadioButtons(binning_type_buttons_axes, ('None', 'square', 'voronoi', 'advanced'), active=2)
def func_binning_type_buttons(label_binning_type_buttons):
	bf.print_cust(label_binning_type_buttons)
radio_binning_type_buttons.on_clicked(func_binning_type_buttons)

binning_significance_type_buttons_axes = plt.axes([0.01, 0.1, 0.1, 0.1])
dict_binning_significance_type = {'snr', 'ew'}
active_binning_significance_type = list(dict_binning_significance_type).index(str(par_dict['voronoi_snr_type']))
#print (active_binning_significance_type)
#print (list(dict_binning_significance_type)[active_binning_significance_type])
radio_binning_significance_type = RadioButtons(binning_significance_type_buttons_axes, ('snr', 'ew'), active=0)
def func_binning_significance_type(label_binning_significance_type):
	bf.print_cust(label_binning_significance_type)
radio_binning_significance_type.on_clicked(func_binning_significance_type)

#work_type_buttons_axes = plt.axes([0.01, 0.65, 0.1, 0.1])
#radio_work_type = RadioButtons(work_type_buttons_axes, ('snr_map', 'fit'), active=0)
#dict_work_type = {'SNR Map', 'Fit'}
#def func_work_type(label_work_type):
#	bf.print_cust(label_work_type)
#radio_work_type.on_clicked(func_work_type)

fit_type_buttons_axes = plt.axes([0.01, 0.5, 0.1, 0.15])
dict_fit_type = {'auto', 'emission', 'absorption'}
active_dict_fit_type = list(dict_fit_type).index(str(par_dict['execution_fit_type']))
#print (active_dict_fit_type)
#print (list(dict_fit_type)[active_dict_fit_type])
radio_fit_type = RadioButtons(fit_type_buttons_axes, ('auto', 'emission', 'absorption'), active=0)
def func_fit_type(label_fit_type):
	bf.print_cust(label_fit_type)
radio_fit_type.on_clicked(func_fit_type)

emission_fit_type_buttons_axes = plt.axes([0.01, 0.4, 0.1, 0.1])
dict_emission_fit_type = {'custom', 'ppxf'}
active_dict_emission_fit_type = list(dict_emission_fit_type).index(str(par_dict['emission_fit_type']))
#print (active_dict_emission_fit_type)
#print (list(dict_emission_fit_type)[active_dict_emission_fit_type])
radio_emission_fit_type = RadioButtons(emission_fit_type_buttons_axes, ('custom', 'ppxf'), active=1)
def func_emission_fit_type(label_emission_fit_type):
	bf.print_cust(label_emission_fit_type)
radio_emission_fit_type.on_clicked(func_emission_fit_type)

def update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated):
	bf.print_cust('Updating Figure...')
	global im_cl
	im = ax.scatter(x_array_binned_updated, y_array_binned_updated, c = snr_revised_updated, cmap='viridis', vmin=snr_vmin_val1, vmax=snr_vmax_val1)
	im_cl = bf.add_colorbar_lin(im)
	bf.print_cust('Figure updated...')
	sys.stdout.flush()
	time.sleep(0.1)
	plt.draw()

bin_data_buttons_axes = plt.axes([0.01, 0.88, 0.1, 0.02])
bin_data_button = Button(bin_data_buttons_axes, 'Bin Data')
def func_bin_data(event):
	bf.print_cust('Binning Data...')
	global im_cl
	binning_quant_updated = int(sbinning_quant_slider.val)
	sspectral_smoothing_updated = int(1)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	radio_work_type_updated = str('snr_map')
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	bf.print_cust(f"Current Status: Spectral Smoothing Factor: {sspectral_smoothing_updated},  Binning Quant: {binning_quant_updated},  Radio Binning Type: {radio_binning_type_updated},  Radio Binning Significance Type: {radio_binning_significance_type_updated},  Radio Work Type: {radio_work_type_updated},  Radio Fit Type: {radio_fit_type_updated},  Radio Emission Fit Type: {radio_emission_fit_type_updated}")
	par_dict_updated = gdff.modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated)
	file_name_rev_binned_updated = cfba.stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict_updated)
	x_axis_binned_rev_updated, y_axis_binned_rev_updated, bin_num_updated = np.loadtxt(file_name_rev_binned_updated).T
	bin_num_unique_updated = np.unique(bin_num_updated)
	main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = cfba.get_output_filenames_and_arrays(dir_name_4, bin_num_unique_updated, wave_rev, par_dict_updated)
	cfba.stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, plot_val=False)
	sys.stdout.flush()
	time.sleep(0.1)
	x_array_binned_updated, y_array_binned_updated, snr_revised_updated = pf.get_snr_map_revised(file_name_rev_linear, file_name_rev_binned_updated, par_dict_updated, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)
	ax.cla()
	im_cl.remove()
	update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated)
bin_data_button.on_clicked(func_bin_data)


fit_data_buttons_axes = plt.axes([0.01, 0.82, 0.1, 0.02])
fit_data_button = Button(fit_data_buttons_axes, 'Fit Data')
def func_fit_data(event):
	bf.print_cust('Binning Data...')
	global im_cl
	binning_quant_updated = int(sbinning_quant_slider.val)
	sspectral_smoothing_updated = int(1)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	radio_work_type_updated = str('snr_map')
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	bf.print_cust(f"Current Status: Spectral Smoothing Factor: {sspectral_smoothing_updated},  Binning Quant: {binning_quant_updated},  Radio Binning Type: {radio_binning_type_updated},  Radio Binning Significance Type: {radio_binning_significance_type_updated},  Radio Work Type: {radio_work_type_updated},  Radio Fit Type: {radio_fit_type_updated},  Radio Emission Fit Type: {radio_emission_fit_type_updated}")
	par_dict_updated = gdff.modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated)
	file_name_rev_binned_updated = cfba.stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict_updated)
	x_axis_binned_rev_updated, y_axis_binned_rev_updated, bin_num_updated = np.loadtxt(file_name_rev_binned_updated).T
	bin_num_unique_updated = np.unique(bin_num_updated)
	main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = cfba.get_output_filenames_and_arrays(dir_name_4, bin_num_unique_updated, wave_rev, par_dict_updated)
	cfba.stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, plot_val=False)
	sys.stdout.flush()
	time.sleep(0.1)
	par_dict_updated['execution_type'] = 'fit'
	cfba.stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, plot_val=False)
	sys.stdout.flush()
	time.sleep(0.1)
	x_array_binned_updated, y_array_binned_updated, snr_revised_updated = pf.get_snr_map_revised(file_name_rev_linear, file_name_rev_binned_updated, par_dict_updated, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)
	ax.cla()
	im_cl.remove()
	update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated)
fit_data_button.on_clicked(func_fit_data)





show_results_buttons_axes = plt.axes([0.01, 0.76, 0.1, 0.02])
show_results_button = Button(show_results_buttons_axes, 'Show Results')
def func_show_results(event):
	bf.print_cust('Show Results...')
	binning_quant_updated = int(sbinning_quant_slider.val)
	sspectral_smoothing_updated = int(1)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	radio_work_type_updated = str('snr_map')
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	bf.print_cust(f"Current Status: Spectral Smoothing Factor: {sspectral_smoothing_updated},  Binning Quant: {binning_quant_updated},  Radio Binning Type: {radio_binning_type_updated},  Radio Binning Significance Type: {radio_binning_significance_type_updated},  Radio Work Type: {radio_work_type_updated},  Radio Fit Type: {radio_fit_type_updated},  Radio Emission Fit Type: {radio_emission_fit_type_updated}")
	par_dict_updated = gdff.modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated)
	file_name_rev_binned_updated = cfba.stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict_updated)
	x_axis_binned_rev_updated, y_axis_binned_rev_updated, bin_num_updated = np.loadtxt(file_name_rev_binned_updated).T
	bin_num_unique_updated = np.unique(bin_num_updated)
	main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = cfba.get_output_filenames_and_arrays(dir_name_4, bin_num_unique_updated, wave_rev, par_dict_updated)
	cfba.stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, plot_val=False)
	sys.stdout.flush()
	time.sleep(0.1)
	par_dict_updated['execution_type'] = 'fit'
	cfba.stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, plot_val=False)
	sys.stdout.flush()
	time.sleep(0.1)
	expanded_filename = str(dir_name_4) + str('/') + str(par_dict_updated['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.hdf5')
	expanded_hdr_filename = str(dir_name_4) + str('/') + str(par_dict_updated['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_hdr.fits')
	if (os.path.exists(main_result_cube_filename) and os.path.exists(main_absorption_result_cube_filename) and os.path.exists(main_emission_result_cube_filename)):
		expanded_hdr_filename1, expanded_filename1 = cfba.stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict_updated, wave_rev, file_name_rev_linear, file_name_rev_binned_updated, expanded_hdr_filename, header_rev, expanded_filename)
		code_name = "/show_results.py "
		str_to_execute = str("python3 ") + str(sys.argv[2]) + str(code_name) + str(expanded_hdr_filename1) + str(" ") + str(expanded_filename1)
		bf.print_cust(f"Executing... {str_to_execute}")
		os.system(str_to_execute)
	else:
		bf.print_cust("Please run the fit first")
show_results_button.on_clicked(func_show_results)




fig.suptitle('A2670-3')
plt.show()











quit()







