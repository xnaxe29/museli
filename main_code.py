import sys
import pyfiglet
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
sys.path.append(str(sys.argv[2]) + '/custom_functions')
from custom_functions_by_axe import *
#
ascii_banner = pyfiglet.figlet_format("MUSELI", font='isometric1', justify='center')
print_cust(ascii_banner)
print_cust('Copyright: 크리스, Martin, Axe company')
#
print_cust(str("Successfully imported all modules...." + sys.argv[2]))
#
print_cust(str("Running script: " + sys.argv[0]))
print_cust(str("Current Working Directory: " + sys.argv[1]))
print_cust(str("Base Working Directory: " + sys.argv[2]))
#
#######################OBTAINING_PARAMETER_INFORMATION#######################
#Obtains parameter information from a given parameter file
#If that is not available, it obtains default parameters from code,
#except for the main filename.
#
parameter_file_string_current = str(sys.argv[1]) + str('/parameter_file.dat')
par_dict_init = initial_guess_from_parfile(parameter_file_string_current)
dir_name_6 = str(sys.argv[2]) + str('/custom_functions/basic_files')
par_dict = revise_dictionary(par_dict_init, dir_name_6)
par_dict['quiet'] = True
#
#Telling the code to not print anything on screen (It will still be there in the logfile regardless)
quiet_val = str2bool(par_dict['quiet'])
if not quiet_val: warnings.simplefilter(action='ignore', category=FutureWarning)
if not quiet_val: np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
if not quiet_val: np.seterr(divide = 'ignore')
import custom_functions_by_axe
custom_functions_by_axe.quiet_val = quiet_val

#######################MAKING_DIRECTORIES#######################

dir_name_1 = str(sys.argv[1]) + str("/results")
dir_name_2 = str(sys.argv[1]) + str("/images")
dir_name_3 = str(sys.argv[1]) + str("/post_process_images")
dir_name_4 = str(sys.argv[1]) + str("/data")
check_directory(dir_name_1)
check_directory(dir_name_2)
check_directory(dir_name_3)
check_directory(dir_name_4)


#######################MAKING_DIRECTORIES#######################

#######################EXECUTING_STAGE_ONE#######################

new_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_revised.fits')
if (os.path.exists(new_filename)):
	prompt_1=input(f"File {new_filename} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_1.lower()):
		print_cust(f"deleting...{new_filename}")
		string_for_deleting = str('rm ') + new_filename
		print_cust(string_for_deleting)
		os.system(string_for_deleting)
		stage_one_analysis(par_dict, new_filename, quiet_val=quiet_val)
else:
	stage_one_analysis(par_dict, new_filename, quiet_val=quiet_val)

#######################EXECUTING_STAGE_ONE#######################

#######################EXECUTING_STAGE_TWO#######################

file_name_eq_width = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_eq_width.fits')
if (os.path.exists(file_name_eq_width)):
	prompt_2=input(f"File {file_name_eq_width} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_2.lower()):
		print_cust(f"deleting...{file_name_eq_width}")
		string_for_deleting = str('rm ') + file_name_eq_width
		print_cust(string_for_deleting)
		os.system(string_for_deleting)
		stage_two_analysis(par_dict, new_filename, file_name_eq_width)
else:
		stage_two_analysis(par_dict, new_filename, file_name_eq_width)

#######################EXECUTING_STAGE_TWO#######################

'''
#To obtain EW information from EW data cube
header_ew_original, header_ew_original_err, data_ew_original_file, err_ew_original_file = open_ifu_fits(file_name_eq_width)
fig = plt.figure()
ax1 = fig.add_subplot(121)
im1 = ax1.imshow(data_ew_original_file[-1, :, :], origin='lower', vmin=-10, vmax=10)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')
ax2 = fig.add_subplot(122)
im2 = ax2.imshow(err_ew_original_file[-1, :, :], origin='lower', vmin=-10, vmax=10)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im2, cax=cax, orientation='vertical');
plt.show()
'''

#To obtain information from revised data
print_cust('Loading information from revised data...')
header_rev, header_rev_err, header_rev_cont, header_rev_fwhm, header_rev_velscale, header_rev_sky, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file, velscale_rev_file, sky_rev_file = open_ifu_fits_custom_file(new_filename)
ra_rev, dec_rev, wave_rev = obtain_physical_axis(header_rev)
header_ew_original, header_ew_original_err, data_ew_original_file, err_ew_original_file = open_ifu_fits(file_name_eq_width)
print_cust('information from revised data loaded...')
#print_cust(f'{data_ew_original_file.shape}')
halpha_eq_width_map = data_ew_original_file[int(par_dict['lick_idx_for_halpha']), :, :]
halpha_eq_width_err_map = err_ew_original_file[int(par_dict['lick_idx_for_halpha']), :, :]

halpha_eq_width_map[np.abs(halpha_eq_width_map)>200.0] = 1e-6
halpha_eq_width_err_map[np.abs(halpha_eq_width_map)>200.0] = 1e-5

'''
plt.imshow(halpha_eq_width_map, origin='lower')
plt.colorbar()
plt.show()
quit()
'''


#######################EXECUTING_STAGE_THREE#######################
#par_dict = revise_dictionary(par_dict)
#par_dict['region_of_binning_interest'] = [50, 50, 250, 250]
file_name_rev_linear = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_linear_snr.dat')
#x_array, y_array, ra_array, dec_array, signal1, noise1, signal2, noise2 = np.loadtxt(file_name_rev_linear).T
if (os.path.exists(file_name_rev_linear)):
	prompt_3=input(f"File {file_name_rev_linear} already exists. Overwrite?(y/n) : ")
	if ('y' in prompt_3.lower()):
		print_cust(f"deleting...{file_name_rev_linear}")
		string_for_deleting = str('rm ') + file_name_rev_linear
		print_cust(string_for_deleting)
		os.system(string_for_deleting)
		create_linear_file_for_binning(file_name_rev_linear, par_dict, header_rev, cont_rev_file, data_rev_file, err_rev_file, data_ew_original_file, err_ew_original_file)
else:
	create_linear_file_for_binning(file_name_rev_linear, par_dict, header_rev, cont_rev_file, data_rev_file, err_rev_file, data_ew_original_file, err_ew_original_file)

file_name_rev_binned = stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict)

#######################EXECUTING_STAGE_THREE#######################

'''
par_dict = revise_dictionary(par_dict)
pos_x = 162
pos_y = 168
#pos_x = 150
#pos_y = 150
#pos_x = 105
#pos_y = 250

flux1 = data_rev_file[:, pos_x, pos_y]
err1 = np.sqrt(err_rev_file[:, pos_x, pos_y])
halpha_ew_val1 = halpha_eq_width_map[pos_x, pos_y]
FWHM_gal1 = np.nanmedian(fwhm_rev_file[:, pos_x, pos_y])
flux_sky1 = sky_rev_file
init_cont1 = cont_rev_file[:, pos_x, pos_y]

wave_fit, flux_fit, flux_err_fit, continuum_fit, result_fit, mask_fit, solution_fit, solution_err_fit = main_fitting_func(wave_rev, flux1, err1, par_dict, halpha_ew_val1, FWHM_gal = FWHM_gal1, flux_sky = flux_sky1, init_cont = init_cont1)

plt.errorbar(wave_fit, flux_fit, yerr=flux_err_fit, ds='steps-mid', color='tab:blue', zorder=1, label='data')
plt.plot(wave_fit, result_fit, 'r--', zorder=2, label='fit')
plt.legend()
plt.show()
quit()
'''

#######################EXECUTING_STAGE_FOUR#######################
#file_name_rev_binned = stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict)
#par_dict = revise_dictionary(par_dict)
x_axis_binned_rev, y_axis_binned_rev, bin_num = np.loadtxt(file_name_rev_binned).T
bin_num_unique = np.unique(bin_num)
main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = get_output_filenames_and_arrays(dir_name_4, bin_num_unique, wave_rev, par_dict)

stage_four_execution(par_dict, file_name_rev_linear, file_name_rev_binned, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file)
#######################EXECUTING_STAGE_FOUR#######################

#######################EXECUTING_STAGE_FIVE#######################

expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.hdf5')
expanded_hdr_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_hdr.npy')

expanded_hdr_filename1, expanded_filename1 = stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_binned, expanded_hdr_filename, expanded_filename)

#######################EXECUTING_STAGE_FIVE#######################





snr_vmin_val1 = 0.1
snr_vmax_val1 = 20
physical_axes1 = False
x_array_binned, y_array_binned, snr_revised = get_snr_map_revised(file_name_rev_linear, file_name_rev_binned, par_dict, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)

fig, ax = plt.subplots()
ax.cla()
im = ax.scatter(x_array_binned, y_array_binned, c = snr_revised, cmap='viridis', vmin=snr_vmin_val1, vmax=snr_vmax_val1)
im_cl = add_colorbar_lin(im)

'''
def update(val):
	print_cust("Current Status: ")
	binning_quant_updated = int(sbinning_quant_slider.val)
	print_cust(binning_quant_updated)
	sspectral_smoothing_updated = int(sspectral_smoothing_slider.val)
	print_cust(binning_quant_updated)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	print_cust(binning_quant_updated)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	print_cust(binning_quant_updated)
	radio_work_type_updated = str(radio_work_type.value_selected)
	print_cust(binning_quant_updated)
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	print_cust(binning_quant_updated)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	print_cust(binning_quant_updated)
'''

binning_quant_slider_axes = plt.axes([0.1, 0.01, 0.1, 0.02])
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
radio_binning_type_buttons = RadioButtons(binning_type_buttons_axes, ('None', 'square', 'voronoi', 'advanced'), active=1)
def func_binning_type_buttons(label_binning_type_buttons):
	print_cust(label_binning_type_buttons)
radio_binning_type_buttons.on_clicked(func_binning_type_buttons)

binning_significance_type_buttons_axes = plt.axes([0.01, 0.1, 0.1, 0.1])
dict_binning_significance_type = {'snr', 'ew'}
active_binning_significance_type = list(dict_binning_significance_type).index(str(par_dict['voronoi_snr_type']))
#print (active_binning_significance_type)
#print (list(dict_binning_significance_type)[active_binning_significance_type])
radio_binning_significance_type = RadioButtons(binning_significance_type_buttons_axes, ('snr', 'ew'), active=0)
def func_binning_significance_type(label_binning_significance_type):
	print_cust(label_binning_significance_type)
radio_binning_significance_type.on_clicked(func_binning_significance_type)

#work_type_buttons_axes = plt.axes([0.01, 0.65, 0.1, 0.1])
#radio_work_type = RadioButtons(work_type_buttons_axes, ('snr_map', 'fit'), active=0)
#dict_work_type = {'SNR Map', 'Fit'}
#def func_work_type(label_work_type):
#	print_cust(label_work_type)
#radio_work_type.on_clicked(func_work_type)

fit_type_buttons_axes = plt.axes([0.01, 0.5, 0.1, 0.15])
dict_fit_type = {'auto', 'emission', 'absorption'}
active_dict_fit_type = list(dict_fit_type).index(str(par_dict['execution_fit_type']))
#print (active_dict_fit_type)
#print (list(dict_fit_type)[active_dict_fit_type])
radio_fit_type = RadioButtons(fit_type_buttons_axes, ('auto', 'emission', 'absorption'), active=0)
def func_fit_type(label_fit_type):
	print_cust(label_fit_type)
radio_fit_type.on_clicked(func_fit_type)

emission_fit_type_buttons_axes = plt.axes([0.01, 0.4, 0.1, 0.1])
dict_emission_fit_type = {'custom', 'ppxf'}
active_dict_emission_fit_type = list(dict_emission_fit_type).index(str(par_dict['emission_fit_type']))
#print (active_dict_emission_fit_type)
#print (list(dict_emission_fit_type)[active_dict_emission_fit_type])
radio_emission_fit_type = RadioButtons(emission_fit_type_buttons_axes, ('custom', 'ppxf'), active=0)
def func_emission_fit_type(label_emission_fit_type):
	print_cust(label_emission_fit_type)
radio_emission_fit_type.on_clicked(func_emission_fit_type)

def update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated):
	print_cust('Updating Figure...')
	im = ax.scatter(x_array_binned_updated, y_array_binned_updated, c = snr_revised_updated, cmap='viridis', vmin=snr_vmin_val1, vmax=snr_vmax_val1)
	im_cl = add_colorbar_lin(im)
	print_cust('Figure updated...')


bin_data_buttons_axes = plt.axes([0.01, 0.88, 0.1, 0.02])
bin_data_button = Button(bin_data_buttons_axes, 'Bin Data')
def func_bin_data(label_func_bin_data):
	print_cust('Binning Data...')
	global im_cl
	binning_quant_updated = int(sbinning_quant_slider.val)
	sspectral_smoothing_updated = int(1)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	radio_work_type_updated = str('snr_map')
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	print_cust("Current Status: Spectral Smoothing Factor: {sspectral_smoothing_updated},  Binning Quant: {binning_quant_updated},  Radio Binning Type: {radio_binning_type_updated},  Radio Binning Significance Type: {radio_binning_significance_type_updated},  Radio Work Type: {radio_work_type_updated},  Radio Fit Type: {radio_fit_type_updated},  Radio Emission Fit Type: {radio_emission_fit_type_updated}")
	par_dict_updated = modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated)
	file_name_rev_binned_updated = stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict_updated)
	x_axis_binned_rev_updated, y_axis_binned_rev_updated, bin_num_updated = np.loadtxt(file_name_rev_binned_updated).T
	bin_num_unique_updated = np.unique(bin_num_updated)
	main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = get_output_filenames_and_arrays(dir_name_4, bin_num_unique_updated, wave_rev, par_dict_updated)
	stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file)
	sys.stdout.flush()
	time.sleep(0.1)
	x_array_binned_updated, y_array_binned_updated, snr_revised_updated = get_snr_map_revised(file_name_rev_linear, file_name_rev_binned_updated, par_dict_updated, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)
	ax.cla()
	#im_cl.remove()
	update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated)
bin_data_button.on_clicked(func_bin_data)
fig.canvas.draw_idle()

fit_data_buttons_axes = plt.axes([0.01, 0.82, 0.1, 0.02])
fit_data_button = Button(fit_data_buttons_axes, 'Fit Data')
def func_fit_data(event):
	print_cust('Binning Data...')
	global im_cl
	binning_quant_updated = int(sbinning_quant_slider.val)
	sspectral_smoothing_updated = int(1)
	radio_binning_type_updated = str(radio_binning_type_buttons.value_selected)
	radio_binning_significance_type_updated = str(radio_binning_significance_type.value_selected)
	radio_work_type_updated = str('snr_map')
	radio_fit_type_updated = str(radio_fit_type.value_selected)
	radio_emission_fit_type_updated = str(radio_emission_fit_type.value_selected)
	print_cust("Current Status: Spectral Smoothing Factor: {sspectral_smoothing_updated},  Binning Quant: {binning_quant_updated},  Radio Binning Type: {radio_binning_type_updated},  Radio Binning Significance Type: {radio_binning_significance_type_updated},  Radio Work Type: {radio_work_type_updated},  Radio Fit Type: {radio_fit_type_updated},  Radio Emission Fit Type: {radio_emission_fit_type_updated}")
	par_dict_updated = modify_dictionary(par_dict, binning_quant_updated, radio_binning_type_updated, radio_binning_significance_type_updated, radio_work_type_updated, radio_fit_type_updated, radio_emission_fit_type_updated)
	file_name_rev_binned_updated = stage_three_analysis(file_name_rev_linear, dir_name_4, par_dict_updated)
	x_axis_binned_rev_updated, y_axis_binned_rev_updated, bin_num_updated = np.loadtxt(file_name_rev_binned_updated).T
	bin_num_unique_updated = np.unique(bin_num_updated)
	main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, main_result_cube, main_absorption_result_cube, main_emission_result_cube = get_output_filenames_and_arrays(dir_name_4, bin_num_unique_updated, wave_rev, par_dict_updated)
	stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file)
	sys.stdout.flush()
	time.sleep(0.1)
	par_dict_updated['execution_type'] = 'fit'
	stage_four_execution(par_dict_updated, file_name_rev_linear, file_name_rev_binned_updated, dir_name_4, wave_rev, halpha_eq_width_map, sky_rev_file, data_rev_file, err_rev_file, cont_rev_file, fwhm_rev_file)
	x_array_binned_updated, y_array_binned_updated, snr_revised_updated = get_snr_map_revised(file_name_rev_linear, file_name_rev_binned_updated, par_dict_updated, physical_axes = physical_axes1, snr_vmin_val = snr_vmin_val1, snr_vmax_val = snr_vmax_val1, quiet_val=False, return_map=True)
	ax.cla()
	#im_cl.remove()
	update_fig(ax, x_array_binned_updated, y_array_binned_updated, snr_revised_updated)
fit_data_button.on_clicked(func_fit_data)
fig.canvas.draw_idle()




show_results_buttons_axes = plt.axes([0.01, 0.76, 0.1, 0.02])
show_results_button = Button(show_results_buttons_axes, 'Show Results')
def func_show_results(event):
	print_cust('Show Results...')
	expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.hdf5')
	#expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.fits')
	expanded_hdr_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_hdr.npy')
	if (os.path.exists(expanded_hdr_filename) and os.path.exists(expanded_filename)):
		expanded_hdr_filename1, expanded_filename1 = stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_binned, expanded_hdr_filename, expanded_filename)
		code_name = "/show_results.py "
		str_to_execute = str("python3 ") + str(sys.argv[2]) + str(code_name) + str(expanded_hdr_filename1) + str(" ") + str(expanded_filename1)
		print_cust(f"Executing... {str_to_execute}")
		os.system(str_to_execute)
	else:
		print_cust("Please run the fit first")
show_results_button.on_clicked(func_show_results)




fig.suptitle('A2670-3')
plt.show()











quit()












#######################EXECUTING_STAGE_FIVE#######################

expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.hdf5')
#expanded_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_data.fits')
expanded_hdr_filename = str(dir_name_4) + str('/') + str(par_dict['muse_data_filename']).split("/")[-1].replace('.fits', '') + str('_temporary_extended_hdr.npy')

expanded_hdr_filename1, expanded_filename1 = stage_five_analysis(main_result_cube_filename, main_absorption_result_cube_filename, main_emission_result_cube_filename, par_dict, wave_rev, file_name_rev_binned, expanded_hdr_filename, expanded_filename)

print_cust("Loading expanded data")
start_time6 = time.time()
header_information_absorption, header_information_err_absorption, header_information_emission, header_information_err_emission = np.load(expanded_hdr_filename1, allow_pickle=True)
with h5py.File(expanded_filename1, 'r') as hf:
	amp_length = hf["amp_length"]
	wave_new2 = hf["wave_new"]
	expanded_data_cube = hf["expanded_data_cube"][:]
	expanded_information_cube_absorption = hf["expanded_information_cube_absorption"][:]
	expanded_information_cube_emission = hf["expanded_information_cube_emission"][:]

print_cust(f"Loading took {float(time.time() - start_time6)} seconds ---")

start_time4 = time.time()
map = np.nanmean(expanded_data_cube[:, :, :, 0], axis=(0))
plt.imshow(map, origin='lower')
plt.title('Flux')
plt.colorbar()
plt.show()

map = expanded_information_cube_absorption[0, :, :, 0]
plt.imshow(map, origin='lower')
plt.title('Absorption')
plt.colorbar()
plt.show()

map = expanded_information_cube_emission[0, :, :, 0]
plt.imshow(map, origin='lower')
plt.title('Emission')
plt.colorbar()
plt.show()

print_cust(f"Plotting took {float(time.time() - start_time4)} seconds ---")
quit()
#######################EXECUTING_STAGE_FIVE#######################









quit()



















