#import wx
#import matplotlib.backends.backend_wxagg
import sys
import time
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pyfiglet
from astropy.cosmology import Planck15 as LCDM
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import h5py
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.io.fits as pyfits
from plotbin.plot_velfield import plot_velfield
from matplotlib.patches import Ellipse
import matplotlib.ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#from kinemetry import kinemetry
#import kinemetry as kin
import time
from os import path
from progress.bar import IncrementalBar
quiet_val = False

#----------------------------------------------------------------------------
def vel_prof(x, centre):
    xnew = c_kms * ((x-centre)/x)
    return (xnew)


#----------------------------------------------------------------------------
def print_cust(print_str, **kwargs):
	sleep_time = kwargs.get('sleep_time', 0.1)  # Sleep Time
	quite_val = kwargs.get('quiet_val', False)  # Quite Val
	if not quite_val:
		print ("\n")
		print (print_str)
		print ("\n")
		sys.stdout.flush()
		time.sleep(sleep_time)


#----------------------------------------------------------------------------
def save_linear_file(header_information_absorption, header_information_emission, expanded_information_cube_absorption, expanded_information_cube_emission, ra_array, dec_array, filename='test_linear_kin_file.dat'):
    #num, xbin, ybin, velbin, er_velbin, sigbin, er_sigbin, h3bin, er_h3bin, h4bin, er_h4bin = np.genfromtxt(file, unpack=True)
	#print (header_information_absorption)
	#print (ra_array)
	#print (dec_array)
	idx_v = np.where(header_information_absorption=='V-1')[0][0]
	idx_sigma = np.where(header_information_absorption=='sigma1')[0][0]
	idx_v_emission = np.where(header_information_emission=='V-1')[0][0]
	idx_sigma_emission = np.where(header_information_emission=='sigma-1')[0][0]
	idx_h3 = np.where(header_information_absorption=='h3-1')[0][0]
	idx_h4 = np.where(header_information_absorption=='h4-1')[0][0]
	#print (idx_v, idx_sigma, idx_h3, idx_h4)
	vel_map = expanded_information_cube_absorption[idx_v, :, :, 0]
	er_vel_map = expanded_information_cube_absorption[idx_v, :, :, 1]
	sigma_map = expanded_information_cube_absorption[idx_sigma, :, :, 0]
	er_sigma_map = expanded_information_cube_absorption[idx_sigma, :, :, 1]
	vel_map_em = expanded_information_cube_emission[idx_v_emission, :, :, 0]
	er_vel_map_em = expanded_information_cube_emission[idx_v_emission, :, :, 1]
	sigma_map_em = expanded_information_cube_emission[idx_sigma_emission, :, :, 0]
	er_sigma_map_em = expanded_information_cube_emission[idx_sigma_emission, :, :, 1]
	h3_map = expanded_information_cube_absorption[idx_h3, :, :, 0]
	er_h3_map = expanded_information_cube_absorption[idx_h3, :, :, 1]
	h4_map = expanded_information_cube_absorption[idx_h4, :, :, 0]
	er_h4_map = expanded_information_cube_absorption[idx_h4, :, :, 1]
	chi_sq_abs = expanded_information_cube_absorption[5, :, :, 1]
	count=1
	linear_array = np.chararray([int(len(ra_array)*len(dec_array))+1, 16], itemsize=100)
	linear_array[0, :] = np.array(['#Count', 'RA', 'Dec', 'V', 'del_V', 'sigma', 'del_sigma', 'V_em', 'del_V_em', 'sigma_em', 'del_sigma_em', 'h3', 'del_h3', 'h4', 'del_h4', 'chi_sq'])
	bar = IncrementalBar('Countdown', max = int(len(ra_array)*len(dec_array)))
	for i in range(len(ra_array)):
		for j in range(len(dec_array)):
			bar.next()
			#linear_array.extend(np.transpose(np.array([int(count), float(ra_array[i]), float(dec_array[j]), vel_map[j,i], er_vel_map[j,i], sigma_map[j,i], er_sigma_map[j,i], h3_map[j,i], er_h3_map[j,i], h4_map[j,i], er_h4_map[j,i]])))
			linear_array[count, 0] = int(count)
			linear_array[count, 1] = float(ra_array[i])
			linear_array[count, 2] = float(dec_array[j])
			linear_array[count, 3] = vel_map[j,i]
			linear_array[count, 4] = er_vel_map[j,i]
			linear_array[count, 5] = sigma_map[j,i]
			linear_array[count, 6] = er_sigma_map[j,i]
			linear_array[count, 7] = vel_map_em[j,i]
			linear_array[count, 8] = er_vel_map_em[j,i]
			linear_array[count, 9] = sigma_map_em[j,i]
			linear_array[count, 10] = er_sigma_map_em[j,i]
			linear_array[count, 11] = h3_map[j,i]
			linear_array[count, 12] = er_h3_map[j,i]
			linear_array[count, 13] = h4_map[j,i]
			linear_array[count, 14] = er_h4_map[j,i]
			linear_array[count, 15] = chi_sq_abs[j,i]
			count+=1

	linear_array1 = linear_array.astype(np.str_)
	np.savetxt(filename, linear_array1, fmt='%s', delimiter=',', encoding='latin1')
	print_cust(f'Linear array: {filename} saved')




c = 299792.458 # Speed in Light in Km/s
c_kms = 299792.458 # Speed in Light in Km/s
#
ascii_banner = pyfiglet.figlet_format("MUSELI DEUX", font='isometric1', justify='center')
print_cust(ascii_banner)
print_cust('Authors: Axe, Saulder, Martin, and Bag')
# 크리스, Martin & Axe corp.')
#
print_cust(str("Successfully imported all modules...." + sys.argv[2]))
#
#Setting up cosmology
cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.308)

#Some random formatting code
underline = '\033[4m'
end_formatting = end_format = reset = '\033[0m'

codename = str(sys.argv[0])
test5 = codename.split('/')
test6 = str(test5[-1])
test7 = codename.replace(test6, "")
test8 = test7 + str("/custom_functions")
sys.path.append(test8)
import basic_functions as bf
import plotting_functions as pf
import handling_fits_files as hff

expanded_hdr_filename2 = str(sys.argv[1])
test1 = expanded_hdr_filename2.split('/')
test2 = str(test1[-2])+str('/')+str(test1[-1])
test3 = expanded_hdr_filename2.replace(test2, "")
test4 = test3 + str("/post_process_images")
figname1d = test4 + str("/fig1d_")
figname2d = test4 + str("/fig2d_")
expanded_filename2 = str(sys.argv[2])

print_cust("Loading expanded data")
start_time6 = time.time()
#header_information_absorption, header_information_err_absorption, header_information_emission, header_information_err_emission = np.load(expanded_hdr_filename2, allow_pickle=True)

with fits.open(expanded_hdr_filename2) as hdul:
	header_original = hdul[0].header
	header_information_absorption = hdul[1].data['header_information_absorption']
	header_information_err_absorption = hdul[1].data['header_information_err_absorption']
	header_information_emission = hdul[1].data['header_information_emission']
	header_information_err_emission = hdul[1].data['header_information_err_emission']

header_information_absorption = header_information_absorption[header_information_absorption != '']
header_information_err_absorption = header_information_err_absorption[header_information_err_absorption != '']
header_information_emission = header_information_emission[header_information_emission != '']
header_information_err_emission = header_information_err_emission[header_information_err_emission != '']

#print (header_information_absorption)
header_information_absorption[10:16] = ['Age <0.1', '0.1 < Age < 0.5', '0.5 < Age < 1.0', '1.0 < Age < 5.0', '5.0 < Age < 10.0', 'Age > 10.0']
header_information_absorption[16:22] = ['Age <0.1', '0.1 < Age < 0.5', '0.5 < Age < 1.0', '1.0 < Age < 5.0', '5.0 < Age < 10.0', 'Age > 10.0']
header_information_absorption[22:28] = ['Age <0.1', '0.1 < Age < 0.5', '0.5 < Age < 1.0', '1.0 < Age < 5.0', '5.0 < Age < 10.0', 'Age > 10.0']
#print (header_information_absorption)
#quit()

#print (header_information_emission)


with h5py.File(expanded_filename2, 'r') as hf:
	#header_original = hf["header"]
	amp_length = hf["amp_length"]
	wave_new2 = hf["wave_new"][:]
	x_axis_unique = hf["ra"][:]
	y_axis_unique = hf["dec"][:]
	sky_spectrum = hf["sky"][:]
	expanded_data_cube = hf["expanded_data_cube"][:]
	expanded_information_cube_absorption = hf["expanded_information_cube_absorption"][:]
	expanded_information_cube_emission = hf["expanded_information_cube_emission"][:]

#'''
significance_number = np.zeros([expanded_data_cube.shape[1], expanded_data_cube.shape[2]])
bar = IncrementalBar('Countdown', max = int(expanded_data_cube.shape[1]*expanded_data_cube.shape[2]))
if (len(sky_spectrum)==int(expanded_data_cube.shape[0])):
	for i in range(expanded_data_cube.shape[1]):
		for j in range(expanded_data_cube.shape[2]):
			bar.next()
			significance_number[i, j] = float(np.nansum(np.abs(expanded_data_cube[:, i, j, 3] - sky_spectrum)) / int(expanded_data_cube.shape[0]))
else:
	for i in range(expanded_data_cube.shape[1]):
		for j in range(expanded_data_cube.shape[2]):
			bar.next()
			significance_number[i, j] = float(np.nansum(np.abs(expanded_data_cube[:, i, j, 3] / expanded_data_cube[:, i, j, 1])) / int(expanded_data_cube.shape[0]))
#'''



muse_wcs = WCS(header_original).celestial
header_information_absorption = np.array(np.hstack(header_information_absorption.flatten()))
header_information_emission = np.array(np.hstack(header_information_emission.flatten()))
#print (header_information_emission)
#print (len(header_information_emission))
#print (expanded_information_cube_emission.shape)
#print (len(header_information_absorption))
#print (expanded_information_cube_absorption.shape)
#print (x_axis_unique)
#print (y_axis_unique)
#quit()
ra_array = x_axis_unique
dec_array = y_axis_unique
#test4 = test3 + str("/post_process_images")


#tmp_linear_filename = test3 + str("/data/tmp_linear_data_rev.dat")
#save_linear_file(header_information_absorption, header_information_emission, expanded_information_cube_absorption, expanded_information_cube_emission, ra_array, dec_array, filename=tmp_linear_filename)


map = expanded_information_cube_emission[0, :, :, 0]
vel_kinemetry_figname = test3 + str("/post_process_images/vel_kinemetry.pdf")
sigma_kinemetry_figname = test3 + str("/post_process_images/sigma_kinemetry.pdf")
h3_kinemetry_figname = test3 + str("/post_process_images/h3_kinemetry.pdf")
h4_kinemetry_figname = test3 + str("/post_process_images/h4_kinemetry.pdf")

#flux_map = data_updated = np.nanmean(expanded_data_cube[:, :, :, 0], axis=(0))
#flux_err_map = data_updated = np.nanmean(expanded_data_cube[:, :, :, 1], axis=(0))
#stellar_map = data_updated = np.nanmean(expanded_data_cube[:, :, :, 2], axis=(0))
#fitted_map = data_updated = np.nanmean(expanded_data_cube[:, :, :, 3], axis=(0))
#mask_map = data_updated = np.nanmean(expanded_data_cube[:, :, :, 4], axis=(0))
#total_3d_data_maps = np.array([flux_map, flux_err_map, stellar_map, fitted_map, mask_map])
#print (total_3d_data_maps.shape)
print_cust(f"Loading took {float(time.time() - start_time6)} seconds ---")

#'''
fig, ax = plt.subplots()
ax.cla()
im = ax.imshow(map, origin='lower', cmap='viridis')
im_cl = bf.add_colorbar(im)
ax.title.set_text(str(header_information_emission[0]))
#'''

plot_file_em = expanded_information_cube_emission
chi_sq_em = expanded_information_cube_emission[-2, :, :, 0]
#significance_number = expanded_information_cube_emission[-1, :, :, 1]
plot_file_abs = expanded_information_cube_absorption
chi_sq_abs = expanded_information_cube_absorption[5, :, :, 1]
plot_file_data = expanded_data_cube
chi_sq_tot_min = np.nanmin(np.append(chi_sq_em, chi_sq_abs))
chi_sq_tot_max = np.nanmax(np.append(chi_sq_em, chi_sq_abs))

'''
plt.imshow(chi_sq_em, origin='lower')
plt.colorbar()
plt.show()
plt.imshow(chi_sq_abs, origin='lower')
plt.colorbar()
plt.show()
quit()
'''

def update(val):
	global im_cl
	ax.cla()
	im_cl.remove()
	emission_val = int(emission_slider.val)
	absorption_val = int(absorption_slider.val)
	significance_val = float(significance_slider.val)
	chi_sq_limit_val = float(chi_sq_slider.val)
	radio_plot_val = int(dict_radio_plot_buttons[radio_plot_buttons.value_selected])
	if (radio_plot_val):
		data_updated = plot_file_abs[absorption_val, :, :, 0]
		err_updated = plot_file_abs[absorption_val, :, :, 1]
		err_updated[err_updated==0] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		plot_data[significance_number < significance_val] = np.nan
		plot_data[chi_sq_abs > chi_sq_limit_val] = np.nan
		im = ax.imshow(plot_data, origin='lower', cmap='viridis')
		ax.title.set_text(str(header_information_absorption[absorption_val]))
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl = bf.add_colorbar(im)
		else:
			im_cl = bf.add_colorbar_lin(im)
	else:
		data_updated = plot_file_em[emission_val, :, :, 0]
		err_updated = plot_file_em[emission_val, :, :, 1]
		err_updated[err_updated==0] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		plot_data[significance_number < significance_val] = np.nan
		plot_data[chi_sq_em > chi_sq_limit_val] = np.nan
		im = ax.imshow(plot_data, origin='lower', cmap='viridis')
		ax.title.set_text(str(header_information_emission[emission_val]))
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl = bf.add_colorbar(im)
		else:
			im_cl = bf.add_colorbar_lin(im)

fig.canvas.draw_idle()

emission_slider_axes = plt.axes([0.1, 0.02, 0.2, 0.05])
emission_slider = Slider(emission_slider_axes, 'Emission', 0, len(header_information_emission)-1, valinit=0, valfmt="%i")
emission_slider.on_changed(update)

absorption_slider_axes = plt.axes([0.4, 0.02, 0.2, 0.05])
absorption_slider = Slider(absorption_slider_axes, 'Absorption', 0, len(header_information_absorption)-1, valinit=0, valfmt="%i")
absorption_slider.on_changed(update)

significance_slider_axes = plt.axes([0.7, 0.02, 0.2, 0.05])
significance_slider = Slider(significance_slider_axes, 'Significance', np.nanmin(significance_number), np.nanmax(significance_number), valinit=np.nanmin(significance_number), valfmt="%2.2f")
significance_slider.on_changed(update)

chi_sq_slider_axes = plt.axes([0.7, 0.9, 0.2, 0.05])
chi_sq_slider = Slider(chi_sq_slider_axes, 'Chi_sq', chi_sq_tot_min, chi_sq_tot_max, valinit=chi_sq_tot_max, valfmt="%2.2f")
chi_sq_slider.on_changed(update)


def toggle_selector(event):
	print(' Key pressed.')
	if event.key in ['H', 'h']:
		global ix, iy
		ix, iy = int(event.xdata), int(event.ydata)
		print_cust(f"Selected: X - {ix}, Y - {iy}")
plt.connect('key_press_event', toggle_selector)


plot_one_d_buttons_axes = plt.axes([0.02, 0.88, 0.15, 0.05])
plot_one_d_button = Button(plot_one_d_buttons_axes, 'Plot 1D')
def func_plot_one_d(event):
	print ('Plot 1D...')
	global ix, iy
	mask = (wave_new2!=0)
	flux = expanded_data_cube[:, iy, ix, 0]
	flux_err = expanded_data_cube[:, iy, ix, 1]
	cont = expanded_data_cube[:, iy, ix, 2]
	fit = expanded_data_cube[:, iy, ix, 3]
	equivalent_width = expanded_information_cube_emission[-1, iy, ix, 0]
	#new_array = np.array(np.vstack([wave_new2, flux, flux_err, cont, fit]))
	#np.savetxt('test.dat', new_array)
	x_pos = int(iy)
	y_pos = int(ix)
	lick_index_species = np.array(['Hbeta_o', 'Mg1', 'NaD', 'Fe6189', 'Halpha'])
	mean_wave_list = np.array([4864.8735, 5120.541667, 5898.791667, 6170.666667, 6540.166667])
	em_index_species = np.array(['Hbeta', '[OIII]4958', '[OIII]5007', '[NII]6547', 'Halpha', '[NII]6583', '[SII]6716', '[SII]6730'])
	em_mean_wave_list = np.array([4861.32, 4958.83, 5006.77, 6547.96, 6562.8, 6583.34, 6716.31, 6730.68])
	lam_gal = wave_new2[mask]
	galaxy = flux[mask]
	noise = flux_err[mask]
	bestfit_solution_array = fit[mask]
	residuals = (fit[mask] - flux[mask]) / flux_err[mask]
	st_age_unique1 = plot_file_abs[10:16, iy, ix, 0]
	st_mass_unique1 = plot_file_abs[16:22, iy, ix, 0]
	st_lum_unique1 = plot_file_abs[22:28, iy, ix, 0]
	redshift_val = 0.07527116015236746
	#redshift_val = 0.050200146
	figname1d_rev = figname1d + str(int(iy)) + str("_") + str(int(ix)) + str(".pdf")
	pf.ppxf_figure(lam_gal, galaxy, noise, residuals, bestfit_solution_array, st_age_unique=st_age_unique1, st_mass_unique=st_mass_unique1, st_lum_unique=st_lum_unique1, figname=figname1d_rev, redshift=redshift_val, plot_type='emission')
plot_one_d_button.on_clicked(func_plot_one_d)


make_figure_buttons_axes = plt.axes([0.02, 0.82, 0.15, 0.05])
make_figure_button = Button(make_figure_buttons_axes, 'Make 2D')
def func_make_figure(event):
	print_cust('Making Figure...')
	fig2 = plt.figure(2, figsize=[10,8])
	ax_save = fig2.add_subplot(111, projection=muse_wcs)
	emission_val = int(emission_slider.val)
	absorption_val = int(absorption_slider.val)
	significance_val = float(significance_slider.val)
	chi_sq_limit_val = float(chi_sq_slider.val)
	radio_plot_val = int(dict_radio_plot_buttons[radio_plot_buttons.value_selected])
	dx = (x_axis_unique[1]-x_axis_unique[0])/2.
	dy = (y_axis_unique[1]-y_axis_unique[0])/2.
	extent1 = [x_axis_unique[0]-dx, x_axis_unique[-1]+dx, y_axis_unique[0]-dy, y_axis_unique[-1]+dy]

	if (radio_plot_val):
		data_updated = plot_file_abs[absorption_val, :, :, 0]
		err_updated = plot_file_abs[absorption_val, :, :, 1]
		err_updated[err_updated<=0.] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		title_name = str(header_information_absorption[absorption_val])
		plot_data[significance_number < significance_val] = np.nan
		plot_data[chi_sq_abs > chi_sq_limit_val] = np.nan
		if ('v' in title_name.lower()):
			plot_data = plot_data - np.nanmedian(plot_data)
		im_save = ax_save.imshow(plot_data, origin='lower', cmap='viridis', vmax=1.0)
		ax_save.set_title(title_name, fontsize=20, fontweight='bold')
		figname2d_2 = figname2d + str("_abs_")
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl_save = bf.add_colorbar(im_save)
		else:
			im_cl_save = bf.add_colorbar_lin(im_save)
	else:
		data_updated = plot_file_em[emission_val, :, :, 0]
		err_updated = plot_file_em[emission_val, :, :, 1]
		err_updated[err_updated<=0.] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		title_name = str(header_information_emission[emission_val])
		plot_data[significance_number < significance_val] = np.nan
		plot_data[chi_sq_em > chi_sq_limit_val] = np.nan
		if ('v' in title_name.lower()):
			plot_data = plot_data - np.nanmedian(plot_data)
		im_save = ax_save.imshow(plot_data, origin='lower', cmap='viridis')
		ax_save.title.set_text(title_name)
		figname2d_2 = figname2d + str("_em_")
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl_save = bf.add_colorbar(im_save)
		else:
			im_cl_save = bf.add_colorbar_lin(im_save)

	figname2d_rev = figname2d_2 + title_name + str(".pdf")
	fig2.savefig(figname2d_rev, dpi=100)
	plt.close(fig2)
	print_cust(f'Figure {figname2d_rev} Saved...')
make_figure_button.on_clicked(func_make_figure)



'''
pafit_buttons_axes = plt.axes([0.02, 0.1, 0.15, 0.05])
pafit_button = Button(pafit_buttons_axes, 'PAFIT')
def func_pafit(event):
	print_cust('PAFIT...')
	num, xbin, ybin, velbin, er_velbin, sigbin, er_sigbin, h3bin, er_h3bin, h4bin, er_h4bin, chi_sq_lin = np.genfromtxt(tmp_linear_filename, unpack=True)
	print_cust('kinemetry on velocity map + fit centre')
	x0_val = xbin[162]
	y0_val = ybin[168]
	scale_val = np.sqrt( (xbin[162]-xbin[161])**2 + (ybin[168]-ybin[167])**2 )*10
	data = velbin
	data_err = er_velbin
	vel_map_centre_kinemetry(xbin, ybin, data, data_err, x0_val, y0_val, scale_val, figname=vel_kinemetry_figname, title_name='A2670-3')
	#print_cust('kinemetry on velocity dispersion map, limited PA, Q, fixed centre')
	#sig_map_kinemetry(xbin, ybin, data, data_err, x0_val, y0_val, scale_val, ntrm_val=8, pa_min=-180., pa_max=180., q_min=0., q_max=1., figname=sigma_kinemetry_figname, title_name='A2670-3')
	#print_cust('kinemetry on h3 map, limited PA, Q, fixed centre')
	#sig_map_kinemetry(xbin, ybin, data, data_err, x0_val, y0_val, scale_val, ntrm_val=8, pa_min=-180., pa_max=180., q_min=0., q_max=1., figname=h3_kinemetry_figname, title_name='A2670-3')
	#print_cust('kinemetry on h4 map, limited PA, Q, fixed centre')
	#sig_map_kinemetry(xbin, ybin, data, data_err, x0_val, y0_val, scale_val, ntrm_val=8, pa_min=-180., pa_max=180., q_min=0., q_max=1., figname=h3_kinemetry_figname, title_name='A2670-3')
	#vel_map_centre_kinemetry(figname=vel_kinemetry_figname)
	#sig_map_kinemetry(figname=sigma_kinemetry_figname)
	#sig_map_kinemetry(figname=h3_kinemetry_figname)
	#sig_map_kinemetry(figname=h4_kinemetry_figname)

pafit_button.on_clicked(func_pafit)
'''








plot_buttons_axes = plt.axes([0.02, 0.5, 0.15, 0.2])
#radio_plot_buttons = RadioButtons(plot_buttons_axes, ('Flux', 'Err', 'Stellar', 'Result', 'Mask', 'Emission', 'Absorption'), active=5)
#dict_radio_plot_buttons = {'Flux': 0, 'Err': 1, 'Stellar': 2, 'Result': 3, 'Mask': 4, 'Emission': 5, 'Absorption': 6}
radio_plot_buttons = RadioButtons(plot_buttons_axes, ('Emission', 'Absorption'), active=0)
dict_radio_plot_buttons = {'Emission': 0, 'Absorption': 1}
def func_plot_buttons(label_plot_buttons):
	print (label_plot_buttons)
radio_plot_buttons.on_clicked(func_plot_buttons)


loglin_buttons_axes = plt.axes([0.02, 0.2, 0.15, 0.15])
radio_loglin_buttons = RadioButtons(loglin_buttons_axes, ('lin', 'log'), active=0)
dict_radio_loglin_buttons = {'lin': 0, 'log': 1}
def func_loglin_buttons(label_loglin_buttons):
	print (label_loglin_buttons)
radio_loglin_buttons.on_clicked(func_loglin_buttons)









fig.suptitle('A2670-3')
plt.show()










