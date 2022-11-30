import sys
import time
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons, RectangleSelector
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pyfiglet
from astropy.cosmology import Planck15 as LCDM
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import h5py
import numpy as np
quiet_val = False
#
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

c = 299792.458 # Speed in Light in Km/s
c_kms = 299792.458 # Speed in Light in Km/s
#
ascii_banner = pyfiglet.figlet_format("MUSELI DEUX", font='isometric1', justify='center')
print_cust(ascii_banner)
print_cust('Copyright: 크리스, Martin & Axe corp.')
#
print_cust(str("Successfully imported all modules...." + sys.argv[2]))
#
#Setting up cosmology
cosmo = FlatLambdaCDM(H0=67.8 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.308)

#Some random formatting code
underline = '\033[4m'
end_formatting = end_format = reset = '\033[0m'

expanded_hdr_filename2 = str(sys.argv[1])
expanded_filename2 = str(sys.argv[2])

print_cust("Loading expanded data")
start_time6 = time.time()
header_information_absorption, header_information_err_absorption, header_information_emission, header_information_err_emission = np.load(expanded_hdr_filename2, allow_pickle=True)
with h5py.File(expanded_filename2, 'r') as hf:
	amp_length = hf["amp_length"]
	wave_new2 = hf["wave_new"][:]
	expanded_data_cube = hf["expanded_data_cube"][:]
	expanded_information_cube_absorption = hf["expanded_information_cube_absorption"][:]
	expanded_information_cube_emission = hf["expanded_information_cube_emission"][:]
print_cust(f"Loading took {float(time.time() - start_time6)} seconds ---")

header_information_absorption = np.array(np.hstack(header_information_absorption))
header_information_emission = np.array(np.hstack(header_information_emission))

map = expanded_information_cube_emission[0, :, :, 0]

fig, ax = plt.subplots()
ax.cla()
im = ax.imshow(map, origin='lower', cmap='viridis')
im_cl = add_colorbar(im)

def update(val):
	global im_cl
	ax.cla()
	im_cl.remove()
	emission_val = int(emission_slider.val)
	absorption_val = int(absorption_slider.val)
	significance_val = int(significance_slider.val)
	radio_plot_val = int(dict_radio_plot_buttons[radio_plot_buttons.value_selected])
	if (radio_plot_val==5):
		plot_file = expanded_information_cube_emission
		data_updated = plot_file[emission_val, :, :, 0]
		err_updated = plot_file[emission_val, :, :, 1]
		err_updated[err_updated==0] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		plot_data[ratio < significance_val] = np.nan
		im = ax.imshow(plot_data, origin='lower', cmap='viridis')
		ax.title.set_text(str(header_information_emission[emission_val]))
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl = add_colorbar(im)
		else:
			im_cl = add_colorbar_lin(im)
	elif (radio_plot_val==6):
		plot_file = expanded_information_cube_absorption
		data_updated = plot_file[absorption_val, :, :, 0]
		err_updated = plot_file[absorption_val, :, :, 1]
		err_updated[err_updated==0] = 1.
		ratio = data_updated / err_updated
		plot_data = data_updated
		plot_data[ratio < significance_val] = np.nan
		im = ax.imshow(plot_data, origin='lower', cmap='viridis')
		ax.title.set_text(str(header_information_absorption[absorption_val]))
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl = add_colorbar(im)
		else:
			im_cl = add_colorbar_lin(im)
	else:
		plot_file = expanded_data_cube
		data_updated = np.nanmean(expanded_data_cube[:, :, :, radio_plot_val], axis=(0))
		im = ax.imshow(data_updated, origin='lower', cmap='viridis')
		if ('log' in str(radio_loglin_buttons.value_selected)):
			im_cl = add_colorbar(im)
		else:
			im_cl = add_colorbar_lin(im)

fig.canvas.draw_idle()

emission_slider_axes = plt.axes([0.1, 0.02, 0.2, 0.05])
emission_slider = Slider(emission_slider_axes, 'Emission', 0, len(header_information_emission)-1, valinit=0, valfmt="%i")
emission_slider.on_changed(update)

absorption_slider_axes = plt.axes([0.4, 0.02, 0.2, 0.05])
absorption_slider = Slider(absorption_slider_axes, 'Absorption', 0, len(header_information_absorption)-1, valinit=0, valfmt="%i")
absorption_slider.on_changed(update)

significance_slider_axes = plt.axes([0.7, 0.02, 0.2, 0.05])
significance_slider = Slider(significance_slider_axes, 'Significance', 1, 10, valinit=3, valfmt="%i")
significance_slider.on_changed(update)

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
	flux = expanded_data_cube[:, ix, iy, 0]
	flux_err = expanded_data_cube[:, ix, iy, 1]
	cont = expanded_data_cube[:, ix, iy, 2]
	fit = expanded_data_cube[:, ix, iy, 3]
	fig_lin, ax_lin = plt.subplots()
	ax_lin.errorbar(wave_new2[mask], flux[mask], yerr=flux_err[mask], ds='steps-mid', color='tab:blue', label='data', zorder=1)
	ax_lin.plot(wave_new2[mask], cont[mask], ls='--', color='green', label='cont/stellar', zorder=2)
	ax_lin.plot(wave_new2[mask], fit[mask], ls='--', color='tab:red', label='fit', zorder=3)
	ax_lin.legend(loc="best")
	fig_lin.show()
plot_one_d_button.on_clicked(func_plot_one_d)

make_figure_buttons_axes = plt.axes([0.02, 0.82, 0.15, 0.05])
make_figure_button = Button(make_figure_buttons_axes, 'Make 2D')
def func_make_figure(event):
	print ('Making Figure...')
make_figure_button.on_clicked(func_make_figure)



plot_buttons_axes = plt.axes([0.02, 0.4, 0.15, 0.4])
radio_plot_buttons = RadioButtons(plot_buttons_axes, ('Flux', 'Err', 'Stellar', 'Result', 'Mask', 'Emission', 'Absorption'), active=5)
dict_radio_plot_buttons = {'Flux': 0, 'Err': 1, 'Stellar': 2, 'Result': 3, 'Mask': 4, 'Emission': 5, 'Absorption': 6}
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










