import numpy as np
import matplotlib.pyplot as plt
from basic_functions import *
quite_val = False

##################################MAKE_REVISED_PPXF_FIGURE##################################
def make_figure_ppxf(fig_x_size, fig_y_size, dpi_val, lam_gal, galaxy, noise, residuals, bestfit_solution_array, st_age_unique, st_mass_unique, st_lum_unique, size_of_font, vel_window, mean_wave_list, lick_index_species, x_pos, y_pos, **kwargs):

	savefig = kwargs.get('savefig', 'display') # Save figure?
	good_index = kwargs.get('good_index', np.ones([len(lam_gal)], dtype=np.bool)) # Goodpixels in fit
	str_fit_val = kwargs.get('str_fit_val', str('Fit')) # Goodpixels in fit
	quiet = kwargs.get('quiet', False) # Goodpixels in fit

	if not (savefig=='display'):
		fig = plt.figure(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig = plt.figure()
		size_of_font = size_of_font/2

	#gs = gridspec.GridSpec(4, 5, height_ratios=[5, 2, 2, 2])
	gs = gridspec.GridSpec(4, len(mean_wave_list), height_ratios=[5, 2, 2, 2])
	ax0 = plt.subplot(gs[0, :])
	ax1 = plt.subplot(gs[1, :])
	ax2 = plt.subplot(gs[3, :])
	ax3 = []
	for i in range(len(mean_wave_list)):
		ax3 = np.append(ax3, plt.subplot(gs[2, i]))

	#lam_gal = lam_gal / (1.+0.07527116015236746)
	ax0.plot(lam_gal, galaxy, color='tab:blue', drawstyle='steps-mid', alpha=0.4, label=r'data', zorder=1)
	ax0.plot(lam_gal[good_index], galaxy[good_index], color='tab:blue', drawstyle='steps-mid', alpha=0.8, label=r'goodpixel', zorder=2)
	ax0.plot(lam_gal, bestfit_solution_array, 'r--', label=str(str_fit_val), zorder=3)
	ax0.legend(loc='best', fontsize=size_of_font)
	ax0.set_xticks([])
	ax0.set_ylabel(r"Relative Flux", fontsize=size_of_font)
	ax1.plot(lam_gal, residuals, 'k.', markersize=size_of_font/4, alpha=0.2, zorder=4)
	ax1.plot(lam_gal[good_index], residuals[good_index], 'k.', markersize=size_of_font/4, alpha=0.8, zorder=5)
	ax1.axhline(0.0, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(-3, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(3, color='green', ls='dashed', alpha=0.5)
	ax1.text(6500, -4.0, s=r'-3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.text(6500, 3.2, s=r'3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.set_ylim(-5,5)
	ax1.set_xlabel(r"Rest Wavelength ($\rm \AA$)", fontsize=size_of_font)
	ax1.set_ylabel(r"Residuals", fontsize=size_of_font)
	label_x_ticks = np.array(["<0.1", "0.1-0.5", "0.5-1.0", "1.0-5.0", "5.0-10.0", ">10.0"])
	ax2.plot(label_x_ticks, np.log10(st_mass_unique), 'bo', markersize=size_of_font, label='Stellar Mass')
	ax2.plot(label_x_ticks, np.log10(st_lum_unique), 'r*', markersize=size_of_font, label='Stellar Luminosity')
	ax2.set_xlabel(r"Age (Gyr)", fontsize=size_of_font)
	ax2.set_ylabel(r"Relative Weights (log)", fontsize=size_of_font)
	ax2.legend(loc='best', fontsize=size_of_font)
	
	for i in range(len(mean_wave_list)):
		vel31 = vel_prof(lam_gal, mean_wave_list[i])
		idx_start = np.searchsorted(vel31, -vel_window)
		idx_end = np.searchsorted(vel31, vel_window)
		ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
		ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
		ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
		ax3[i].set_title(lick_index_species[i], fontsize=size_of_font)
		ax3[i].tick_params(axis='both', labelsize=size_of_font)

	ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)
	ax0.tick_params(axis='both', labelsize=size_of_font)
	ax1.tick_params(axis='both', labelsize=size_of_font)
	ax2.tick_params(axis='both', labelsize=size_of_font)
	ax0.tick_params(size=size_of_font)
	ax1.tick_params(size=size_of_font)
	ax2.tick_params(size=size_of_font)
	fig_name_cust = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure.pdf")
	if not (savefig=='display'):
		fig.tight_layout()
		plt.savefig(fig_name_cust)
		print_cust(f'{fig_name_cust}, saved', quiet_val=quiet)
		plt.close()
	else:
		plt.show()
##################################MAKE_REVISED_PPXF_FIGURE##################################


##################################MAKE_FIGURE##################################
def make_figure(x_slot, y_slot, fig_x_size, fig_y_size, dpi_val, size_of_font, qso_name, fig_name, name_list_sorted, center_list_sorted, wave_fit, flux_fit, err_fit, cont_array_fitted, result_fitted, reddening_array_fit, redshift_val, vel_window, center_list_init, comments_on_balmer, number_of_narrow_components_init, number_of_wide_components_init, amp_array_rev, center_array, sigma_array, quiet_val=False):
	f, axarr = plt.subplots(x_slot, y_slot, sharex=False, sharey=False, figsize=((fig_x_size), (fig_y_size)), dpi=dpi_val)
	f.subplots_adjust(wspace=0.4, hspace=0.4)
	for u in range(x_slot):
		for v in range(y_slot):
			if (len(center_list_sorted) > ((x_slot*v)+u)):
				test_x3 = vel_prof(wave_fit/(1.+redshift_val), center_list_sorted[(x_slot*v)+u])
				idx1 = find_nearest_idx(test_x3,-vel_window-100)
				idx2 = find_nearest_idx(test_x3,vel_window+100)
				if (u>1):
					axarr[u,v].errorbar(wave_fit[idx1:idx2], flux_fit[idx1:idx2], yerr=err_fit[idx1:idx2], color='black', drawstyle='steps-mid', zorder=1)
					axarr[u,v].plot(wave_fit[idx1:idx2], result_fitted[idx1:idx2], color='red', ls='--', zorder=2)
				else:
					axarr[v].errorbar(wave_fit[idx1:idx2], flux_fit[idx1:idx2], yerr=err_fit[idx1:idx2], color='black', drawstyle='steps-mid', zorder=1)
					axarr[v].plot(wave_fit[idx1:idx2], result_fitted[idx1:idx2], color='red', ls='--', zorder=2)
				count_amp = 0
				for z in range(len(center_list_init)):
					centre = center_list_init[z]
					vel_array = vel_prof(wave_fit/(1.+redshift_val), centre)
					count = 0
					if (comments_on_balmer[z]):
						for k in range(number_of_narrow_components_init):
							#group_prof += gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)

							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
						for l in range(number_of_wide_components_init):
							#group_prof += gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'm--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='magenta', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'm--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='magenta', alpha=0.5)
							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
					else:
						for m in range(number_of_narrow_components_init):
							#group_prof += gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count])
							prof_fit = gaus_prof_vel(vel_array, amp_array_rev[count_amp], center_array[count], sigma_array[count]) + cont_array_fitted
							prof_fit_adv = func_6(wave_fit/(1.+redshift_val), prof_fit, *reddening_array_fit)
							if (u>1):
								axarr[u,v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[u,v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							else:
								axarr[v].plot(wave_fit[idx1:idx2]*(1.+redshift_val), prof_fit_adv[idx1:idx2], 'g--')
								axarr[v].axvline(wave_prof(center_array[count], centre)*(1.+redshift_val), ls='dashed', color='green', alpha=0.5)
							prof_fit = []
							prof_fit_adv = []
							count+=1
							count_amp+=1
				if (u>1):
					axarr[u,v].set_title(name_list_sorted[(x_slot*v)+u], fontsize=size_of_font)
					axarr[u,v].set_xlim(wave_fit[idx1], wave_fit[idx2])
					axarr[u,v].tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font/2, width=2, colors='k')
					axarr[u,v].tick_params(axis = 'both', which = 'minor', direction='in', length=size_of_font/4, width=1, colors='k')
					axarr[u,v].tick_params(axis='both', labelsize=size_of_font)
				else:
					axarr[v].set_title(name_list_sorted[(x_slot*v)+u], fontsize=size_of_font)
					axarr[v].set_xlim(wave_fit[idx1], wave_fit[idx2])
					axarr[v].tick_params(axis = 'both', which = 'major', direction='in', length=size_of_font/2, width=2, colors='k')
					axarr[v].tick_params(axis = 'both', which = 'minor', direction='in', length=size_of_font/4, width=1, colors='k')
					axarr[v].tick_params(axis='both', labelsize=size_of_font)
			else:
				if (u>1):
					axarr[u,v].set_visible(False)
				else:
					axarr[v].set_visible(False)
	f.text(0.48, 0.94, qso_name, fontsize=1.2*size_of_font)
	f.text(0.52, 0.06, r'Wavelength ($\rm \AA$)', ha='center', va='center', fontsize=1.2*size_of_font)
	f.text(0.06, 0.5, r'Flux (10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ AA$^{-1}$)', ha='center', va='center', rotation='vertical', fontsize=1.2*size_of_font)
	f.patch.set_linewidth(2)
	f.patch.set_edgecolor('black')
	plt.savefig(fig_name, edgecolor=f.get_edgecolor())
	print_cust(f'{fig_name} saved...', quiet_val=quiet_val)

def plot_quick_figure(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init):
	fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots()
	ax_secondary_1d_data.errorbar(wave, flux, yerr=err, color='tab:blue', alpha=0.5, label='Data', zorder=1)
	ax_secondary_1d_data.errorbar(wave_fitted, flux_fitted, yerr=flux_err_fitted, color='tab:blue', label='Fitted Region', zorder=2)
	ax_secondary_1d_data.plot(wave_fitted, result_fitted, 'r.-', label='Total Fit', zorder=4)
	ax_secondary_1d_data.plot(wave_fitted, continuum_fitted, 'g--', label='Continuum', zorder=6)
	count_amp = 0
	for j in range(len(center_list_init)):
		centre = center_list_init[j]*(1.+redshift_val)
		vel_array = vel_prof(wave_fitted, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='magenta', alpha=0.5)
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
	ax_secondary_1d_data.set_xlabel(r"Wavelength ($\rm \AA$)")
	ax_secondary_1d_data.set_ylabel(r"Relative Flux")
	plt.legend()
	plt.show()

def save_quick_figure_rev(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, residuals, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init, fig_x_size, fig_y_size, dpi_val, x_pos, y_pos, size_of_font, str_fit, mean_wave_list, lick_index_species, vel_window, save_fig='display', **kwargs):
	quite_val = kwargs.get('quite_val', False)  # length of amplitude array to be fitted
	if not (save_fig=='display'):
		fig = plt.figure(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig = plt.figure()

	gs = gridspec.GridSpec(3, len(mean_wave_list), height_ratios=[5, 2, 2])
	ax3 = []
	for i in range(len(mean_wave_list)):
		ax3 = np.append(ax3, plt.subplot(gs[2, i]))
	ax0 = plt.subplot(gs[0, :])
	ax1 = plt.subplot(gs[1, :], sharex=ax0)
	lam_gal = wave_fitted / (1.+redshift_val)
	wave = wave / (1.+redshift_val)
	galaxy = flux_fitted
	noise = flux_err_fitted
	bestfit_solution_array = result_fitted
	ax0.errorbar(wave, flux, yerr=err, color='tab:blue', drawstyle='steps-mid', alpha=0.5, label=r'data', zorder=1)
	ax0.plot(lam_gal, flux_fitted, color='tab:blue', drawstyle='steps-mid', alpha=0.8, label=r'to fit', zorder=2)
	ax0.plot(lam_gal, result_fitted, 'r--', label=str(str_fit), zorder=3)
	ax0.legend(loc='best', fontsize=size_of_font)
	ax0.set_ylabel(r"Relative Flux", fontsize=size_of_font)
	ax0.set_xlim(np.nanmin(lam_gal)-10., np.nanmax(lam_gal)+10.)
	ax0.set_ylim(np.nanmin(flux_fitted)-5., np.nanmax(flux_fitted)+5.)
	ax1.plot(lam_gal, residuals, 'k.', markersize=size_of_font/4, alpha=0.8, zorder=4)
	ax1.axhline(0.0, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(-3, color='green', ls='dashed', alpha=0.5)
	ax1.axhline(3, color='green', ls='dashed', alpha=0.5)
	ax1.text(6500, -4.0, s=r'-3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.text(6500, 3.2, s=r'3$\rm \sigma$', fontsize=size_of_font, color='green')
	ax1.set_ylim(-5,5)
	ax1.set_xlabel(r"Rest Wavelength ($\rm \AA$)", fontsize=size_of_font)
	ax1.set_ylabel(r"Residuals", fontsize=size_of_font)
	for i in range(len(mean_wave_list)):
		vel31 = vel_prof(lam_gal, mean_wave_list[i])
		idx_start = np.searchsorted(vel31, -vel_window)
		idx_end = np.searchsorted(vel31, vel_window)
		ax3[i].errorbar(vel31[idx_start:idx_end], galaxy[idx_start:idx_end], yerr=noise[idx_start:idx_end], ds='steps-mid', color='tab:blue')
		ax3[i].plot(vel31[idx_start:idx_end], bestfit_solution_array[idx_start:idx_end], 'r--')
		ax3[i].set_xlabel(r"Rel. Velocity (kms$^{-1}$)", fontsize=size_of_font)
		ax3[i].set_title(lick_index_species[i], fontsize=size_of_font)
		ax3[i].tick_params(axis='both', labelsize=size_of_font)

	ax0.tick_params(axis='both', labelsize=size_of_font)
	ax1.tick_params(axis='both', labelsize=size_of_font)
	ax3[0].set_ylabel(r"Rel. Flux", fontsize=size_of_font)
	ax0.tick_params(size=size_of_font)
	ax1.tick_params(size=size_of_font)
	fig_name_extension = kwargs.get('fig_name_extension', '')  # length of amplitude array to be fitted
	fig_name_cust_orig = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure") + str(fig_name_extension) + str(".pdf")
	if not (save_fig=='display'):
		fig.tight_layout()
		plt.savefig(fig_name_cust_orig)
		print_cust(f'{fig_name_cust_orig} saved', quiet_val=quite_val)
	else:
		plt.show()

	plt.close('all')

def save_quick_figure(wave, flux, err, wave_fitted, flux_fitted, flux_err_fitted, result_fitted, continuum_fitted, center_list_init, redshift_val, comments_on_balmer, center_array, number_of_narrow_components_init, number_of_wide_components_init, fig_x_size, fig_y_size, dpi_val, x_pos, y_pos, save_fig='display', quite_val=False):
	if not (save_fig=='display'):
		fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots(figsize=[fig_x_size, fig_y_size], dpi=dpi_val)
	else:
		fig_secondary_1d_data, (ax_secondary_1d_data) = plt.subplots()

	ax_secondary_1d_data.errorbar(wave, flux, yerr=err, color='tab:blue', alpha=0.5, label='Data', zorder=1)
	ax_secondary_1d_data.errorbar(wave_fitted, flux_fitted, yerr=flux_err_fitted, color='tab:blue', label='Fitted Region', zorder=2)
	ax_secondary_1d_data.plot(wave_fitted, result_fitted, 'r.-', label='Total Fit', zorder=4)
	ax_secondary_1d_data.plot(wave_fitted, continuum_fitted, 'g--', label='Continuum', zorder=6)
	count_amp = 0
	for j in range(len(center_list_init)):
		centre = center_list_init[j]*(1.+redshift_val)
		vel_array = vel_prof(wave_fitted, centre)
		count = 0
		if (comments_on_balmer[j]):
			for k in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
			for l in range(number_of_wide_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='magenta', alpha=0.5)
				count+=1
				count_amp+=1
		else:
			for m in range(number_of_narrow_components_init):
				ax_secondary_1d_data.axvline(wave_prof(center_array[count], centre), ls='dashed', color='green', alpha=0.5)
				count+=1
				count_amp+=1
	ax_secondary_1d_data.set_xlabel(r"Wavelength ($\rm \AA$)")
	ax_secondary_1d_data.set_ylabel(r"Relative Flux")
	fig_secondary_1d_data.legend()
	fig_name_cust = str("./all_pixel_images/Pos_") + str(int(x_pos)) + str("_") + str(int(y_pos)) + str("_fit_figure.pdf")
	if not (save_fig=='display'):
		plt.savefig(fig_name_cust)
		print_cust(f'{fig_name_cust} saved', quiet_val=quite_val)
	else:
		plt.show()
	plt.close('all')
##################################MAKE_FIGURE##################################




####################GET_SNR_MAP_FOR_BINNED_DATA####################
def get_snr_map_revised(file_name_rev_linear, file_name_rev_binned, par_dict, physical_axes = False, snr_vmin_val = 0.1, snr_vmax_val = 20, quiet_val=False, return_map=False):
	print_cust('Plotting SNR Map...', quiet_val=quiet_val)
	assert ('voronoi_snr_type' in par_dict), f"voronoi_snr_type required..."
	y_array, x_array, ra_array, dec_array, signal1, noise1, signal2, noise2 = np.loadtxt(file_name_rev_linear).T
	x_array_binned, y_array_binned, bin_num = np.loadtxt(file_name_rev_binned).T
	bin_unique = np.unique(bin_num)
	signal_total = np.zeros_like(bin_num)
	noise_total = np.ones_like(bin_num)
	if 'ew' in par_dict['voronoi_snr_type']:
		signal = signal1
		noise = noise1
	else:
		signal = signal2
		noise = noise2
	snr = signal / noise
	for i in range(len(bin_unique)):
		mask = np.isin(bin_num, bin_unique[i])
		signal_tmp = np.nansum(signal[mask])
		noise_tmp = np.sqrt(np.nansum(noise[mask]**2))
		signal_total[mask] = signal_tmp
		noise_total[mask] = noise_tmp
		signal_tmp = 0.0
		noise_tmp = 1.0
	snr_revised = signal_total / noise_total
	if (return_map):
		return (x_array_binned, y_array_binned, snr_revised)
	else:
		print_cust('Plotting SNR Map...', quiet_val=quiet_val)
		fig2, axs2 = plt.subplots(2, figsize=(6,10), sharex=True, sharey=True)
		if (physical_axes):
			im2 = axs2[0].scatter(dec, ra, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
		else:
			im2 = axs2[0].scatter(x_array, y_array, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
		axs2[0].set_title("SNR")
		add_colorbar(im2)
		#snr_revised = np.log10(snr_revised)
		im22 = axs2[1].scatter(x_array_binned, y_array_binned, c = snr_revised, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val*100)
		#im2 = axs[1].scatter(x_array_binned, y_array_binned, c = np.log10(snr_revised), cmap='viridis')
		axs2[1].set_title("Binning")
		add_colorbar(im22)
		#plt.pause(1)
		plt.show()
		print_cust('SNR Map Plotted...', quiet_val=quiet_val)
		#return None
		#continue
####################GET_SNR_MAP_FOR_BINNED_DATA####################





