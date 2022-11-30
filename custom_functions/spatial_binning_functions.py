from vorbin.voronoi_2d_binning import voronoi_2d_binning
import numpy as np
from progress.bar import IncrementalBar
from basic_functions import *
from handling_fits_files import *
quite_val = False



##############################NO_SPATIAL_BINNING##################################
def no_spatial_binning_func_main(input_filename, output_filename, physical_axes=False, quiet_val=False):
	print_cust(f"Loading linear data file: {input_filename} - ", quiet_val=quiet_val)
	y, x, ra, dec, signal1, noise1, signal2, noise2 = np.loadtxt(input_filename).T
	if (physical_axes):
		x_lin = ra
		y_lin = dec
	else:
		x_lin = x
		y_lin = y
	count = 1
	final_pos_x = np.array([])
	final_pos_y = np.array([])
	final_pos_count = np.array([])
	for i in range(0, len(x_lin), 1):
		final_pos_x = np.append(final_pos_x, x_lin[i])
		final_pos_y = np.append(final_pos_y, y_lin[i])
		final_pos_count = np.append(final_pos_count, count)
		count+=1
	print_cust(f"Number of bins created: {count}", quiet_val=quiet_val)
	final_array_to_save = np.transpose(np.array([final_pos_x, final_pos_y, final_pos_count]))
	np.savetxt(output_filename, final_array_to_save, fmt='%s', delimiter='\t')
	print_cust(f"Binning Saved in file: {output_filename}", quiet_val=quiet_val)
def no_spatial_binning_func(input_filename, output_filename, par_dict, physical_axes_val=False, quiet_val=False):
	if (os.path.isfile(output_filename)):
		save_prompt6= input(f"Binning file: {output_filename} already exists. Overwrite?(y/n) : ")
		if (save_prompt6=='y'):
			print_cust('Setting linear file as binned file', quiet_val=quiet_val)
			no_spatial_binning_func_main(input_filename, output_filename, physical_axes=physical_axes_val)
			print_cust('Saved...', quiet_val=quiet_val)
		else:
			print_cust('Exiting binning process....', quiet_val=quiet_val)
	else:
		print_cust('Setting linear file as binned file', quiet_val=quiet_val)
		no_spatial_binning_func_main(input_filename, output_filename, physical_axes=physical_axes_val)
		print_cust('Saved...', quiet_val=quiet_val)
##############################NO_SPATIAL_BINNING##################################


##############################SQUARE_BINNING##################################
def square_binning_func_main(input_filename, output_filename, sqaure_unit, physical_axes=False, quiet_val=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	print_cust(f"Loading linear data file: {input_filename} - ", quiet_val=quiet_val)
	y1, x1, ra1, dec1, signal11, noise11, signal12, noise12 = np.loadtxt(input_filename).T

	region_of_binning_interest = np.array(region_of_binning_interest, dtype=np.int)
	if (region_of_binning_interest[2]!=-1):
		if ('logical' in region_of_binning_interest_type):
			mask = (region_of_binning_interest[0]<=x1) & (x1<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=y1) & (y1<=region_of_binning_interest[3])
		else:
			mask = (region_of_binning_interest[0]<=dec1_u) & (dec1_u<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=ra1_u) & (ra1_u<=region_of_binning_interest[3])
	else:
		mask = np.ones([len(x1)], dtype=np.bool)
	
	x = x1[mask]; y = y1[mask]; ra = ra1[mask]; dec = dec1[mask]; signal1 = signal11[mask]; noise1 = noise11[mask]; signal2 = signal12[mask]; noise2 = noise12[mask]

	if (physical_axes):
		x_lin = ra
		y_lin = dec
		x_lin2 = ra1[~mask]
		y_lin2 = dec1[~mask]
	else:
		x_lin = x
		y_lin = y
		x_lin2 = x1[~mask]
		y_lin2 = y1[~mask]
		
	count = 1
	final_pos_x = x_lin
	final_pos_y = y_lin
	final_pos_count = np.zeros_like(x_lin)
	x_lin_unique = np.sort(np.unique(x_lin))
	y_lin_unique = np.sort(np.unique(y_lin))
	for i in range(0, len(x_lin_unique), int(2*sqaure_unit)):
		for j in range(0, len(y_lin_unique), int(2*sqaure_unit)):
			mask = np.isin(x_lin, x_lin_unique[i-sqaure_unit:i+sqaure_unit]) & np.isin(y_lin, y_lin_unique[j-sqaure_unit:j+sqaure_unit])
			#final_pos_x = np.append(final_pos_x, x_lin[mask])
			#final_pos_y = np.append(final_pos_y, y_lin[mask])
			#final_pos_count = np.append(final_pos_count, np.full([len(x_lin[mask])], count))
			final_pos_count[mask] = count
			count+=1
			
	binned_array2 = np.column_stack([final_pos_x, final_pos_y, final_pos_count])
	binned_array1 = np.column_stack([x_lin2, y_lin2, [int(count)]*len(x_lin2)])
	binned_array3 = np.row_stack([binned_array2, binned_array1])
	count+=1
	print_cust(f"Number of bins created: {count}", quiet_val=quiet_val)
	#final_array_to_save = np.transpose(np.array([final_pos_x, final_pos_y, final_pos_count]))
	final_array_to_save = binned_array3
	np.savetxt(output_filename, final_array_to_save, fmt='%s', delimiter='\t')
	print_cust(f"Binning Saved in file: {output_filename}", quiet_val=quiet_val)

def square_binning_func(input_filename, output_filename, par_dict, physical_axes_val=False, quiet_val=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	sqaure_unit = int(par_dict['binning_quant'])
	if (os.path.isfile(output_filename)):
		save_prompt6= input(f"Binning file: {output_filename} already exists. Overwrite?(y/n) : ")
		if (save_prompt6=='y'):
			square_binning_func_main(input_filename, output_filename, sqaure_unit, physical_axes=physical_axes_val, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
			print_cust('Saved...', quiet_val=quiet_val)
		else:
			print_cust('Exiting binning process....', quiet_val=quiet_val)
	else:
		square_binning_func_main(input_filename, output_filename, sqaure_unit, physical_axes=physical_axes_val, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
		print_cust('Saved...', quiet_val=quiet_val)
##############################SQUARE_BINNING##################################


##############################ADVANCED_BINNING##################################
#Function created with the help of Dr. Christoph Saulder (13/10/2022)
#####Function to get the angle in degrees between 0,0 and a point in the x,y plane
def calcAngleDegrees(x, y):
	return np.arctan2(y, x) * 180.0 / np.pi

#####Function to get the region convered by an angle, using center - center_x, center_y; points in x,y plane - position_x_all, position_y_all; and between the start and end angles - angle_start, angle_end(in degrees)
def region_inside_angle(center_x, center_y, position_x_all, position_y_all, angle_start, angle_end):
	angle = calcAngleDegrees((position_x_all-np.full([len(position_x_all)], center_x)), (position_y_all-np.full([len(position_y_all)], center_y)))
	mask_custom = (angle_start<angle) & (angle<=angle_end)
	return (mask_custom)

#####Function to get an ellipse, using center - h,k; points in x,y plane - x,y; tilted by an angle - theta (in degrees) and with - a,b semi minor, major axes
def ellipse_function_cust(x,y,h,k,a,b,theta):
	p1 = ((x-h)*np.cos(np.radians(theta)) + (y-k)*np.sin(np.radians(theta)))**2 / a**2
	p2 = ((x-h)*np.sin(np.radians(theta)) + (y-k)*np.cos(np.radians(theta)))**2 / b**2
	return (p1+p2)

#####Function to get the region convered by an ellipse, using center - x0,y0; points in x,y plane - x_array,y_array; tilted by an angle - theta (in degrees) and with - radius r as semi minor axis and r*2 as semi major axis
def points_inside_ellipse(x_array, y_array, radius, theta=0, x0=0, y0=0):
	mask_custom = (ellipse_function_cust(x_array,y_array,x0,y0,radius,radius*2,theta) <= 1.0)
	return mask_custom

#####Complex function combining the region of an ellipse using function - 'points_inside_ellipse' and arc region within a specific angle between angle_start_for_arc and angle_end_for_arc(in degrees)
def region_selector(x_axis, y_axis, center_x, center_y, ellipse_radius_end, ellipse_angle, angle_end_for_arc, ellipse_radius_init=0, angle_start_for_arc=0):
	if (ellipse_radius_init):
		ellipse_start_mask = points_inside_ellipse(x_axis, y_axis, ellipse_radius_init, theta=ellipse_angle, x0=center_x, y0=center_y)
		ellipse_end_mask = points_inside_ellipse(x_axis, y_axis, ellipse_radius_end, theta=ellipse_angle, x0=center_x, y0=center_y)
		doughnut_mask = (np.invert(ellipse_start_mask)) & (ellipse_end_mask)

	else:
		doughnut_mask = points_inside_ellipse(x_axis, y_axis, ellipse_radius_end, theta=ellipse_angle, x0=center_x, y0=center_y)

	regioned_mask = region_inside_angle(center_x, center_y, x_axis, y_axis, angle_start_for_arc, angle_end_for_arc)
	mask_final = (doughnut_mask & regioned_mask)
	return (mask_final)

def advanced_binning_func_basic(input_filename, output_filename, ellipse_angle, bin_chopping_angle, binning_radius_list, center_coordinates, physical_axes=False, quiet_val=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	print_cust(f"Loading linear data file: {input_filename} - ", quiet_val=quiet_val)
	y1, x1, ra1, dec1, signal11, noise11, signal12, noise12 = np.loadtxt(input_filename).T

	region_of_binning_interest = np.array(region_of_binning_interest, dtype=np.int)
	if (region_of_binning_interest[2]!=-1):
		if ('logical' in region_of_binning_interest_type):
			mask = (region_of_binning_interest[0]<=x1) & (x1<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=y1) & (y1<=region_of_binning_interest[3])
		else:
			mask = (region_of_binning_interest[0]<=dec1_u) & (dec1_u<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=ra1_u) & (ra1_u<=region_of_binning_interest[3])
	else:
		mask = np.ones([len(x1)], dtype=np.bool)
	
	x = x1[mask]; y = y1[mask]; ra = ra1[mask]; dec = dec1[mask]; signal1 = signal11[mask]; noise1 = noise11[mask]; signal2 = signal12[mask]; noise2 = noise12[mask]

	if (physical_axes):
		x_lin = ra
		y_lin = dec
		x_lin2 = ra1[~mask]
		y_lin2 = dec1[~mask]
	else:
		x_lin = x
		y_lin = y
		x_lin2 = x1[~mask]
		y_lin2 = y1[~mask]
		
	count = 1
	final_pos_x = np.array([])
	final_pos_y = np.array([])
	final_pos_count = np.array([])
	center_x = center_coordinates[0]
	center_y = center_coordinates[1]
	if (np.abs(float(ellipse_angle))<180):
		ellipse_angle = 180. - float(ellipse_angle)
	if (bin_chopping_angle!=0 and binning_radius_list):
		angle_array = np.arange(-180, 180.0001, int(bin_chopping_angle))
		radius_array = np.array(binning_radius_list)
		for i in range(len(radius_array)-1):
			for j in range(len(angle_array)-1):
				#print_cust('ellipse_angle', quiet_val=quiet_val)
				mask = region_selector(y_lin, x_lin, center_x, center_y, radius_array[i+1], ellipse_angle, angle_array[j+1], ellipse_radius_init=radius_array[i], angle_start_for_arc=angle_array[j])
				final_pos_x = np.append(final_pos_x, x_lin[mask])
				final_pos_y = np.append(final_pos_y, y_lin[mask])
				final_pos_count = np.append(final_pos_count, np.full([len(x_lin[mask])], count))
				count+=1
	elif (bin_chopping_angle!=0 and not binning_radius_list):
		angle_array = np.arange(-180, 180.0001, int(bin_chopping_angle))
		for j in range(len(angle_array)-1):
			mask = region_selector(y_lin, x_lin, center_x, center_y, 1000, ellipse_angle, angle_array[j+1], ellipse_radius_init=0, angle_start_for_arc=angle_array[j])
			final_pos_x = np.append(final_pos_x, x_lin[mask])
			final_pos_y = np.append(final_pos_y, y_lin[mask])
			final_pos_count = np.append(final_pos_count, np.full([len(x_lin[mask])], count))
			count+=1
	elif (bin_chopping_angle==0 and binning_radius_list):
		radius_array = np.array(binning_radius_list)
		for i in range(len(radius_array)-1):
			mask = region_selector(y_lin, x_lin, center_x, center_y, radius_array[i+1], ellipse_angle, 180.0001, ellipse_radius_init=radius_array[i], angle_start_for_arc=-180)
			final_pos_x = np.append(final_pos_x, x_lin[mask])
			final_pos_y = np.append(final_pos_y, y_lin[mask])
			final_pos_count = np.append(final_pos_count, np.full([len(x_lin[mask])], count))
			count+=1
	else:
		print_cust(f"Warning! No specific bin chipping angle or radius list for ellipse defined.", quiet_val=quiet_val)
		mask = region_selector(y_lin, x_lin, center_x, center_y, 1000, ellipse_angle, 180.0001, ellipse_radius_init=0, angle_start_for_arc=-180)
		final_pos_x = np.append(final_pos_x, x_lin[mask])
		final_pos_y = np.append(final_pos_y, y_lin[mask])
		final_pos_count = np.append(final_pos_count, np.full([len(x_lin[mask])], count))
		count+=1

	binned_array2 = np.column_stack([final_pos_x, final_pos_y, final_pos_count])
	binned_array1 = np.column_stack([x_lin2, y_lin2, [int(count)]*len(x_lin2)])
	binned_array3 = np.row_stack([binned_array2, binned_array1])
	count+=1
	print_cust(f"Number of bins created: {count}", quiet_val=quiet_val)
	#final_array_to_save = np.transpose(np.array([final_pos_x, final_pos_y, final_pos_count]))
	final_array_to_save = binned_array3
	np.savetxt(output_filename, final_array_to_save, fmt='%s', delimiter='\t')
	print_cust(f"Binning Saved in file: {output_filename}", quiet_val=quiet_val)

def advanced_binning_func(input_filename, output_filename, par_dict, physical_axes_val=False, quiet_val=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	if ('binning_radius_list' not in par_dict):
		binning_radius_list = False
	elif (par_dict['binning_radius_list']=='[]'):
		binning_radius_list = False
	else:
		binning_radius_list = par_dict['binning_radius_list']

	if ('ellipse_angle_for_advanced_binning' in par_dict):
		ellipse_angle = par_dict['ellipse_angle_for_advanced_binning']
	else:
		ellipse_angle = 0

	if ('center_coordinate_for_advanced_binning' in par_dict):
		center_coordinates = par_dict['center_coordinate_for_advanced_binning']
	else:
		center_coordinates = [0, 0]

	bin_chopping_angle = par_dict['binning_quant']
	if (os.path.isfile(output_filename)):
		save_prompt5= input(f"Binning file: {output_filename} already exists. Overwrite?(y/n) : ")
		if (save_prompt5=='y'):
			advanced_binning_func_basic(input_filename, output_filename, ellipse_angle, bin_chopping_angle, binning_radius_list, center_coordinates, physical_axes=physical_axes_val, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
			print_cust('Saved...', quiet_val=quiet_val)
		else:
			print_cust('Exiting binning process....', quiet_val=quiet_val)
	else:
		advanced_binning_func_basic(input_filename, output_filename, ellipse_angle, bin_chopping_angle, binning_radius_list, center_coordinates, physical_axes=physical_axes_val, region_of_binning_interest=region_of_binning_interest, region_of_binning_interest_type=region_of_binning_interest_type)
		print_cust('Saved...', quiet_val=quiet_val)
##############################ADVANCED_BINNING##################################

######################VORONOI_BINNING##################################
def voronoi_binning_func(input_filename, output_filename, targetSN, snr_type, plot_val=False, quiet_val=True, physical_axes=False, region_of_binning_interest=[0, 0, -1, -1], region_of_binning_interest_type='logical'):
	y1, x1, ra1, dec1, signal11, noise11, signal12, noise12 = np.loadtxt(input_filename).T

	region_of_binning_interest = np.array(region_of_binning_interest, dtype=np.int)
	if (region_of_binning_interest[2]!=-1):
		if ('logical' in region_of_binning_interest_type):
			mask = (region_of_binning_interest[0]<=x1) & (x1<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=y1) & (y1<=region_of_binning_interest[3])
		else:
			mask = (region_of_binning_interest[0]<=dec1_u) & (dec1_u<=region_of_binning_interest[2]) & (region_of_binning_interest[1]<=ra1_u) & (ra1_u<=region_of_binning_interest[3])
	else:
		mask = np.ones([len(x1)], dtype=np.bool)
	
	x = x1[mask]; y = y1[mask]; ra = ra1[mask]; dec = dec1[mask]; signal1 = signal11[mask]; noise1 = noise11[mask]; signal2 = signal12[mask]; noise2 = noise12[mask]

	print_cust(f"Target SNR: {targetSN}", quiet_val=quiet_val)
	PIXSIZE = np.sqrt( (x[1]-x[0])**2. + (y[1]-y[0])**2. )
	if ('ew' in snr_type.lower()):
		signal = signal2
		noise = noise2
	else:
		signal = signal1
		noise = noise1
	snr = signal/noise
	print_cust(f"SNR/EW Limits: {np.min(snr)} to {np.max(snr)}", quiet_val=quiet_val)
	if (os.path.isfile(output_filename)):
		save_prompt4= input(f"Binning file: {output_filename} already exists. Overwrite?(y/n) : ")
		if (save_prompt4=='y'):
			if (physical_axes):
				binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(dec, ra, signal, noise, targetSN, plot=plot_val, quiet=quiet_val, pixelsize=PIXSIZE)
				binned_array2 = np.column_stack([dec, ra, binNum])
				binned_array1 = np.column_stack([dec[~mask], ra[~mask], [int(max(binNum)+1)]*len(ra[~mask])])
				binned_array3 = np.row_stack([binned_array2, binned_array1])
				np.savetxt(output_filename, binned_array3, fmt=b'%10.6f %10.6f %8i')
				print_cust('Saved...', quiet_val=quiet_val)
			else:
				binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(x, y, signal, noise, targetSN, plot=plot_val, quiet=quiet_val, pixelsize=PIXSIZE)
				binned_array2 = np.column_stack([x, y, binNum])
				binned_array1 = np.column_stack([x[~mask], y[~mask], [int(max(binNum)+1)]*len(x[~mask])])
				binned_array3 = np.row_stack([binned_array2, binned_array1])
				np.savetxt(output_filename, binned_array3, fmt=b'%10.6f %10.6f %8i')
				print_cust('Saved...', quiet_val=quiet_val)
		else:
			print_cust('Exiting binning process....', quiet_val=quiet_val)
	else:
		if (physical_axes):
			binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(dec, ra, signal, noise, targetSN, plot=plot_val, quiet=quiet_val, pixelsize=PIXSIZE)
			binned_array2 = np.column_stack([dec, ra, binNum])
			binned_array1 = np.column_stack([dec[~mask], ra[~mask], [int(max(binNum)+1)]*len(ra[~mask])])
			binned_array3 = np.row_stack([binned_array2, binned_array1])
			np.savetxt(output_filename, binned_array3, fmt=b'%10.6f %10.6f %8i')
			print_cust('Saved...', quiet_val=quiet_val)
		else:
			binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(x, y, signal, noise, targetSN, plot=plot_val, quiet=quiet_val, pixelsize=PIXSIZE)
			binned_array2 = np.column_stack([x, y, binNum])
			binned_array1 = np.column_stack([x[~mask], y[~mask], [int(max(binNum)+1)]*len(x[~mask])])
			binned_array3 = np.row_stack([binned_array2, binned_array1])
			np.savetxt(output_filename, binned_array3, fmt=b'%10.6f %10.6f %8i')
			print_cust('Saved...', quiet_val=quiet_val)
######################VORONOI_BINNING##################################

######################DISPLAY_BINNING_RESULT##################################
def binning_result_display(input_filename, output_filename, par_dict, snr_type='snr', physical_axes=False, snr_vmin_val=0, snr_vmax_val=20, size_of_font=10, quiet_val=False):
	print_cust(f"Loading file: {output_filename} for plotting...", quiet_val=quiet_val)
	x_new, y_new, signal_new = np.loadtxt(output_filename).T
	signal_new_unique = np.unique(signal_new)
	x_new_unique = np.zeros_like(signal_new_unique)
	y_new_unique = np.zeros_like(signal_new_unique)
	idx_unique = np.zeros_like(signal_new_unique)
	#Plotting figure 1 only for binned data, because the loops take incredibly long time for No spatial binning
	if ('none' not in str(par_dict['binning_type'].lower())):
		for i in range(len(signal_new)):
			for j in range(len(signal_new_unique)):
				if signal_new_unique[j]==signal_new[i]:
					idx_unique[j] = int(i)
					x_new_unique[j] = x_new[i]
					y_new_unique[j] = y_new[i]
		print_cust('Plotting Figure 1...', quiet_val=quiet_val)
		plt.scatter(x_new, y_new, c = signal_new, cmap='viridis')
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		for i in range(len(x_new_unique)):
			plt.text(x_new_unique[i], y_new_unique[i], str(signal_new_unique[i]), fontsize=5, verticalalignment='top', bbox=props, alpha=1)
		plt.show()
		print_cust('Plotted Figure 1...', quiet_val=quiet_val)
	print_cust(f"Loading file: {input_filename} for plotting...", quiet_val=quiet_val)
	y, x, ra, dec, signal1, noise1, signal2, noise2 = np.loadtxt(input_filename).T
	if ('ew' in snr_type.lower()):
		signal = signal2
		noise = noise2
	else:
		signal = signal1
		noise = noise1
	snr = signal/noise
	print_cust(f"SNR/EW Limits: {np.min(snr)} to {np.max(snr)}", quiet_val=quiet_val)
	print_cust('Plotting Figure 2...', quiet_val=quiet_val)
	fig, axs = plt.subplots(2, figsize=(6,10), sharex=True, sharey=True)
	if (physical_axes):
		im1 = axs[0].scatter(dec, ra, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
	else:
		im1 = axs[0].scatter(x, y, c = snr, cmap='viridis', vmin=snr_vmin_val, vmax=snr_vmax_val)
	axs[0].set_title("SNR")
	add_colorbar_lin(im1)
	im2 = axs[1].scatter(x_new, y_new, c = signal_new, cmap='RdBu')
	axs[1].set_title("Binning")
	add_colorbar_lin(im2)
	plt.show()
	print_cust('Plotted Figure 2...', quiet_val=quiet_val)

######################DISPLAY_BINNING_RESULT##################################




####################CREATING_LINEAR_DATA_FILE_FOR_BINNING####################
def create_linear_file_for_binning(file_name_rev_linear, par_dict, header_rev, cont_rev_file2, data_rev_file2, err_rev_file2, data_ew_original_file2, err_ew_original_file2, quiet_val=False):
	print_cust('Creating linear data file for binning...', quiet_val=quiet_val)
	ra_rev, dec_rev, wave_rev = obtain_physical_axis(header_rev)
	center_position_idx = find_nearest_idx(cont_rev_file2[:, int(cont_rev_file2.shape[1]/2), int(cont_rev_file2.shape[2]/2)],  np.nanmax(cont_rev_file2[:, int(cont_rev_file2.shape[1]/2), int(cont_rev_file2.shape[2]/2)]))
	idx1 = center_position_idx - int(par_dict['window_for_choosing_snr_in_binning'])
	idx2 = center_position_idx + int(par_dict['window_for_choosing_snr_in_binning'])
	lin_x_axis = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_y_axis = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_ra = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_dec = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_signal = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_noise = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_ew = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	lin_ew_err = np.zeros([data_rev_file2.shape[1], data_rev_file2.shape[2]])
	data_ew_original_file = clean_data(data_ew_original_file2)
	err_ew_original_file = clean_data(err_ew_original_file2, type_of_data='err')
	data_ew_original_file[np.abs(data_ew_original_file)>par_dict['max_ew_width_abs']]=1e-4
	err_ew_original_file[np.abs(data_ew_original_file)>par_dict['max_ew_width_abs']]=1e-3
	bar = IncrementalBar('Countdown', max = int(data_rev_file2.shape[1]*data_rev_file2.shape[2]))
	for i in range(data_rev_file2.shape[1]):
		for j in range(data_rev_file2.shape[2]):
			bar.next()
			idx_new = arg_median(data_rev_file2[idx1:idx2,i,j])
			lin_x_axis[i,j] = i+1
			lin_y_axis[i,j] = j+1
			lin_ra[i,j] = dec_rev[i]
			lin_dec[i,j] = ra_rev[j]
			#lin_signal[i,j] = np.abs(data_rev_file2[idx_new,i,j] - sky_rev_file2[idx_new])
			#lin_noise[i,j] = np.sqrt(np.abs(err_rev_file2[idx_new,i,j] + flux_err_sky_full[idx_new]**2))
			lin_signal[i,j] = np.abs(data_rev_file2[idx_new,i,j])
			lin_noise[i,j] = np.sqrt(np.abs(err_rev_file2[idx_new,i,j]))
			lin_ew[i,j] = np.abs(data_ew_original_file[int(par_dict['lick_idx_for_halpha']), i, j])
			lin_ew_err[i,j] = np.abs(err_ew_original_file[int(par_dict['lick_idx_for_halpha']), i, j])
	lin_x_axis_flattened = lin_x_axis.flatten()
	lin_y_axis_flattened = lin_y_axis.flatten()
	lin_ra_flattened = lin_ra.flatten()
	lin_dec_flattened = lin_dec.flatten()
	lin_signal_flattened = lin_signal.flatten()
	lin_noise_flattened = lin_noise.flatten()
	lin_ew_flattened = lin_ew.flatten()
	lin_ew_err_flattened = lin_ew_err.flatten()
	lin_data_to_save = np.transpose(np.array([lin_x_axis_flattened, lin_y_axis_flattened, lin_ra_flattened, lin_dec_flattened, lin_signal_flattened, lin_noise_flattened, lin_ew_flattened, lin_ew_err_flattened]))
	np.savetxt(file_name_rev_linear, lin_data_to_save, fmt=b'%10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f')
	print_cust('Linear data file for binning created...', quiet_val=quiet_val)
####################CREATING_LINEAR_DATA_FILE_FOR_BINNING####################


