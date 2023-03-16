import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io.fits import getdata as astropygetdata
from astropy.wcs import WCS
quite_val = False
#from basic_functions import *
import basic_functions as bf


########################################HANDLING_IFU_FITS########################################
def open_ifu_fits(filename, instrument='muse', data_type='cube', quiet_val=False):
	instrument=instrument.lower()
	data_type = data_type.lower()  # data_type: cube, image, spectra
	assert os.path.exists(filename), f"File: {filename} not found..."
	if ('muse' in instrument):
		bf.print_cust('Loading MUSE-IFU data...', quiet_val=quiet_val)
		with fits.open(str(filename)) as hdul:
			header_original = hdul[1].header
			header_original_err = hdul[2].header
			data_original_file = hdul[1].data
			err_original_file = hdul[2].data
		bf.print_cust('MUSE-IFU data loaded...', quiet_val=quiet_val)
		return (header_original, header_original_err, data_original_file, err_original_file)
	elif ('sdss' in instrument):
		if ('spectra' in data_type):
			bf.print_cust('Loading SDSS spectra...', quiet_val=quiet_val)
			with fits.open(str(filename)) as hdul:
				header_original = hdul[0].header
				header_original_err = hdul[1].header
				data_original_tmp = hdul[1].data
			data_original_file = data_original_tmp['flux']
			err_original_file = np.sqrt(1./data_original_tmp['ivar']) #CHECK THIS ONCE
			bf.print_cust('SDSS spectra loaded...', quiet_val=quiet_val)
		elif ('image' in data_type):
			bf.print_cust('Loading SDSS image...', quiet_val=quiet_val)
			with fits.open(str(filename)) as hdul:
				header_original = hdul[0].header
				header_original_err = hdul[1].header
				if (hdul[0].data) is not None:
					data_original_file = hdul[0].data
				else:
					data_original_file = hdul[1].data
			err_original_file = np.zeros_like(data_original_file)
			bf.print_cust('SDSS image loaded...', quiet_val=quiet_val)
		else:
			assert False, f"Data type: {data_type} not defined for SDSS"
		return (header_original, header_original_err, data_original_file, err_original_file)
	else:
		assert False, f"Please define opening method for the instrument: {instrument} first"

def open_ifu_fits_custom_file(filename, instrument='muse', quiet_val=False):
	instrument=instrument.lower()
	if ('muse' in instrument):
		bf.print_cust('Custom loading IFU data...', quiet_val=quiet_val)
		assert os.path.exists(filename), f"File: {filename} not found..."
		data_original_file, header_original = astropygetdata(filename, 1, header=True)
		err_original_file, header_original_err = astropygetdata(filename, 2, header=True)
		cont_original_file, header_original_cont = astropygetdata(filename, 3, header=True)
		fwhm_original_file, header_original_fwhm = astropygetdata(filename, 4, header=True)
		velscale_original_file, header_original_velscale = astropygetdata(filename, 5, header=True)
		sky_original_file, header_original_sky = astropygetdata(filename, 6, header=True)
		#sky_original_file, header_original_sky = astropygetdata(filename, 6, header=True)
		bf.print_cust('Custom IFU data loaded...', quiet_val=quiet_val)
		return (header_original, header_original_err, header_original_cont, header_original_fwhm, header_original_velscale, header_original_sky, data_original_file, err_original_file, cont_original_file, fwhm_original_file, velscale_original_file, sky_original_file)
	else:
		assert False, f"Please define opening method for the instrument: {instrument} first"

def copy_ifu_fits(filename_old, file_name_new, instrument='muse', quiet_val=False):
	instrument=instrument.lower()
	if ('muse' in instrument):
		if (os.path.exists(file_name_new)):
			prompt_1=input(f"File {file_name_new} already exists. Overwrite?(y/n) : ")
			if ('y' in prompt_1.lower()):
				bf.print_cust(f"deleting...{file_name_new}", quiet_val=quiet_val)
				string_for_deleting = str('rm ') + file_name_new
				bf.print_cust(string_for_deleting, quiet_val=quiet_val)
				os.system(string_for_deleting)
				with fits.open(filename_old) as hdu_list:
					hdu_list.writeto(file_name_new)  # to write all HDUs, including the updated one, to a new file
			else:
				bf.print_cust('Overwrite skipped. Returning saved file information', quiet_val=quiet_val)
		else:
			with fits.open(filename_old) as hdu_list:
				hdu_list.writeto(file_name_new)  # to write all HDUs, including the updated one, to a new file
		#Reading back the file
		header_new, header_new_err, data_new_file, err_new_file = open_ifu_fits(file_name_new, instrument=instrument)
		return (header_new, header_new_err, data_new_file, err_new_file)
	else:
		bf.print_cust(f"Please define opening method for the instrument: {instrument} first", quiet_val=quiet_val)
########################################HANDLING_IFU_FITS########################################




########################################MODIFYING_IFU_FITS########################################
def modify_ifu_fits(filename_old, file_name_new, wave_new, data_new, **kwargs):
	quiet_val = kwargs.get('quiet_val', False)  # quiet_val
	instrument = kwargs.get('instrument', 'muse')  # Initializing Observing instrument
	instrument = instrument.lower()
	data_err_init = np.ones([data_new.shape[0], data_new.shape[1], data_new.shape[2]])
	data_err_new = kwargs.get('data_err_cube', data_err_init)  # Initializing data uncertainty
	cont_init = np.ones([data_new.shape[0], data_new.shape[1], data_new.shape[2]])
	cont_new = kwargs.get('data_cont_cube', cont_init)  # Initializing continuum
	fwhm_init = np.ones([data_new.shape[0], data_new.shape[1], data_new.shape[2]])
	fwhm_cube = kwargs.get('data_fwhm_cube', fwhm_init)  # Initializing FWHM
	velscale_init = np.ones([data_new.shape[1], data_new.shape[2]])
	velscale_2d_image = kwargs.get('data_velscale_image', velscale_init)  # Initializing Velscale
	sky_array_init = np.full((data_new.shape[0]), np.nanmin(data_new))
	sky_1d_array = kwargs.get('sky_array', sky_array_init)  # Initializing sky
	crval_cust, naxis_cust, cdn_n, cr_pix = bf.retrieve_wave_header_info(wave_new, wave_type='log')
	if ('muse' in instrument):
		if (os.path.exists(file_name_new)):
			prompt_1=input(f"File {file_name_new} already exists. Overwrite?(y/n) : ")
			if ('y' in prompt_1.lower()):
				bf.print_cust(f"deleting...{file_name_new}", quiet_val=quiet_val)
				string_for_deleting = str('rm ') + file_name_new
				bf.print_cust(string_for_deleting, quiet_val=quiet_val)
				os.system(string_for_deleting)
				header_new, header_new_err, data_new_file, err_new_file = copy_ifu_fits(filename_old, file_name_new, instrument=instrument)
				bf.print_cust('-Saving data and err', quiet_val=quiet_val)
				with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
					hdu_list[1].data = data_new
					hdu_list[2].data = data_err_new
					hdu_list[1].header['NAXIS3'] = (naxis_cust, 'length of data axis 3 ')
					hdu_list[1].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
					hdu_list[1].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
					hdu_list[1].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
					hdu_list[1].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
					hdu_list[2].header['NAXIS3'] = (naxis_cust, 'length of data axis 3 ')
					hdu_list[2].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
					hdu_list[2].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
					hdu_list[2].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
					hdu_list[2].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
				bf.print_cust('-Saved data and err', quiet_val=quiet_val)
				bf.print_cust('-Saving continuum', quiet_val=quiet_val)
				fits.append(file_name_new, cont_new, header_new)
				with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
					hdu_list[3].header['OBJECT'] = ('a2670-3 (CONT)', 'Continuum cube')
					hdu_list[3].header['EXTNAME'] = ('CONT    ', 'This extension contains data continuum values  ')
					hdu_list[3].header['HDUCLAS2'] = ('CONT    ', 'this extension contains the continuum  ')
					hdu_list[3].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
					hdu_list[3].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
					hdu_list[3].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
					hdu_list[3].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
				bf.print_cust('-Saved continuum', quiet_val=quiet_val)
				bf.print_cust('-Saving FWHM', quiet_val=quiet_val)
				fits.append(file_name_new, fwhm_cube, header_new)
				with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
					hdu_list[4].header['OBJECT'] = ('a2670-3 (fwhm)', 'FWHM cube')
					hdu_list[4].header['EXTNAME'] = ('FWHM    ', 'This extension contains data FWHM  ')
					hdu_list[4].header['HDUCLAS2'] = ('FWHM    ', 'this extension contains the FWHM  ')
					hdu_list[4].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
					hdu_list[4].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
					hdu_list[4].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
					hdu_list[4].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
				bf.print_cust('-Saved FWHM', quiet_val=quiet_val)
				bf.print_cust('-Saving Velocity Scale', quiet_val=quiet_val)
				fits.append(file_name_new, velscale_2d_image, header_new)
				with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
					hdu_list[5].header['NAXIS'] = (2, 'number of data axes  ')
					hdu_list[5].header['OBJECT'] = ('a2670-3 (velscale)', 'Velocity Scale  ')
					hdu_list[5].header['EXTNAME'] = ('VELSCALE ', 'This extension contains Velocity scale')
					hdu_list[5].header['HDUCLAS2'] = ('VELSCALE ', 'this extension contains Velocity scale')
					hdu_list[5].header['BUNIT'] = 'km/s'
					del hdu_list[5].header['ERRDATA']; del hdu_list[5].header['CTYPE3']; del hdu_list[5].header['CUNIT3']; del hdu_list[5].header['CD3_3']; del hdu_list[5].header['CRPIX3']; del hdu_list[5].header['CRVAL3']; del hdu_list[5].header['CD1_3']; del hdu_list[5].header['CD2_3']; del hdu_list[5].header['CD3_1']; del hdu_list[5].header['CD3_2']
				bf.print_cust('-Saved Velocity Scale', quiet_val=quiet_val)
				bf.print_cust('-Saving Sky Spectra', quiet_val=quiet_val)
				fits.append(file_name_new, sky_1d_array, header_new)
				with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
					hdu_list[6].header['NAXIS'] = (1, 'number of data axes  ')
					hdu_list[6].header['NAXIS1'] = (naxis_cust, 'length of data axis 1 ')
					hdu_list[6].header['CRPIX1'] = header_new['CRPIX3']
					hdu_list[6].header['CD1_1'] = header_new['CD3_3']
					hdu_list[6].header['CRVAL1'] = header_new['CRVAL3']
					hdu_list[6].header['CUNIT1'] = header_new['CUNIT3']
					hdu_list[6].header['CTYPE1'] = header_new['CTYPE3']
					hdu_list[6].header['OBJECT'] = ('a2670-3 (SKY)', 'Sky Spectra  ')
					hdu_list[6].header['EXTNAME'] = ('SKY    ', 'This extension contains the Sky spectra  ')
					hdu_list[6].header['HDUCLAS2'] = ('SKY    ', 'this extension contains the Sky spectra  ')
					hdu_list[6].header['CRPIX1'] = (cr_pix, 'Pixel coordinate of reference point  ')
					hdu_list[6].header['CD1_1'] = (cdn_n, 'Coordinate transformation matrix element')
					hdu_list[6].header['CRVAL1'] = (crval_cust, 'Value at reference point pixel')
					hdu_list[6].header['CTYPE1'] = ('AWAV    ', 'log Axis data type')
					del hdu_list[6].header['CRPIX2']; del hdu_list[6].header['CRPIX3']; del hdu_list[6].header['CD1_2']; del hdu_list[6].header['CD1_3']; del hdu_list[6].header['CD2_1']; del hdu_list[6].header['CD2_2']; del hdu_list[6].header['CD2_3']; del hdu_list[6].header['CD3_1']; del hdu_list[6].header['CD3_2']; del hdu_list[6].header['CD3_3']; del hdu_list[6].header['CUNIT2']; del hdu_list[6].header['CUNIT3']; del hdu_list[6].header['CTYPE2']; del hdu_list[6].header['CTYPE3']; del hdu_list[6].header['CSYER1']; del hdu_list[6].header['CSYER2']; del hdu_list[6].header['CRVAL2']; del hdu_list[6].header['CRVAL3']
				bf.print_cust('-Saved Sky Spectra', quiet_val=quiet_val)
				bf.print_cust(f"New file: {file_name_new} saved", quiet_val=quiet_val)
			else:
				bf.print_cust('Overwrite skipped. Returning saved file information', quiet_val=quiet_val)
		else:
			header_new, header_new_err, data_new_file, err_new_file = copy_ifu_fits(filename_old, file_name_new, instrument=instrument)
			bf.print_cust('-Saving data and err', quiet_val=quiet_val)
			with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
				hdu_list[1].data = data_new
				hdu_list[2].data = data_err_new
				hdu_list[1].header['NAXIS3'] = (naxis_cust, 'length of data axis 3 ')
				hdu_list[1].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
				hdu_list[1].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
				hdu_list[1].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
				hdu_list[1].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
				hdu_list[2].header['NAXIS3'] = (naxis_cust, 'length of data axis 3 ')
				hdu_list[2].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
				hdu_list[2].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
				hdu_list[2].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
				hdu_list[2].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
			bf.print_cust('-Saved data and err', quiet_val=quiet_val)
			bf.print_cust('-Saving continuum', quiet_val=quiet_val)
			fits.append(file_name_new, cont_new, header_new)
			with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
				hdu_list[3].header['OBJECT'] = ('a2670-3 (CONT)', 'Continuum cube')
				hdu_list[3].header['EXTNAME'] = ('CONT    ', 'This extension contains data continuum values  ')
				hdu_list[3].header['HDUCLAS2'] = ('CONT    ', 'this extension contains the continuum  ')
				hdu_list[3].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
				hdu_list[3].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
				hdu_list[3].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
				hdu_list[3].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
			bf.print_cust('-Saved continuum', quiet_val=quiet_val)
			bf.print_cust('-Saving FWHM', quiet_val=quiet_val)
			fits.append(file_name_new, fwhm_cube, header_new)
			with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
				hdu_list[4].header['OBJECT'] = ('a2670-3 (fwhm)', 'FWHM cube')
				hdu_list[4].header['EXTNAME'] = ('FWHM    ', 'This extension contains data FWHM  ')
				hdu_list[4].header['HDUCLAS2'] = ('FWHM    ', 'this extension contains the FWHM  ')
				hdu_list[4].header['CRPIX3'] = (cr_pix, 'Pixel coordinate of reference point  ')
				hdu_list[4].header['CD3_3'] = (cdn_n, 'Coordinate transformation matrix element')
				hdu_list[4].header['CRVAL3'] = (crval_cust, 'Value at reference point pixel')
				hdu_list[4].header['CTYPE3'] = ('AWAV    ', 'log Axis data type')
			bf.print_cust('-Saved FWHM', quiet_val=quiet_val)
			bf.print_cust('-Saving Velocity Scale', quiet_val=quiet_val)
			fits.append(file_name_new, velscale_2d_image, header_new)
			with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
				hdu_list[5].header['NAXIS'] = (2, 'number of data axes  ')
				hdu_list[5].header['OBJECT'] = ('a2670-3 (velscale)', 'Velocity Scale  ')
				hdu_list[5].header['EXTNAME'] = ('VELSCALE ', 'This extension contains Velocity scale')
				hdu_list[5].header['HDUCLAS2'] = ('VELSCALE ', 'this extension contains Velocity scale')
				hdu_list[5].header['BUNIT'] = 'km/s'
				del hdu_list[5].header['ERRDATA']; del hdu_list[5].header['CTYPE3']; del hdu_list[5].header['CUNIT3']; del hdu_list[5].header['CD3_3']; del hdu_list[5].header['CRPIX3']; del hdu_list[5].header['CRVAL3']; del hdu_list[5].header['CD1_3']; del hdu_list[5].header['CD2_3']; del hdu_list[5].header['CD3_1']; del hdu_list[5].header['CD3_2']
			bf.print_cust('-Saved Velocity Scale', quiet_val=quiet_val)
			bf.print_cust('-Saving Sky Spectra', quiet_val=quiet_val)
			fits.append(file_name_new, sky_1d_array, header_new)
			with fits.open(file_name_new, mode='update', output_verify='fix') as hdu_list:
				hdu_list[6].header['NAXIS'] = (1, 'number of data axes  ')
				hdu_list[6].header['NAXIS1'] = (naxis_cust, 'length of data axis 1 ')
				hdu_list[6].header['CRPIX1'] = header_new['CRPIX3']
				hdu_list[6].header['CD1_1'] = header_new['CD3_3']
				hdu_list[6].header['CRVAL1'] = header_new['CRVAL3']
				hdu_list[6].header['CUNIT1'] = header_new['CUNIT3']
				hdu_list[6].header['CTYPE1'] = header_new['CTYPE3']
				hdu_list[6].header['OBJECT'] = ('a2670-3 (SKY)', 'Sky Spectra  ')
				hdu_list[6].header['EXTNAME'] = ('SKY    ', 'This extension contains the Sky spectra  ')
				hdu_list[6].header['HDUCLAS2'] = ('SKY    ', 'this extension contains the Sky spectra  ')
				hdu_list[6].header['CRPIX1'] = (cr_pix, 'Pixel coordinate of reference point  ')
				hdu_list[6].header['CD1_1'] = (cdn_n, 'Coordinate transformation matrix element')
				hdu_list[6].header['CRVAL1'] = (crval_cust, 'Value at reference point pixel')
				hdu_list[6].header['CTYPE1'] = ('AWAV    ', 'log Axis data type')
				del hdu_list[6].header['CRPIX2']; del hdu_list[6].header['CRPIX3']; del hdu_list[6].header['CD1_2']; del hdu_list[6].header['CD1_3']; del hdu_list[6].header['CD2_1']; del hdu_list[6].header['CD2_2']; del hdu_list[6].header['CD2_3']; del hdu_list[6].header['CD3_1']; del hdu_list[6].header['CD3_2']; del hdu_list[6].header['CD3_3']; del hdu_list[6].header['CUNIT2']; del hdu_list[6].header['CUNIT3']; del hdu_list[6].header['CTYPE2']; del hdu_list[6].header['CTYPE3']; del hdu_list[6].header['CSYER1']; del hdu_list[6].header['CSYER2']; del hdu_list[6].header['CRVAL2']; del hdu_list[6].header['CRVAL3']
			bf.print_cust('-Saved Sky Spectra', quiet_val=quiet_val)
			bf.print_cust(f"New file: {file_name_new} saved", quiet_val=quiet_val)
	else:
		bf.print_cust(f"Please define opening method for the instrument: {instrument} first", quiet_val=quiet_val)
########################################MODIFYING_IFU_FITS########################################



########################################OBTAIN_PHYSICAL_AXIS_FROM_HEADER########################################
def obtain_physical_axis(header_original, **kwargs):
	obtain_type = kwargs.get('obtain_type','wcs')  # Method of fit
	instrument = kwargs.get('instrument','MUSE')  # instrument: MUSE, SDSS
	data_type = kwargs.get('data_type','cube')  # data_type: cube, image, spectra
	quiet_val = kwargs.get('quiet_val', False)  # quiet_val
	if (('sdss' in instrument) and ('image' in data_type)):
		bf.print_cust(r'Assuming SDSS-IMAGE coordinate system(RA(deg), Dec(deg))', quiet_val=quiet_val)
		sdss_wcs = WCS(header_original).celestial
		axis_physical_x = np.zeros([sdss_wcs._naxis[1]])
		axis_physical_y = np.zeros([sdss_wcs._naxis[0]])
		for i in range(len(axis_physical_x)):
			axis_physical_x[i] = sdss_wcs.array_index_to_world_values([i],[0])[1][0]
		for j in range(len(axis_physical_y)):
			axis_physical_y[j] = sdss_wcs.array_index_to_world_values([0],[j])[0][0]
		axis_physical_z = [0]
		return (axis_physical_y, axis_physical_x, axis_physical_z)

	if (('sdss' in instrument) and ('spectra' in data_type)):
		bf.print_cust(r'Assuming SDSS-spectral coordinate system(Wavelength(Ang))', quiet_val=quiet_val)
		data_original = kwargs.get('data_original', np.arange(1000))  # secondary_header
		c0 = header_original['coeff0']
		c1 = header_original['coeff1']
		#npix = header_original['naxis1']
		npix = len(data_original)
		axis_physical_z = 10.**(c0 + c1 * np.arange(npix))
		axis_physical_x = [header_original['PLUG_DEC']]
		axis_physical_y = [header_original['PLUG_RA']]
		return (axis_physical_y, axis_physical_x, axis_physical_z)

	if (('MUSE' in instrument) and ('cube' in data_type)):
		bf.print_cust(r'Assuming VLT-MUSE coordinate system(RA(deg), Dec(deg), Wavelength(Ang))', quiet_val=quiet_val)
		if (obtain_type=='wcs'):
			bf.print_cust('Retrieving physical axis from astropy.wcs.WCS', quiet_val=quiet_val)
			assert WCS(header_original).has_celestial, f"Header does not have celestial coordinate information..."
			assert WCS(header_original).has_spectral, f"Header does not have a spectral coordinate information..."
			muse_wcs = WCS(header_original).celestial
			axis_physical_x = np.zeros([muse_wcs._naxis[1]])
			axis_physical_y = np.zeros([muse_wcs._naxis[0]])
			for i in range(len(axis_physical_x)):
				axis_physical_x[i] = muse_wcs.array_index_to_world_values([i],[0])[1][0]
			for j in range(len(axis_physical_y)):
				axis_physical_y[j] = muse_wcs.array_index_to_world_values([0],[j])[0][0]
			muse_spectral = WCS(header_original).spectral
			axis_physical_z = np.zeros([muse_spectral._naxis[0]])
			for z in range(len(axis_physical_z)):
				axis_physical_z[z] = muse_spectral.array_index_to_world_values([z])[0]*1e10
			if ('log' in header_original.comments['CTYPE3']):
				axis_physical_z = 10**(axis_physical_z)
			return (axis_physical_y, axis_physical_x, axis_physical_z)
		elif (obtain_type=='manual'):
			bf.print_cust('Retrieving physical axis from manual code by Axe(05/10/2022)', quiet_val=quiet_val)
			assert header_original['NAXIS']<=3, f"Header does not have 3 dimensions. Found only {header_original['NAXIS']}"
			start_value_x_axis_physical = (float(header_original['CRVAL1'])) - ((float(header_original['CRPIX1']-1.)) * (float(header_original['CD1_1'])))
			start_value_y_axis_physical = (float(header_original['CRVAL2'])) - ((float(header_original['CRPIX2']-1.)) * (float(header_original['CD2_2'])))
			start_value_z_axis_physical = (float(header_original['CRVAL3'])) - ((float(header_original['CRPIX3']-1.)) * (float(header_original['CD3_3'])))
			x_axis = np.arange(header_original['NAXIS1'])
			y_axis = np.arange(header_original['NAXIS2'])
			z_axis = np.arange(header_original['NAXIS3'])
			x_axis_real = start_value_x_axis_physical
			y_axis_real = start_value_y_axis_physical
			z_axis_real = start_value_z_axis_physical
			for i in range(len(x_axis)-1):
				x_axis_real = np.append(x_axis_real, (start_value_x_axis_physical + ((i+1.)*(float(header_original['CD1_1'])))))
			for j in range(len(y_axis)-1):
				y_axis_real = np.append(y_axis_real, (start_value_y_axis_physical + ((j+1.)*(float(header_original['CD2_2'])))))
			for k in range(len(z_axis)-1):
				z_axis_real = np.append(z_axis_real, (start_value_z_axis_physical + ((k+1.)*(float(header_original['CD3_3'])))))
			if ('log' in header_original.comments['CTYPE3']):
				z_axis_real = 10**(z_axis_real)
			return(x_axis_real, y_axis_real, z_axis_real)
		else:
			bf.print_cust('Type for obtaining physical axes unidentified. Quitting...', quiet_val=quiet_val)
			quit()
########################################OBTAIN_PHYSICAL_AXIS_FROM_HEADER########################################


