import numpy as np
from scipy import interpolate
from os import path
#from basic_functions import *
import basic_functions as bf
quite_val = False


##################################EXTINCTION_FITTING##################################
##################################REDENNING_FUNCTION##################################
#Never use SMC Wing for anything - Tip from JK, hence 'SMC Wing' removed
#raj = np.array(['Galactic', 'SMC_Bar', 'LMC_Supershell', 'LMC_Average'])

def func_6(xdata, ydata, *t):
	E_bv_dla, flux_redu, E_bv_qso, del_beta=t
	extinction_class_string = 'Galactic'
	x0_dla, gamma_dla, c1_dla, c2_dla, c3_dla, c4_dla, c5_dla, O1_dla, O2_dla, O3_dla, R_v_dla, k_IR_dla = extinction(extinction_class_string)
	x0_qso, gamma_qso, c1_qso, c2_qso, c3_qso, c4_qso, c5_qso, O1_qso, O2_qso, O3_qso, R_v_qso, k_IR_qso = extinction('SMC_Bar')
	#model = ((reddening_func2(ydata, wave, E_bv_dla, R_v_dla, spline_part_new_fit(xdata, extinction_class_string, c3_dla), 0., R_v_qso, spline_part_new_fit_qso(xdata), 0.))/flux_redu)
	model = ((reddening_func2(ydata, xdata, E_bv_dla, R_v_dla, spline_part_new_fit(xdata, extinction_class_string, c3_dla), E_bv_qso, R_v_qso, spline_part_new_fit_qso(xdata), del_beta))/flux_redu)
	return (model)

##################################REDENNING_FUNCTION##################################

#######################################DEFINING_FIXED_PARAMETERS#######################################

def extinction(extinction_class):
	if (extinction_class=='Galactic'):
		x0 = 4.592; gamma = 0.922; c1 = -0.175;	c2 = 0.807; c3 = 2.991; c4 = 0.319; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 3.001; k_IR = 1.057; x0_err = 0.00; gamma_err = 0.00; c1_err = 0.00; c2_err = 0.00; c3_err = 0.00; c4_err = 0.00
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='SMC_Bar'):
		x0 = 4.600; gamma = 1.000; c1 = -4.959; c2 = 2.264; c3 = 0.389; c4 = 0.461; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.74; k_IR = 1.057; x0_err = 0.00; gamma_err = 0.00; c1_err = 0.197; c2_err = 0.040; c3_err = 0.110; c4_err = 0.079
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='SMC_Wing'):
		x0 = 4.703; gamma = 1.212; c1 = -0.856; c2 = 1.038; c3 = 3.215;	c4 = 0.107; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.05; k_IR = 1.057; x0_err = 0.018; gamma_err = 0.019; c1_err = 0.246; c2_err = 0.074; c3_err = 0.439; c4_err = 0.038
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='LMC_Supershell'):
		x0 = 4.558; gamma = 0.945; c1 = -1.475; c2 = 1.132; c3 = 1.463; c4 = 0.294; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 2.76; k_IR = 1.057; x0_err = 0.021; gamma_err = 0.026; c1_err = 0.152; c2_err = 0.029; c3_err = 0.121; c4_err = 0.057
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00

	elif (extinction_class=='LMC_Average'):
		x0 = 4.579; gamma = 0.934; c1 = -0.890; c2 = 0.998; c3 = 2.719; c4 = 0.400; c5 = 6.097; O1 = 2.055; O2 = 1.322; O3 = 0.000
		R_v = 3.41; k_IR = 1.057; x0_err = 0.007; gamma_err = 0.016; c1_err = 0.142; c2_err = 0.027; c3_err = 0.137; c4_err = 0.036
		c5_err = 0.00; O1_err = 0.00; O2_err = 0.00; O3_err = 0.00; R_v_err = 0.00; k_IR_err = 0.00


	return (x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR)

#######################################DEFINING_FIXED_PARAMETERS#######################################

####################################REDDENING_FUNCTION####################################
def reddening_func2(F_rest_lambda, wave, E_bv_dla, R_v_dla, initial_result_dla, E_bv_qso, R_v_qso, initial_result_qso, del_beta):
	constant = np.zeros([len(F_rest_lambda)])
	F_lambda = np.zeros([len(F_rest_lambda)])
	for i in range(len(F_rest_lambda)):
		constant[i] = -0.4 * ((E_bv_dla * (initial_result_dla[i] + R_v_dla)) + (E_bv_qso * (initial_result_qso[i] + R_v_qso)))
		F_lambda[i] = ((F_rest_lambda[i]*((wave[i]/5510.)**(del_beta))) * (10.**(constant[i])))
	return (F_lambda)
	
def k_lambda_V(x, x0, gamma, c1, c2, c3, c4, c5):
	D_func = np.zeros([len(x)])
	result = np.zeros([len(x)])
	for i in range(len(x)):
		D_func[i] = x[i]**2. / ( ( (x[i]**2.) - (x0**2.) )**2. + (x[i] * gamma)**2. )
		if (x[i] <= c5):
			result[i] = c1 + (c2*x[i]) + c3*D_func[i]
		else:
			result[i] = c1 + (c2*x[i]) + c3*D_func[i] + c4*((x[i]-c5)**2)
	return (result)
	
#GLOBAL PARAMETER DEFINITIONS
U1_pos = 0.27; U2_pos = 0.26; O1_pos = 0.33; O2_pos = 0.4; O3_pos = 0.553; O4_pos = 0.7; I1_pos = 0; I2_pos = 1/0.25; I3_pos = 1/0.50; I4_pos = 1/0.75; I5_pos = 1/1.00
####################################REDDENING_FUNCTION####################################

#########################################REDDENING_BACKBONE_DLA#########################################
def spline_part_new_fit(wave, extinction_class_string, c3_fit_new):
	z = 0.0
	x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR = extinction(extinction_class_string)
	c3=c3_fit_new
	I1_pos_new = 0.00000000001
	I2_pos_new = 1/I2_pos
	I3_pos_new = 1/I3_pos
	I4_pos_new = 1/I4_pos
	I5_pos_new = 1/I5_pos
	O1_val = O1
	O2_val = O2
	O3_val = O3
	def I_n(k_IR, pos, R_v):
		return (((k_IR * (pos**(-1.84))) - R_v))
	I1_val = I_n(k_IR, I1_pos_new, R_v)
	I2_val = I_n(k_IR, I2_pos_new, R_v)
	I3_val = I_n(k_IR, I3_pos_new, R_v)
	I4_val = I_n(k_IR, I4_pos_new, R_v)
	I5_val = I_n(k_IR, I5_pos_new, R_v)
	x_new = np.array([1/O1_pos, 1/O2_pos, 1/O3_pos, 1/I4_pos, 1/I3_pos, 1/I2_pos, 0])
	y_new = np.array([O1_val, O2_val, O3_val, I4_val, I3_val, I2_val, I1_val])
	func = interpolate.interp1d(x_new, y_new)
	wave_new = (wave*1e-4)/(1+z)
	x = np.sort(1/(wave_new))
	x1 = np.array([])
	x2 = np.array([])
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			x1 = np.append(x1, x[i])
		else:
			x2 = np.append(x2, x[i])
	y1 = func(x1)
	position = np.searchsorted(x, (1/O1_pos))
	y = np.zeros([len(x)])
	y2 = k_lambda_V(x2, x0, gamma, c1, c2, c3, c4, c5)
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			y[i] = y1[i]
		else:
			y[i] = y2[i-position]
	y_k_lambda_V = y
	y_k_lambda_V = bf.flip(y_k_lambda_V, 0)
	return (y_k_lambda_V)
#########################################REDDENING_BACKBONE_DLA#########################################

#########################################REDDENING_BACKBONE_QSO#########################################

def spline_part_new_fit_qso(wave):
	z = 0.0
	extinction_class_string = 'SMC_Bar'
	x0, gamma, c1, c2, c3, c4, c5, O1, O2, O3, R_v, k_IR = extinction(extinction_class_string)
	#c3=c3_fit_new
	I1_pos_new = 0.00000000001
	I2_pos_new = 1/I2_pos
	I3_pos_new = 1/I3_pos
	I4_pos_new = 1/I4_pos
	I5_pos_new = 1/I5_pos
	O1_val = O1
	O2_val = O2
	O3_val = O3
	def I_n(k_IR, pos, R_v):
		return (((k_IR * (pos**(-1.84))) - R_v))
	I1_val = I_n(k_IR, I1_pos_new, R_v)
	I2_val = I_n(k_IR, I2_pos_new, R_v)
	I3_val = I_n(k_IR, I3_pos_new, R_v)
	I4_val = I_n(k_IR, I4_pos_new, R_v)
	I5_val = I_n(k_IR, I5_pos_new, R_v)
	x_new = np.array([1/O1_pos, 1/O2_pos, 1/O3_pos, 1/I4_pos, 1/I3_pos, 1/I2_pos, 0])
	y_new = np.array([O1_val, O2_val, O3_val, I4_val, I3_val, I2_val, I1_val])
	func = interpolate.interp1d(x_new, y_new)
	wave_new = (wave*1e-4)/(1+z)
	x = np.sort(1/(wave_new))
	x1 = np.array([])
	x2 = np.array([])
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			x1 = np.append(x1, x[i])
		else:
			x2 = np.append(x2, x[i])
	y1 = func(x1)
	position = np.searchsorted(x, (1/O1_pos))
	y = np.zeros([len(x)])
	y2 = k_lambda_V(x2, x0, gamma, c1, c2, c3, c4, c5)
	for i in range(len(x)):
		if (x[i] < (1/O1_pos)):
			y[i] = y1[i]
		else:
			y[i] = y2[i-position]
	y_k_lambda_V = y
	y_k_lambda_V = bf.flip(y_k_lambda_V, 0)
	return (y_k_lambda_V)

#########################################REDDENING_BACKBONE_QSO#########################################
##################################EXTINCTION_FITTING##################################



