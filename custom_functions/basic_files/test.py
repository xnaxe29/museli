import numpy as np
import pandas as pd

file_path = 'default_lick_indexes.dat'
data = pd.read_csv(file_path, sep="\t", comment='#', encoding=None)
#print(data['wave11'])
#print(data['wave21'])
#print(data['wave31'])
#print(data['sign'])
#print(data['species'])
#print(data['species'][65])
#print(data['species'][69])
#print(data['species'][80])
#print(data['species'][83])
#print(data['species'][88])
#quit()

file_path = 'default_emission_lines.dat'
data = pd.read_csv(file_path, sep="\t", comment='#', encoding=None)
idx_new = [15, 23, 24, 25, 26, 27]
#print (data['Index'][idx_new])
#print (dir(data['Name'][idx_new]))
#print (type(data['Name'][idx_new].to_numpy()))
#print (data['Lambda'][idx_new].to_numpy().astype(np.int32))
A = (data['Name'][idx_new].to_numpy().astype(np.str_))
B = (data['Lambda'][idx_new].to_numpy().astype(np.int32))
C = np.vstack([A,B])
name_new = np.chararray([len(A)], itemsize=20)
for i in range(len(A)):
	name_new[i] = str(A[i]) + str(B[i])
name_new = name_new.decode("utf-8")
#print (data['Lambda'][idx_new].to_numpy().astype(np.float32))
#print (data['comments_on_balmer'][idx_new].to_numpy().astype(np.bool_))
#print (data['comments_on_tied'][idx_new].to_numpy().astype(np.str_))
#print (data['factor_for_tied'][idx_new].to_numpy().astype(np.float32))
#print (data['factor_fixed_for_tying'][idx_new].to_numpy().astype(np.bool_))

redshift = 0.00
wave_redshifted = (data['Lambda'][idx_new].to_numpy().astype(np.float32))*(1.+redshift)
print (wave_redshifted)
wave_min = 6000.
wave_max = 10000.
mask = (wave_min<wave_redshifted) & (wave_redshifted<wave_max)
print (wave_redshifted[mask])
print (name_new[mask])


quit()

#print(data['comments_on_balmer'])
#print(data['comments_on_tied'])
#print(data['factor_for_tied'])
#print(data['factor_fixed_for_tying'])


file_path = 'default_sky.dat'
data = pd.read_csv(file_path, sep="\t", comment='#', encoding=None)
print(data['wave'])
print(data['sky_flux'])
print(data['sky_flux_err'])
