       ___           ___           ___           ___           ___             
      /\__\         /\__\         /\  \         /\  \         /\__\      ___   
     /::|  |       /:/  /        /::\  \       /::\  \       /:/  /     /\  \  
    /:|:|  |      /:/  /        /:/\ \  \     /:/\:\  \     /:/  /      \:\  \ 
   /:/|:|__|__   /:/  /  ___   _\:\~\ \  \   /::\~\:\  \   /:/  /       /::\__\
  /:/ |::::\__\ /:/__/  /\__\ /\ \:\ \ \__\ /:/\:\ \:\__\ /:/__/     __/:/\/__/
  \/__/~~/:/  / \:\  \ /:/  / \:\ \:\ \/__/ \:\~\:\ \/__/ \:\  \    /\/:/  /   
        /:/  /   \:\  /:/  /   \:\ \:\__\    \:\ \:\__\    \:\  \   \::/__/    
       /:/  /     \:\/:/  /     \:\/:/  /     \:\ \/__/     \:\  \   \:\__\    
      /:/  /       \::/  /       \::/  /       \:\__\        \:\__\   \/__/    
      \/__/         \/__/         \/__/         \/__/         \/__/            


DOCUMENTATION

DATE: 29/11/2022
DEVELOPERS: Adarsh Ranjan, Garreth Martin, Christophe Saulder, and Satadru Bag


INTRODUCTION:

MUSELI is a graphical user interface (GUI) tool python to perform stellar population synthesis and gas dynamical analysis in MUSE IFU spectroscopic data. It is also very easy to extend this program to analyze any other kind of IFU spectroscopic data. The code is completely written in python3 and includes the usage of multiple tools developed by Michele Capellari such as ppxf for stellar absorption fitting, Vorbin for voronoi spatial binning. The code also uses continuum estimation tool inspired from the method by Martin+2021. In future, the code will include 'pafit' routine to fit the position angles developed by  


REQUIREMENTS GUIDE:

python3: 3.8.5
The following python3 modules are required:
pip=22.2.2
requests=2.24.0
pyfiglet=0.8.post1
matplotlib=3.3.2
h5py=3.7.0
progress=1.6
astropy=4.2.1
pandas=1.1.3
scipy=1.7.3

An easy way to do this is through conda using the help of the following commands:

1. Remove conda environment with the specific name, "museli_env", if any - 
->conda env remove -n museli_env

2. Create conda environment - 
->conda env create --file museli_env.yml

3. Check whether the conda environment, "museli_env" is created
->conda env list

4. Activate the environment, "museli_env"
->conda activate museli_env

4.1 Make Sure pip is installed 
-> python3 -m ensurepip

5. Install additional pip requirements using file - "museli_pip_requirements.txt" on the conda environment, "museli_env" - 
->python3 -m pip install -r museli_pip_requirements.txt

6. While running from any directory, activate the conda environment, "museli_env" - 
->conda activate museli_env


INSTALLATION GUIDE:

Since this code does not have an executable, it can be easily run by setting an alias in the bash or csh file in Mac or Linux. In bash, this would be - 
->alias museli="python3 [PATH TO MUSELI DIRECTORY]/museli/main_exec.py"

Then type, 
->museli
into the command line from any directory to run the code

Code requires a parameter file (default name: "parameter_file.dat"). Although a customized parameter file, (e.g. "custom_parameter_file.dat") can be included as a argument for running museli as - 
->museli custom_parameter_file.dat

There is only one mandatory key to be included in the parameter file, "muse_data_filename" as - 
->muse_data_filename: MUSE_DATA_FILE.fits

Please see section~{PARAMETER FILE KEYS} for further information on different type of keys to be included.


PARAMETER FILE KEYS:
This section contains keys and the corresponding acceptable values for the keys following a short description about the key.  

1. muse_data_filename: MUSE_DATA_FILE.fits
The muse data file that is to be analyzed.

2. sky_spectra_file: default_sky.dat
Corresponding sky file in the binary format.
Please note to include the header in the file as - 
"wave	sky_flux	sky_flux_err"
Observed wavelength (key: wave)
Observed sky flux (key: sky_flux)
Corresponding uncertainty on the observed sky flux (key: sky_flux_err)

3. lick_index_file: default_lick_indexes.dat
The lick index file that is used for galaxy absorption line estimation.
Please note to include the header in the file as - 
"wave11	wave12	wave21	wave22	wave31	wave32	sign	species	reference"
Rest wavelength11 (key: wave11)
Rest wavelength12 (key: wave12)
Rest wavelength21 (key: wave21)
Rest wavelength22 (key: wave22)
Rest wavelength31 (key: wave31)
Rest wavelength32 (key: wave32)
Sign (key: sign) (A or M)
Ion Species name and wavelength (key: species)
Reference for the line (key: reference)

4. lick_index_species_list: List of indexes to be used from the lick_index_file for equivalent width estimation. Please include indexes of all important stellar absorption lines.
lick_index_species_list: [65, 69, 80, 83, 88]

5. lick_idx_for_halpha: Lick index specifying the position of Halpha in the lick_index_species_list. Since Halpha is the most important line both in emission and absorption, it will be important to extracting basic information for halpha for binning and other analysis. If Halpha is not in observing range, please include the index of a species with very prominently strong absorption or emission feature. Please note that this index is relative to lick_index_species_list, not the lick_index_file   
lick_idx_for_halpha: 4


6. emission_file: default_emission_lines.dat
The emission index file that is used for gas emission line estimation.
Please note to include the header in the file as - 
"Index	Name	Lambda	action	l_kind	A_i	V_g_i	sig_g_i	fit_kind	AoN	Comment	comments_on_balmer	comments_on_tied	factor_for_tied	factor_fixed_for_tying"
Index (key: Index)
Name of the species (key: Name)
Rest wavelength (key: Lambda)
Action (key: action) f/i/m
l_kind (key: l_kind) l/d#
line_ratio (key: A_i) 1/relative_ratio
V_g_i (key: V_g_i) in km/s 
sig_g_i (key: V_g_i) in km/s 
fit_kind (key: fit_kind) f/h/t
AoN (key: AoN) 
Comment on the line (key: Comment) 
comments_on_balmer (key: comments_on_balmer) See section~EMISSION_FITTING_DETAILS
comments_on_tied (key: comments_on_tied) See section~EMISSION_FITTING_DETAILS
factor_for_tied (key: factor_for_tied) See section~EMISSION_FITTING_DETAILS
factor_fixed_for_tying (key: factor_fixed_for_tying) See section~EMISSION_FITTING_DETAILS

7. emission_index_species_list: List of indexes to be used from the emission_file for emission line fitting. Please include indexes of all important gas emission lines.
(key: emission_index_species_list) 
(Default: [15, 23, 24, 25, 26, 27])

8. quiet: Supresses printing all comments. 
Allowed values: {True or False}; Default: False

9. spectral_smoothing: Smoothing in the dispersive/spectral dimension
Allowed values: {above 1}; Default: 1(no smoothing)

10. redshift_val: Rough estimate of redshift of the galaxy
Allowed values: {above 0}; Default: 0.0(no redshift)

11. binning_type: Type of binning
Allowed values: {'None', 'square', 'advanced', 'voronoi'}; Default: None. See section~BINNING_DETAILS for more information.

12. binning_quant: Quantification of the binning value. See section~BINNING_DETAILS for more information.
Allowed values: {above 0}; Default: 10

13. region_of_binning_interest: Region on which the binning needs to be applied. See section~BINNING_DETAILS for more information.
Allowed values: {[start_x, start_y, end_x, end_y] }; Default: [0, 0, -1, -1] 

14. region_of_binning_interest_type: Type of quantification for the binning. Logical (numbers) or physical(ra/dec) indexes. See section~BINNING_DETAILS for more information.
Allowed values: {['logical' or 'physical'] }; Default: [logical] 

15. voronoi_snr_type: Is binning to be done based on SNR or Equivalent Width?. See section~BINNING_DETAILS for more information.
Allowed values: {['snr' or 'ew'] }; Default: [snr] 

16. binning_radius_list: For advanced binning, the circular regions in which the arcs needs to be cutout in logical units. See section~BINNING_DETAILS for more information.
Allowed values: {arc sector radius list}; Default: [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 500]

17. center_coordinate_for_advanced_binning: Center coordinates for advanced binning in logical units. See section~BINNING_DETAILS for more information.
Allowed values: {center logical coordinates: x, y}; Default: [0, 0]

18. ellipse_angle_for_advanced_binning: Angle in which the ellipse will be tilted for advanced binning in degrees. See section~BINNING_DETAILS for more information.
Allowed values: {0 to 360}; Default: 0

19. window_for_choosing_snr_in_binning: The code will automatically decide a center spectral position in logical units. This indicates the number 'n' for window such that data and noise will be picked up from median(data[center-n: center+n]) and median(noise[center-n: center+n])
Allowed values: {1 to logical length of the wavelength array}; Default: 100

20. execution_type: What would one like the code to do? One can do complete fitting ('fit') or stop to check binning ('snr_map'). 
Allowed values: {'snr_map', 'fit'}; Default: snr_map

21. execution_fit_type: Type of fitting to be done: stellar absorption ('absorption'), gas emission ('emission') or automatic estimation using equivalent width sign ('auto')
Allowed values: {'auto', 'absorption', 'emission'}; Default: auto

21.2 decider_ew: If execution_fit_type is 'auto', then decider_ew decides whether a spaxel has to be fitted for emission (<decider_ew) or absorption (>decider_ew)
Allowed values: {float}; Default: 0.0

22. emission_fit_type: Type of emission fitting to be done: using ppxf ('ppxf') or in-house code ('custom')
Allowed values: {'ppxf', 'custom'}; Default: custom

23. ppxf_stars_comp: Number of stellar components to be fitted using ppxf
Allowed values: {0 or above}; Default: 1

23. ppxf_gas_comp: Number of gas components to be fitted using ppxf
Allowed values: {0 or above}; Default: 1

24. ppxf_emission_tie_balmer: Tying Balmer lines for emission fit in ppxf
Allowed values: {True or False}; Default: True

25. ppxf_emission_limit_doublets: Limit doublets for emission fit in ppxf
Allowed values: {True or False}; Default: True

26. st_age_list_start: To estimate morphology from stellar absorption fit, the code needs a list of stellar ages (in Gyr) in which the population can be segregated. Then the code will output the stellar morphology for each individual bins, within st_age_list_start[index] to st_age_list_end[index].  
Allowed values: {start list of all stellar ages in Gyr}; Default: [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]

27. st_age_list_end: To estimate morphology from stellar absorption fit, the code needs a list of stellar ages (in Gyr) in which the population can be segregated. Then the code will output the stellar morphology for each individual bins, within st_age_list_start[index] to st_age_list_end[index].  
Allowed values: {end list of all stellar ages in Gyr}; Default: [0.1, 0.5, 1.0, 5.0, 10.0, 20.0]

28. number_of_narrow_components: Number of narrow emission components (<100 km/s) to be fitted using custom emission line fit.
Allowed values: {0 or above}; Default: 1

29. number_of_wide_components: Number of wide emission components (>100 km/s) to be fitted using custom emission line fit.
Allowed values: {0 or above}; Default: 0

30. stopping_number_for_continuum: While estimating order of the continuum, a number needs to be fixed such that the chebyshev polynomial outputs a decent fit.
Allowed values: {1 or above}; Default: 5e6

31. minimal_tying_factor_variance: While the ratio for two lines of the same species is tied and fixed, this value would represent how much the polynomial can move. This is required as scipy gives garbage uncertainty to a fit if a fitted variable is tied in absolute terms (variance=0).
Allowed values: {1e-3 or above}; Default: 1e-3

32. minimum_amplitude_val: Minimum amplitude for emission line of any species given in log. 
Allowed values: {any value except 0}; Default: -1

33. maximum_amplitude_val: Maximum amplitude for emission line of any species given in log. 
Allowed values: {any value except 0}; Default: 6

34. minimum_reddening_val: Minimum reddening (E(B-V)) value for emission line fitting. Same reason as 31.
Allowed values: {1e-4 or above}; Default: 1e-4

35. maximum_reddening_val: Maximum reddening (E(B-V)) value for emission line fitting. This is important as, in low SNR spectra, the fit could suggest high reddening (garbage value) as this would be consistent with 0 flux in the model. Hence, it is important to check how much reddening can actually occur in each spaxel. 
Allowed values: {Any maximum reasonable reddening value}; Default: 2.0

36. sigma_bound_narrow_min: Minimum bound for velocity dispersion (in km/s) in the narrow line region of the emission line fit. This makes sure that the fit does not give unphysical dispersion values. 
Allowed values: {Of the order of the FWHM}; Default: 10.0

37. sigma_bound_narrow_max: Maximum bound for velocity dispersion (in km/s) in the narrow line region of the emission line fit. This makes sure that the fit does not give unphysical dispersion values. 
Allowed values: {Of the order of the FWHM}; Default: 200.0

38. sigma_bound_wide_max: Maximum bound for velocity dispersion (in km/s) in the wide line region of the emission line fit. This makes sure that the fit does not give unphysical dispersion values. 
Allowed values: {Of the order of the FWHM}; Default: 5000.0

39. poly_bound_val: +/- Value for binding the chebyshev polynomial fit for estimating the order for continuum. This bound is necessary for the customized chebyshev polynomial functions. 
Allowed values: {1 or above}; Default: 1000.0

40. maximum_accepted_reddening: Maximum accepted value for reddening (E(B-V)) value for emission line fitting. This is important as, in low SNR spectra, the fit could suggest high reddening (garbage value) as this would be consistent with 0 flux in the model. There is another parameter that sets the upper bound: maximum_reddening_val for curve fit, but that can be violated. This parameter forces reddening in the code to be below the mentioned values.
Allowed values: {1e-4 or above}; Default: 2.0

41. e_b_minus_v_init: Initial value for reddening (E(B-V)) in the emission line fitting.
Allowed values: {0.0 to maximum_reddening_val}; Default: 0.0

42. window_for_fit: Velocity window (in km/s) within which the fit for each single emission line species will take place. This ensures efficient fitting within limited time.
Allowed values: {10.0 or above}; Default: 2000.0

43. multiplicative_factor_for_plot: For efficient, time saving spectral line emission fitted plots, one can chose to bin the data separately (fitting would use all spectral points). This binning factor can be defined as the multiplicative_factor_for_plot.
Allowed values: {1 or above}; Default: 10

43.1 max_ew_width_abs: Maximum absolute equivalent width of lines in the galaxy spectra. This value will be used to limit failed equivalents width estimations from low SNR and bad spaxels.
Allowed values: {above 0.0}; Default: 100

43.2 allowed_ew_significance: Minimum absolute allowed equivalent width significance (data/err) of lines in the galaxy spectra. This value will be used to limit failed equivalents width estimations from low SNR and bad spaxels.
Allowed values: {above 0.0}; Default: 0.1

44. fit_velocity_val: Should the emission line fitting be done with fixed (False) or variable (True) central velocity?
Allowed values: {True or False}; Default: True

45. fit_vel_disp_val: Should the emission line fitting be done with fixed (False) or variable (True) central velocity dispersion?
Allowed values: {True or False}; Default: True

46. fit_continuum_val: Should the continuum be fitted again during custom emission line fitting? False will indicate an emission line fit with the initially estimated continuum from the method inspired by Martin+2021. 
Allowed values: {True or False}; Default: True

47. fit_reddening_val: Should the reddening be fitted again during custom emission line fitting? 
Allowed values: {True or False}; Default: True

48. plot_fit_refined: Should the emission line fitted results be plotted with binned spectral points (binned by multiplicative_factor_for_plot)? 
Allowed values: {True or False}; Default: False

49. ppxf_emission_tie_balmer: Tie Balmer line for ppxf emission fit for the entire cube
Allowed values: {True or False}; Default: True

50. ppxf_emission_limit_doublets: Limit doublets in the emission line fit for the entire cube
Allowed values: {True or False}; Default: True

51. center_init_array: Array for initial values (km/s) of central velocity for the emission line fitting
Allowed values: {Allowed values for central velocity}; Default: [array of values between -200 to 200]

52. sigma_init_array: Array for initial values (km/s) of central velocity dispersion for the emission line fitting
Allowed values: {Allowed values for central velocity dispersion}; Default: [array of center_init_array length with all values=100.0]

53. observed_instrument: Instrument used for observations
Allowed values: {This will depend on whether the function has been written to analyze the specific instrument data}; Default: VLT-MUSE

54. contour_count: Count for the number of contours to be plotted in MUSELI


EMISSION_FITTING_DETAILS:

There are additional keys that can be input for emission line fitting. The value for these keys are automatically estimated from the emission line file and using the get_emission_keys function from the get_data_from_file module. However these can be changed directly in the code as well. Although, we caution users to understand their importance before changing anything.  

1. center_list_names: This key provides the label for each lines from the emission line file. Default: species_name + int(lambda)
2. center_list: This key provides the lambda (rest wavelength) for each lines from the emission line file. 
3. comments_on_balmer: This key defines whether a specific line is a Balmer line (True) or not (False).
4. comments_on_tied: This key defines whether a specific line is tied to another. If so, the fit assumes a single tied line flux (rescaled according to line ratios).
5. factor_for_tied: This key defines the tied factor for tied lines.
6. factor_fixed_for_tying: This key indicates whether the tied factor can be fiddled around during the fit. If enabled, the fitting factor is also taken as another variable along with the specific line flux.



PYTHON CUSTOMIZED MODULES FOR MUSELI:

This section describes the different modules written for MUSELI. We will also describe the functions within the different modules at a later stage. There are short comments with the functions to describe what they do.

A. basic_functions: As the name suggests, this includes a lot of different basic functions that are used within the code. We explain the use of different functions below:
B. custom_functions_by_axe: This module includes all the main stage executing functions. Functions within this is called directly by the main code. Functions within this modules call the functions from all other python modules.
C. emission_fit_function: Module handling all the emission line fitting code (developed in-house)
D. extinction_fitting_functions: Module handling all the extinction curve fitting codes 
E. get_data_from_file: Module handling all the secondary data analysis like binning, refining, continuum estimation, data retrieval from different file types, saving data to different file types, sky spectra retrieval, initialising parameters and dictionary keys. 
F. handling_fits_files: Module handling all the fits file related functions (opening, modifying, saving).
G. plotting_functions: Module including all different functions used for plotting data.
H. ppxf_fitting_functions: Module including all advanced functions that call ppxf routines.
I. spatial_binning_functions: Module including all different routines to execute various types of spatial binning.
J. ppxf_customized directory; This directory contains all the ppxf modules from Michele Cappellari but modified just enough to fix some bugs with gas emission line fitting and adjusted for MUSELI. Due to these modifications, MUSELI includes a customized ppxf directory and calls this instead of the original ppxf available online.  
 

BASIC FILES:

This section lists are the basic (non-python) files included with MUSELI:

If some specific parameter needs to be changed, we recommend copying the file and changing within the copied version and calling it using the relevant dictionary key (See section.~PARAMETER FILE KEYS)

A. default_emission_lines.dat: Default file for emission line details. If you would like to use your own file, it has to be in this specific format.
B. default_lick_indexes.dat: Default file for lick index line details. If you would like to use your own file, it has to be in this specific format.
C. default_sky.dat: Default file for indicating sky spectra. If you would like to use your own file, it has to be in this specific format.
D. museli_env.yml: yml file for creating a conda environment within specific modules required for museli
E. museli_pip_requirements.txt: complementary txt file for installing all pip modules (not available using conda) related to museli
F. README.dat: Readme file
G. basic_parameter_file.dat: File showing all basic parameter keys (for additional convenience to users). Users can copy and modify this to their needs.


IMPORTANT NOTES:

1. 'clean_data' function defaults nan and empty elements as 1e-6 for data and 1e-5 for error. If this needs to be changed, please change it in 'basic_functions' module.
2. 'FWHM' value less than 4.0 (in km/s) gives error while running custom emission line fitting. For now, the code by default sets FWHM to 4.0 if the value drops before 4.0. We are investigating the cause for this and will hopefully solve this in the coming updates.
3. 


 
