;	start-up file for KL14 analysis

 !PROMPT = 'KL14> '

 PRINT,' '
 PRINT,'Loading all necessary procedures and functions....'


 COMMON numbers,tiny,huge,dynrange
 COMMON physconst,kboltz,ggravity,hplanck,clight,mhydrogen,sigmaboltz
 COMMON astrconst,kpc,pc,msun,lsun,yr,myr,kms

 tiny = 1.d-99
 huge = 9.d99
 dynrange = 1.d10

                                ; set physical constants (SI)

 kboltz = 1.3807d-23            ; Boltzmann constant	[J K-1]
 ggravity = 6.67259d-11         ; Gravity constant	[N m2 kg-2]
 hplanck = 6.62607d-34          ; Planck constant	[J s]
 clight = 2.99792458d8          ; Light speed		[m s-1]
 mhydrogen = 1.6726d-27         ; Hydrogen mass		[kg]
 sigmaboltz = 5.67d-8

                                ; astronomical constants

 kpc = 3.08572d19               ; 1000 parsec		[m]
 pc = kpc/1.d3
 msun = 1.989d30                ; solar mass		   [kg]
 lsun = 3.862d26                ; solar luminosity	[J s-1]
 yr = 86400.*365.25           ; sidereal year		[s]
 myr = 1.d6*yr
 kms = 1.d3
                                ; compile programs
.com derivefunc
.com clfind2d
.com clplot2d
.com clstats2d
.com smoothfits
.com peaks2d
.com derivefunc
.com fitKL14
.com unity_symmetric_sigmoidal
.com fourier_flux_loss_correction
.com lin_unity_intercept
.com fourier_zeta_correction
.com derivephys
.com plotfits
.com plotfits_files
.com strsplit
.com strnumber
.com astrometry_equal
.com mask_inside_circle
.com mask_inside_ellipse
.com mask_ds9_box_vertices
.com mask_ds9_file_convert
.com mask_ds9_file_mask
.com mask_tool.pro
.com is_an_integer
.com fourier_lowpass_ideal
.com fourier_highpass_ideal
.com fourier_lowpass_gaussian
.com fourier_highpass_gaussian
.com fourier_lowpass_butterworth
.com fourier_highpass_butterworth
.com fourier_filter_tool
.com fourier_diffuse_fraction
.com flux_fraction_calc
.com fsingle_array_make
.com fourier_power_fraction_vec
.com fourier_power_fraction
.com power_fraction_calc
.com get_clfind_peaks_filename
.com fourier_overlap_correction
.com med_peak_relative_nearest_neighbour_dist_calc

.com calc_levels
.com ds9_display_peaks
.com peak_var_read
.com peak_find
.com interactive_peak_find

.com f_error.pro
.com f_string.pro

 .com histoplot_update.pro


if keyword_set(COMMAND_LINE_ARGS()) then inputfile=COMMAND_LINE_ARGS() $
                                    else read,' please specify the full/absolute path of the input file (enter 0 to stop autorun): ',inputfile

.run tuningfork
