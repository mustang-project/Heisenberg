pro make_input_file, input_file_filepath   $ ; general variables
                   , datadir, galaxy, starfile, gasfile   $  ; File names compulsory variables

                   , distance, inclination, posangle, centrex, centrey, minradius, maxradius, Fs1_Fs2_min, max_sample, astr_tolerance, nbins $ ; # INPUT PARAMETERS 1 (basic map analysis)
                   , lapmin, lapmax, naperture, peak_res, max_res $  ; # INPUT PARAMETERS 2 (apertures)
                  ; ################################################################################################################################################################################################################################################################

                   , starfile2 = starfile2, gasfile2 = gasfile2, starfile3 = starfile3  $  ; File names keywords
                   , maskdir = maskdir, star_ext_mask1 = star_ext_mask, star_int_mask1 = star_int_mask, gas_ext_mask1 = gas_ext_mask, gas_int_mask1 = gas_int_mask, star_ext_mask2 = star_ext_mask2, star_int_mask2 = star_int_mask2, gas_ext_mask2 = gas_ext_mask2, gas_int_mask2 = gas_int_mask2, star_ext_mask3 = star_ext_mask3, star_int_mask3 = star_int_mask3  $  ; Mask file names keywords ; note star_ext_mask1 = star_ext_mask etc. prevents ambigious keyword error
                   , unfiltdir = unfiltdir, star_unfilt_file = star_unfilt_file, gas_unfilt_file = gas_unfilt_file $ ; unfiltered image keywords

                   , mask_images = mask_images, regrid = regrid, smoothen = smoothen, sensitivity = sensitivity, id_peaks = id_peaks, calc_ap_flux = calc_ap_flux, generate_plot = generate_plot, get_distances = get_distances, calc_obs = calc_obs, calc_fit = calc_fit, diffuse_frac = diffuse_frac, derive_phys = derive_phys, write_output = write_output, cleanup = cleanup, autoexit = autoexit $ ; FLAGS 1 (keywords)'
                   , use_star2 = use_star2, use_gas2 = use_gas2, use_star3 = use_star3 $ ;  FLAGS 2 keywords
                   , mstar_ext1 = mstar_ext, mstar_int1 = mstar_int, mgas_ext1 = mgas_ext, mgas_int1 = mgas_int, mstar_ext2 = mstar_ext2, mstar_int2 = mstar_int2, mgas_ext2 = mgas_ext2, mgas_int2 = mgas_int2, mstar_ext3 = mstar_ext3, mstar_int3 = mstar_int3, convert_masks = convert_masks, cut_radius = cut_radius $ ; # FLAGS 3 (masking-related options) keywords ; note mstar_ext1 = mstar_ext etc. prevents ambigious keyword error
                   , set_centre = set_centre, tophat = tophat, loglevels = loglevels, peak_find_tui = peak_find_tui, flux_weight = flux_weight, calc_ap_area = calc_ap_area, tstar_incl = tstar_incl, peak_prof = peak_prof, map_units = map_units, star_tot_mode = star_tot_mode, gas_tot_mode = gas_tot_mode,  use_X11 = use_X11, log10_output = log10_output $ ; # FLAGS 4 (choose analysis options)
                   , npixmin = npixmin, nsigma = nsigma, logrange_s = logrange_s, logspacing_s = logspacing_s, logrange_g = logrange_g, logspacing_g = logspacing_g, nlinlevel_s = nlinlevel_s, nlinlevel_g = nlinlevel_g $ ; # INPUT PARAMETERS 3 (peak identification)
                   , tstariso_val = tstariso, tstariso_errmin = tstariso_errmin, tstariso_errmax = tstariso_errmax, tgasmini = tgasmini, tgasmaxi = tgasmaxi, tovermini = tovermini $ ; # INPUT PARAMETERS 4 (timeline) ; note tstariso_val = tstariso prevents ambigious keyword error
                   , nmc = nmc, ndepth = ndepth, ntry = ntry, nphysmc = nphysmc $ ; # INPUT PARAMETERS 5 (fitting)
                   , use_unfilt_ims, diffuse_quant, f_filter_type, bw_order, filter_len_conv, emfrac_cor_mode, rpeak_cor_mode = rpeak_cor_mode, rpeaks_cor_val = rpeaks_cor_val, rpeaks_cor_emin = rpeaks_cor_emin, rpeaks_cor_emax = rpeaks_cor_emax, rpeakg_cor_val = rpeakg_cor_val, rpeakg_cor_emin = rpeakg_cor_emin, rpeakg_cor_emax = rpeakg_cor_emax $ ; # INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)
                   , convstar_val = convstar, convstar_rerr = convstar_rerr, convgas_val = convgas, convgas_rerr = convgas_rerr, convstar3_val = convstar3, convstar3_rerr = convstar3_rerr, lighttomass = lighttomass, momratetomass = momratetomass, star_tot_val = star_tot_val, star_tot_err = star_tot_err, gas_tot_val = gas_tot_val, gas_tot_err = gas_tot_err $; # INPUT PARAMETERS 6 (conversions and constants to calculate derived quantities) ; note convgas_val = convgas avoids % Ambiguous keyword abbreviation: CONVGAS.
                   , use_stds = use_stds, std_star1 = std_star1, std_star3 = std_star3, std_gas = std_gas $ ; # INPUT PARAMETERS 8 (sensitivity) ; star_std1 prevents ambigious keywords
                   , use_noisecut = use_noisecut, noisethresh_s = noisethresh_s, noisethresh_g = noisethresh_g $ ; # INPUT PARAMETERS 9 (noise threshold)
                   ; ################################################################################################################################################################################################################################################################
                   , use_guess = use_guess, initial_guess = initial_guess, iter_criterion = iter_criterion, iter_crit_len = iter_crit_len, iter_nmax = iter_nmax, iter_filter = iter_filter, iter_bwo = iter_bwo, iter_len_conv = iter_len_conv, iter_rpeak_mode = iter_rpeak_mode, iter_tot_mode_s = iter_tot_mode_s, iter_tot_mode_g = iter_tot_mode_g, iter_autoexit = iter_autoexit, use_nice = use_nice, nice_value = nice_value ; # INPUT PARAMETERS 9 (Fourier diffuse removal iteration) # note that these parameters are ignored if the call sequence idl kl14 -arg [full/absolute path of input file] is used. They are only used if the call sequence idl iterate_kl14 -arg [full/absolute path of input file] is used.



  openw, inp_lun, input_file_filepath, width = 450, /get_lun


  pad_string = ' '
  comment_sep_length = 9



  ; #############################################
  ; # FILE NAMES
  ; #############################################
  if n_elements(starfile2) ne 1 then starfile2 = '-'
  if n_elements(gasfile2) ne 1 then gasfile2 = '-'
  if n_elements(starfile3) ne 1 then starfile3 = '-'

  datadir_str = string(datadir)
  galaxy_str = string(galaxy)
  starfile_str = string(starfile)
  starfile2_str = string(starfile2)
  gasfile_str = string(gasfile)
  gasfile2_str = string(gasfile2)
  starfile3_str = string(starfile3)

  max_string_len = max([strlen(datadir_str), strlen(galaxy_str), strlen(starfile_str), strlen(gasfile_str), strlen(starfile2_str), strlen(gasfile2_str), strlen(starfile3_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(datadir_str) lt max_string_len) do datadir_str = datadir_str + pad_string
  while(strlen(galaxy_str) lt max_string_len) do galaxy_str = galaxy_str + pad_string
  while(strlen(starfile_str) lt max_string_len) do starfile_str = starfile_str + pad_string
  while(strlen(starfile2_str) lt max_string_len) do starfile2_str = starfile2_str + pad_string
  while(strlen(gasfile_str) lt max_string_len) do gasfile_str = gasfile_str + pad_string
  while(strlen(gasfile2_str) lt max_string_len) do gasfile2_str = gasfile2_str + pad_string
  while(strlen(starfile3_str) lt max_string_len) do starfile3_str = starfile3_str + pad_string





  ; #############################################
  ; # FLAGS 1 (keywords)'
  ; #############################################
  if n_elements(maskdir) ne 1 then maskdir = '-'
  if n_elements(star_ext_mask) ne 1 then star_ext_mask = '-'
  if n_elements(star_int_mask) ne 1 then star_int_mask = '-'
  if n_elements(gas_ext_mask) ne 1 then gas_ext_mask = '-'
  if n_elements(gas_int_mask) ne 1 then gas_int_mask = '-'
  if n_elements(star_ext_mask2) ne 1 then star_ext_mask2 = '-'
  if n_elements(star_int_mask2) ne 1 then star_int_mask2 = '-'
  if n_elements(gas_ext_mask2) ne 1 then gas_ext_mask2 = '-'
  if n_elements(gas_int_mask2) ne 1 then gas_int_mask2 = '-'
  if n_elements(star_ext_mask3) ne 1 then star_ext_mask3 = '-'
  if n_elements(star_int_mask3) ne 1 then star_int_mask3 = '-'

  maskdir_str = string(maskdir)
  star_ext_mask_str = string(star_ext_mask)
  star_int_mask_str = string(star_int_mask)
  gas_ext_mask_str = string(gas_ext_mask)
  gas_int_mask_str = string(gas_int_mask)
  star_ext_mask2_str = string(star_ext_mask2)
  star_int_mask2_str = string(star_int_mask2)
  gas_ext_mask2_str = string(gas_ext_mask2)
  gas_int_mask2_str = string(gas_int_mask2)
  star_ext_mask3_str = string(star_ext_mask3)
  star_int_mask3_str = string(star_int_mask3)

  max_string_len = max([strlen(maskdir_str), strlen(star_ext_mask_str), strlen(star_int_mask_str), strlen(gas_ext_mask_str), strlen(gas_int_mask_str), strlen(star_ext_mask2_str), strlen(star_int_mask2_str), strlen(gas_ext_mask2_str), strlen(gas_int_mask2_str), strlen(star_ext_mask3_str), strlen(star_int_mask3_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while( strlen(maskdir_str) lt max_string_len) do maskdir_str = maskdir_str + pad_string
  while( strlen(star_ext_mask_str) lt max_string_len) do star_ext_mask_str = star_ext_mask_str + pad_string
  while( strlen(star_int_mask_str) lt max_string_len) do star_int_mask_str = star_int_mask_str + pad_string
  while( strlen(gas_ext_mask_str) lt max_string_len) do gas_ext_mask_str = gas_ext_mask_str + pad_string
  while( strlen(gas_int_mask_str) lt max_string_len) do gas_int_mask_str = gas_int_mask_str + pad_string
  while( strlen(star_ext_mask2_str) lt max_string_len) do star_ext_mask2_str = star_ext_mask2_str + pad_string
  while( strlen(star_int_mask2_str) lt max_string_len) do star_int_mask2_str = star_int_mask2_str + pad_string
  while( strlen(gas_ext_mask2_str) lt max_string_len) do gas_ext_mask2_str = gas_ext_mask2_str + pad_string
  while( strlen(gas_int_mask2_str) lt max_string_len) do gas_int_mask2_str = gas_int_mask2_str + pad_string
  while( strlen(star_ext_mask3_str) lt max_string_len) do star_ext_mask3_str = star_ext_mask3_str + pad_string
  while( strlen(star_int_mask3_str) lt max_string_len) do star_int_mask3_str = star_int_mask3_str + pad_string


  ; #############################################
  ; UNFILTERED FILE NAMES
  ; #############################################
  if n_elements(unfiltdir) ne 1 then maskdir = '-'
  if n_elements(star_unfilt_file) ne 1 then star_unfilt_file = '-'
  if n_elements(gas_unfilt_file) ne 1 then gas_unfilt_file = '-'


  unfiltdir_str = string(unfiltdir)
  star_unfilt_file_str = string(star_unfilt_file)
  gas_unfilt_file_str = string(gas_unfilt_file)

  max_string_len = max([strlen(unfiltdir_str), strlen(star_unfilt_file_str), strlen(gas_unfilt_file_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while( strlen(unfiltdir_str) lt max_string_len) do unfiltdir_str = unfiltdir_str + pad_string
  while( strlen(star_unfilt_file_str) lt max_string_len) do star_unfilt_file_str = star_unfilt_file_str + pad_string
  while( strlen(gas_unfilt_file_str) lt max_string_len) do gas_unfilt_file_str = gas_unfilt_file_str + pad_string


  ; #############################################
  ; # FLAGS 1 (module switches)
  ; #############################################



  if n_elements(mask_images) ne 1 then mask_images = 1  ; defaults
  if n_elements(regrid) ne 1 then regrid = 1  ; defaults
  if n_elements(smoothen) ne 1 then smoothen = 1  ; defaults
  if n_elements(sensitivity) ne 1 then sensitivity = 1  ; defaults
  if n_elements(id_peaks) ne 1 then id_peaks = 1  ; defaults
  if n_elements(calc_ap_flux) ne 1 then calc_ap_flux = 1  ; defaults
  if n_elements(generate_plot) ne 1 then generate_plot = 1  ; defaults
  if n_elements(get_distances) ne 1 then get_distances = 1  ; defaults
  if n_elements(calc_obs) ne 1 then calc_obs = 1  ; defaults
  if n_elements(calc_fit) ne 1 then calc_fit = 1  ; defaults
  if n_elements(diffuse_frac) ne 1 then diffuse_frac = 1  ; defaults
  if n_elements(derive_phys) ne 1 then derive_phys = 1  ; defaults
  if n_elements(write_output) ne 1 then write_output = 1  ; defaults
  if n_elements(cleanup) ne 1 then cleanup = 1  ; defaults
  if n_elements(autoexit) ne 1 then autoexit = 0  ; defaults

  mask_images_str = strcompress(string(mask_images))
  regrid_str = strcompress(string(regrid))
  smoothen_str = strcompress(string(smoothen))
  sensitivity_str = strcompress(string(sensitivity))
  id_peaks_str = strcompress(string(id_peaks))
  calc_ap_flux_str = strcompress(string(calc_ap_flux))
  generate_plot_str = strcompress(string(generate_plot))
  get_distances_str = strcompress(string(get_distances))
  calc_obs_str = strcompress(string(calc_obs))
  calc_fit_str = strcompress(string(calc_fit))
  diffuse_frac_str = strcompress(string(diffuse_frac))
  derive_phys_str = strcompress(string(derive_phys))
  write_output_str = strcompress(string(write_output))
  cleanup_str = strcompress(string(cleanup))
  autoexit_str = strcompress(string(autoexit))


  max_string_len = max([strlen(mask_images_str), strlen(regrid_str), strlen(smoothen_str), strlen(sensitivity_str), strlen(id_peaks_str), strlen(calc_ap_flux_str), strlen(generate_plot_str), strlen(get_distances_str), strlen(calc_obs_str), strlen(calc_fit_str), strlen(derive_phys_str), strlen(write_output_str), strlen(cleanup_str), strlen(autoexit_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while (strlen(mask_images_str) lt max_string_len) do mask_images_str = mask_images_str + pad_string
  while (strlen(regrid_str) lt max_string_len) do regrid_str = regrid_str + pad_string
  while (strlen(smoothen_str) lt max_string_len) do smoothen_str = smoothen_str + pad_string
  while (strlen(sensitivity_str) lt max_string_len) do sensitivity_str = sensitivity_str + pad_string
  while (strlen(id_peaks_str) lt max_string_len) do id_peaks_str = id_peaks_str + pad_string
  while (strlen(calc_ap_flux_str) lt max_string_len) do calc_ap_flux_str = calc_ap_flux_str + pad_string
  while (strlen(generate_plot_str) lt max_string_len) do generate_plot_str = generate_plot_str + pad_string
  while (strlen(get_distances_str) lt max_string_len) do get_distances_str = get_distances_str + pad_string
  while (strlen(calc_obs_str) lt max_string_len) do calc_obs_str = calc_obs_str + pad_string
  while (strlen(calc_fit_str) lt max_string_len) do calc_fit_str = calc_fit_str + pad_string
  while (strlen(diffuse_frac_str) lt max_string_len) do diffuse_frac_str = diffuse_frac_str + pad_string
  while (strlen(derive_phys_str) lt max_string_len) do derive_phys_str = derive_phys_str + pad_string
  while (strlen(write_output_str) lt max_string_len) do write_output_str = write_output_str + pad_string
  while (strlen(cleanup_str) lt max_string_len) do cleanup_str = cleanup_str + pad_string
  while (strlen(autoexit_str) lt max_string_len) do autoexit_str = autoexit_str + pad_string


  ; #############################################
  ;'# FLAGS 2 (use ancillary files)''
  ; #############################################

  if n_elements(use_star2) ne 1 then use_star2 = 0
  if n_elements(use_gas2) ne 1 then use_gas2 = 0
  if n_elements(use_star3) ne 1 then use_star3 = 0

  use_star2_str = strcompress(string(use_star2))
  use_gas2_str = strcompress(string(use_gas2))
  use_star3_str = strcompress(string(use_star3))

  max_string_len = max([strlen(use_star2_str), strlen(use_gas2_str), strlen(use_star3_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(use_star2_str) lt max_string_len) do use_star2_str = use_star2_str + pad_string
  while(strlen(use_gas2_str) lt max_string_len) do use_gas2_str = use_gas2_str + pad_string
  while(strlen(use_star3_str) lt max_string_len) do use_star3_str = use_star3_str + pad_string




  ; #############################################
  ; '# FLAGS 3 (masking-related options)'
  ; #############################################

  if n_elements(mstar_ext) ne 1 then mstar_ext = 0
  if n_elements(mstar_int) ne 1 then mstar_int = 0
  if n_elements(mgas_ext) ne 1 then mgas_ext = 0
  if n_elements(mgas_int) ne 1 then mgas_int = 0
  if n_elements(mstar_ext2) ne 1 then mstar_ext2 = 0
  if n_elements(mstar_int2) ne 1 then mstar_int2 = 0
  if n_elements(mgas_ext2) ne 1 then mgas_ext2 = 0
  if n_elements(mgas_int2) ne 1 then mgas_int2 = 0
  if n_elements(mstar_ext3) ne 1 then mstar_ext3 = 0
  if n_elements(mstar_int3) ne 1 then mstar_int3 = 0
  if n_elements(convert_masks) ne 1 then convert_masks = 0
  if n_elements(cut_radius) ne 1 then cut_radius = 0

  mstar_ext_str = strcompress(string(mstar_ext))
  mstar_int_str = strcompress(string(mstar_int))
  mgas_ext_str = strcompress(string(mgas_ext))
  mgas_int_str = strcompress(string(mgas_int))
  mstar_ext2_str = strcompress(string(mstar_ext2))
  mstar_int2_str = strcompress(string(mstar_int2))
  mgas_ext2_str = strcompress(string(mgas_ext2))
  mgas_int2_str = strcompress(string(mgas_int2))
  mstar_ext3_str = strcompress(string(mstar_ext3))
  mstar_int3_str = strcompress(string(mstar_int3))
  convert_masks_str = strcompress(string(convert_masks))
  cut_radius_str = strcompress(string(cut_radius))

  max_string_len = max([strlen(mstar_ext_str), strlen(mstar_int_str), strlen(mgas_ext_str), strlen(mgas_int_str), strlen(mstar_ext2_str), strlen(mstar_int2_str), strlen(mgas_ext2_str), strlen(mgas_int2_str), strlen(mstar_ext3_str), strlen(mstar_int3_str), strlen(convert_masks_str), strlen(cut_radius_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(mstar_ext_str) lt max_string_len) do mstar_ext_str = mstar_ext_str + pad_string
  while(strlen(mstar_int_str) lt max_string_len) do mstar_int_str = mstar_int_str + pad_string
  while(strlen(mgas_ext_str) lt max_string_len) do mgas_ext_str = mgas_ext_str + pad_string
  while(strlen(mgas_int_str) lt max_string_len) do mgas_int_str = mgas_int_str + pad_string
  while(strlen(mstar_ext2_str) lt max_string_len) do mstar_ext2_str = mstar_ext2_str + pad_string
  while(strlen(mstar_int2_str) lt max_string_len) do mstar_int2_str = mstar_int2_str + pad_string
  while(strlen(mgas_ext2_str) lt max_string_len) do mgas_ext2_str = mgas_ext2_str + pad_string
  while(strlen(mgas_int2_str) lt max_string_len) do mgas_int2_str = mgas_int2_str + pad_string
  while(strlen(mstar_ext3_str) lt max_string_len) do mstar_ext3_str = mstar_ext3_str + pad_string
  while(strlen(mstar_int3_str) lt max_string_len) do mstar_int3_str = mstar_int3_str + pad_string
  while(strlen(convert_masks_str) lt max_string_len) do convert_masks_str = convert_masks_str + pad_string
  while(strlen(cut_radius_str) lt max_string_len) do cut_radius_str = cut_radius_str + pad_string




  ; #############################################
  ; # FLAGS 4 (choose analysis options)
  ; #############################################

  if n_elements(set_centre) ne 1 then set_centre = 1
  if n_elements(tophat) ne 1 then tophat = 1
  if n_elements(loglevels) ne 1 then loglevels = 1
  if n_elements(peak_find_tui) ne 1 then peak_find_tui = 0
  if n_elements(flux_weight) ne 1 then flux_weight = 0
  if n_elements(calc_ap_area) ne 1 then calc_ap_area = 1
  if n_elements(tstar_incl) ne 1 then tstar_incl = 0
  if n_elements(peak_prof) ne 1 then peak_prof = 2
  if n_elements(map_units) ne 1 then map_units = 1
  if n_elements(star_tot_mode) ne 1 then star_tot_mode = 0
  if n_elements(gas_tot_mode) ne 1 then star_tot_mode = 0
  if n_elements(use_X11) ne 1 then use_X11 = 1
  if n_elements(log10_output) ne 1 then log10_output = 1

  set_centre_str = strcompress(string(set_centre))
  tophat_str = strcompress(string(tophat))
  loglevels_str = strcompress(string(loglevels))
  peak_find_tui_str = strcompress(string(peak_find_tui))
  flux_weight_str = strcompress(string(flux_weight))
  calc_ap_area_str = strcompress(string(calc_ap_area))
  tstar_incl_str = strcompress(string(tstar_incl))
  peak_prof_str = strcompress(string(peak_prof))
  map_units_str = strcompress(string(map_units))
  star_tot_mode_str = strcompress(string(star_tot_mode))
  gas_tot_mode_str = strcompress(string(gas_tot_mode))
  use_X11_str = strcompress(string(use_X11))
  log10_output_str = strcompress(string(log10_output))


  max_string_len = max([strlen(set_centre_str), strlen(tophat_str), strlen(loglevels_str), strlen(peak_find_tui_str), strlen(flux_weight_str), strlen(calc_ap_area_str), strlen(tstar_incl_str), strlen(peak_prof_str), strlen(map_units_str), strlen(star_tot_mode_str), strlen(gas_tot_mode_str), strlen(use_X11_str), strlen(log10_output_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(set_centre_str) lt max_string_len) do set_centre_str = set_centre_str + pad_string
  while(strlen(tophat_str) lt max_string_len) do tophat_str = tophat_str + pad_string
  while(strlen(loglevels_str) lt max_string_len) do loglevels_str = loglevels_str + pad_string
  while(strlen(peak_find_tui_str) lt max_string_len) do peak_find_tui_str = peak_find_tui_str + pad_string
  while(strlen(flux_weight_str) lt max_string_len) do flux_weight_str = flux_weight_str + pad_string
  while(strlen(calc_ap_area_str) lt max_string_len) do calc_ap_area_str = calc_ap_area_str + pad_string
  while(strlen(tstar_incl_str) lt max_string_len) do tstar_incl_str = tstar_incl_str + pad_string
  while(strlen(peak_prof_str) lt max_string_len) do peak_prof_str = peak_prof_str + pad_string
  while(strlen(map_units_str) lt max_string_len) do map_units_str = map_units_str + pad_string
  while(strlen(star_tot_mode_str) lt max_string_len) do star_tot_mode_str = star_tot_mode_str + pad_string
  while(strlen(gas_tot_mode_str) lt max_string_len) do gas_tot_mode_str = gas_tot_mode_str + pad_string
  while(strlen(use_X11_str) lt max_string_len) do use_X11_str = use_X11_str + pad_string
  while(strlen(log10_output_str) lt max_string_len) do log10_output_str = log10_output_str + pad_string










  ; #############################################
  ; # INPUT PARAMETERS 1 (basic map analysis)
  ; #############################################
  if n_elements(distance) ne 1 then stop
  if n_elements(inclination) ne 1 then stop
  if n_elements(posangle) ne 1 then stop
  if n_elements(centrex) ne 1 then stop
  if n_elements(centrey) ne 1 then stop
  if n_elements(minradius) ne 1 then stop
  if n_elements(maxradius) ne 1 then stop
  if n_elements(Fs1_Fs2_min) ne 1 then stop
  if n_elements(max_sample) ne 1 then stop
  if n_elements(astr_tolerance) ne 1 then stop
  if n_elements(nbins) ne 1 then stop


  distance_str = strcompress(string(distance))
  inclination_str = strcompress(string(inclination))
  posangle_str = strcompress(string(posangle))
  centrex_str = strcompress(string(centrex))
  centrey_str = strcompress(string(centrey))
  minradius_str = strcompress(string(minradius))
  maxradius_str = strcompress(string(maxradius))
  Fs1_Fs2_min_str = strcompress(string(Fs1_Fs2_min))
  max_sample_str = strcompress(string(max_sample))
  astr_tolerance_str = strcompress(string(astr_tolerance))
  nbins_str = strcompress(string(nbins))

  max_string_len = max([strlen(distance_str), strlen(inclination_str), strlen(posangle_str), strlen(centrex_str), strlen(centrey_str), strlen(minradius_str), strlen(maxradius_str), strlen(Fs1_Fs2_min_str), strlen(max_sample_str), strlen(astr_tolerance_str), strlen(nbins_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(distance_str) lt max_string_len) do distance_str = distance_str + pad_string
  while(strlen(inclination_str) lt max_string_len) do inclination_str = inclination_str + pad_string
  while(strlen(posangle_str) lt max_string_len) do posangle_str = posangle_str + pad_string
  while(strlen(centrex_str) lt max_string_len) do centrex_str = centrex_str + pad_string
  while(strlen(centrey_str) lt max_string_len) do centrey_str = centrey_str + pad_string
  while(strlen(minradius_str) lt max_string_len) do minradius_str = minradius_str + pad_string
  while(strlen(maxradius_str) lt max_string_len) do maxradius_str = maxradius_str + pad_string
  while(strlen(Fs1_Fs2_min_str) lt max_string_len) do Fs1_Fs2_min_str = Fs1_Fs2_min_str + pad_string
  while(strlen(max_sample_str) lt max_string_len) do max_sample_str = max_sample_str + pad_string
  while(strlen(astr_tolerance_str) lt max_string_len) do astr_tolerance_str = astr_tolerance_str + pad_string
  while(strlen(nbins_str) lt max_string_len) do nbins_str = nbins_str + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 2 (apertures)
  ; #############################################

  if n_elements(lapmin) ne 1 then stop
  if n_elements(lapmax) ne 1 then stop
  if n_elements(naperture) ne 1 then stop
  if n_elements(peak_res) ne 1 then stop
  if n_elements(max_res) ne 1 then stop

  lapmin_str = strcompress(string(lapmin))
  lapmax_str = strcompress(string(lapmax))
  naperture_str = strcompress(string(naperture))
  peak_res_str = strcompress(string(peak_res))
  max_res_str = strcompress(string(max_res))

  max_string_len = max([strlen(lapmin_str), strlen(lapmax_str), strlen(naperture_str), strlen(peak_res_str), strlen(max_res_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(lapmin_str) lt max_string_len) do lapmin_str = lapmin_str + pad_string
  while(strlen(lapmax_str) lt max_string_len) do lapmax_str = lapmax_str + pad_string
  while(strlen(naperture_str) lt max_string_len) do naperture_str = naperture_str + pad_string
  while(strlen(peak_res_str) lt max_string_len) do peak_res_str = peak_res_str + pad_string
  while(strlen(max_res_str) lt max_string_len) do max_res_str = max_res_str + pad_string

  ; #############################################
  ; # INPUT PARAMETERS 3 (peak identification)
  ; #############################################
  if n_elements(npixmin) ne 1 then npixmin = 20
  if n_elements(nsigma) ne 1 then nsigma = 5
  if n_elements(logrange_s) ne 1 then logrange_s = 2.
  if n_elements(logspacing_s) ne 1 then logspacing_s = 0.5
  if n_elements(logrange_g) ne 1 then logrange_g = 2.
  if n_elements(logspacing_g) ne 1 then logspacing_g = 0.5
  if n_elements(nlinlevel_s) ne 1 then nlinlevel_s = 11
  if n_elements(nlinlevel_g) ne 1 then nlinlevel_g = 11

  npixmin_str = strcompress(string(npixmin))
  nsigma_str = strcompress(string(nsigma))
  logrange_s_str = strcompress(string(logrange_s))
  logspacing_s_str = strcompress(string(logspacing_s))
  logrange_g_str = strcompress(string(logrange_g))
  logspacing_g_str = strcompress(string(logspacing_g))
  nlinlevel_s_str = strcompress(string(nlinlevel_s))
  nlinlevel_g_str = strcompress(string(nlinlevel_g))

  max_string_len = max([strlen(npixmin_str), strlen(nsigma_str), strlen(logrange_s_str), strlen(logspacing_s_str), strlen(logrange_g_str), strlen(logspacing_g_str), strlen(nlinlevel_s_str), strlen(nlinlevel_g_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(npixmin_str) lt max_string_len) do npixmin_str = npixmin_str + pad_string
  while(strlen(nsigma_str) lt max_string_len) do nsigma_str = nsigma_str + pad_string
  while(strlen(logrange_s_str) lt max_string_len) do logrange_s_str = logrange_s_str + pad_string
  while(strlen(logspacing_s_str) lt max_string_len) do logspacing_s_str = logspacing_s_str + pad_string
  while(strlen(logrange_g_str) lt max_string_len) do logrange_g_str = logrange_g_str + pad_string
  while(strlen(logspacing_g_str) lt max_string_len) do logspacing_g_str = logspacing_g_str + pad_string
  while(strlen(nlinlevel_s_str) lt max_string_len) do nlinlevel_s_str = nlinlevel_s_str + pad_string
  while(strlen(nlinlevel_g_str) lt max_string_len) do nlinlevel_g_str = nlinlevel_g_str + pad_string

  ; #############################################
  ; # INPUT PARAMETERS 4 (timeline)
  ; #############################################

  if  n_elements(tstariso) ne 1 then tstariso = 1.
  if  n_elements(tstariso_errmin) ne 1 then tstariso_errmin = 0.
  if  n_elements(tstariso_errmax) ne 1 then tstariso_errmax = 0.
  if  n_elements(tgasmini) ne 1 then tgasmini = 0.1
  if  n_elements(tgasmaxi) ne 1 then tgasmaxi = 1000.
  if  n_elements(tovermini) ne 1 then tovermini = 0.01

  tstariso_str = strcompress(string(tstariso))
  tstariso_errmin_str = strcompress(string(tstariso_errmin))
  tstariso_errmax_str = strcompress(string(tstariso_errmax))
  tgasmini_str = strcompress(string(tgasmini))
  tgasmaxi_str = strcompress(string(tgasmaxi))
  tovermini_str = strcompress(string(tovermini))

  max_string_len = max([strlen(tstariso_str), strlen(tstariso_errmin_str), strlen(tstariso_errmax_str), strlen(tgasmini_str), strlen(tgasmaxi_str), strlen(tovermini_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(tstariso_str) lt max_string_len) do tstariso_str = tstariso_str + pad_string
  while(strlen(tstariso_errmin_str) lt max_string_len) do tstariso_errmin_str = tstariso_errmin_str + pad_string
  while(strlen(tstariso_errmax_str) lt max_string_len) do tstariso_errmax_str = tstariso_errmax_str + pad_string
  while(strlen(tgasmini_str) lt max_string_len) do tgasmini_str = tgasmini_str + pad_string
  while(strlen(tgasmaxi_str) lt max_string_len) do tgasmaxi_str = tgasmaxi_str + pad_string
  while(strlen(tovermini_str) lt max_string_len) do tovermini_str = tovermini_str + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 5 (fitting)
  ; #############################################
  if n_elements(nmc) ne 1 then nmc = 1000
  if n_elements(ndepth) ne 1 then ndepth = 4
  if n_elements(ntry) ne 1 then ntry = 101
  if n_elements(nphysmc) ne 1 then nphysmc = 1000000

  nmc_str = strcompress(string(nmc))
  ndepth_str = strcompress(string(ndepth))
  ntry_str = strcompress(string(ntry))
  nphysmc_str = strcompress(string(nphysmc))

  max_string_len = max( [strlen(nmc_str), strlen(ndepth_str), strlen(ntry_str), strlen(nphysmc_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(nmc_str) lt max_string_len) do nmc_str = nmc_str + pad_string
  while(strlen(ndepth_str) lt max_string_len) do ndepth_str = ndepth_str + pad_string
  while(strlen(ntry_str) lt max_string_len) do ntry_str = ntry_str + pad_string
  while(strlen(nphysmc_str) lt max_string_len) do nphysmc_str = nphysmc_str + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)
  ; #############################################
  if n_elements(use_unfilt_ims) ne 1 then use_unfilt_ims = 0
  if n_elements(diffuse_quant) ne 1 then diffuse_quant = 1
  if n_elements(f_filter_type) ne 1 then f_filter_type = 2
  if n_elements(bw_order) ne 1 then bw_order = 2
  if n_elements(filter_len_conv) ne 1 then filter_len_conv = 1.0
  if n_elements(emfrac_cor_mode) ne 1 then emfrac_cor_mode = 0

  if n_elements(rpeak_cor_mode) ne 1 then rpeak_cor_mode = 0
  if n_elements(rpeaks_cor_val) ne 1 then rpeaks_cor_val = 1.0
  if n_elements(rpeaks_cor_emin) ne 1 then rpeaks_cor_emin = 0.5
  if n_elements(rpeaks_cor_emax) ne 1 then rpeaks_cor_emax = 0.5
  if n_elements(rpeakg_cor_val) ne 1 then rpeakg_cor_val = 1.0
  if n_elements(rpeakg_cor_emin) ne 1 then rpeakg_cor_emin = 0.5
  if n_elements(rpeakg_cor_emax) ne 1 then rpeakg_cor_emax = 0.5



  use_unfilt_ims_str = strcompress(string(use_unfilt_ims))
  diffuse_quant_str = strcompress(string(diffuse_quant))
  f_filter_type_str = strcompress(string(f_filter_type))
  bw_order_str = strcompress(string(bw_order))
  filter_len_conv_str = strcompress(string(filter_len_conv))
  emfrac_cor_mode_str = strcompress(string(emfrac_cor_mode))
  rpeak_cor_mode_str = strcompress(string(rpeak_cor_mode))
  rpeaks_cor_val_str = strcompress(string(rpeaks_cor_val))
  rpeaks_cor_emin_str = strcompress(string(rpeaks_cor_emin))
  rpeaks_cor_emax_str = strcompress(string(rpeaks_cor_emax))
  rpeakg_cor_val_str = strcompress(string(rpeakg_cor_val))
  rpeakg_cor_emin_str = strcompress(string(rpeakg_cor_emin))
  rpeakg_cor_emax_str = strcompress(string(rpeakg_cor_emax))

  max_string_len = max([strlen(use_unfilt_ims_str),strlen(diffuse_quant_str), strlen(f_filter_type_str), strlen(bw_order_str), strlen(filter_len_conv_str), strlen(emfrac_cor_mode), strlen(rpeak_cor_mode_str), strlen(rpeaks_cor_val_str), strlen(rpeaks_cor_emin_str), strlen(rpeaks_cor_emax_str), strlen(rpeakg_cor_val_str), strlen(rpeakg_cor_emin_str), strlen(rpeakg_cor_emax_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(use_unfilt_ims_str) lt max_string_len) do use_unfilt_ims_str = use_unfilt_ims_str + pad_string
  while(strlen(diffuse_quant_str) lt max_string_len) do diffuse_quant_str = diffuse_quant_str + pad_string
  while(strlen(f_filter_type_str) lt max_string_len) do f_filter_type_str = f_filter_type_str + pad_string
  while(strlen(bw_order_str) lt max_string_len) do bw_order_str = bw_order_str + pad_string
  while(strlen(filter_len_conv_str) lt max_string_len) do filter_len_conv_str = filter_len_conv_str + pad_string
  while(strlen(emfrac_cor_mode_str) lt max_string_len) do emfrac_cor_mode_str = emfrac_cor_mode_str + pad_string
  while(strlen(rpeak_cor_mode_str) lt max_string_len) do rpeak_cor_mode_str = rpeak_cor_mode_str + pad_string
  while(strlen(rpeaks_cor_val_str) lt max_string_len) do rpeaks_cor_val_str = rpeaks_cor_val_str + pad_string
  while(strlen(rpeaks_cor_emin_str) lt max_string_len) do rpeaks_cor_emin_str = rpeaks_cor_emin_str + pad_string
  while(strlen(rpeaks_cor_emax_str) lt max_string_len) do rpeaks_cor_emax_str = rpeaks_cor_emax_str + pad_string
  while(strlen(rpeakg_cor_val_str) lt max_string_len) do rpeakg_cor_val_str = rpeakg_cor_val_str + pad_string
  while(strlen(rpeakg_cor_emin_str) lt max_string_len) do rpeakg_cor_emin_str = rpeakg_cor_emin_str + pad_string
  while(strlen(rpeakg_cor_emax_str) lt max_string_len) do rpeakg_cor_emax_str = rpeakg_cor_emax_str + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 7 (conversions and constants to calculate derived quantities)
  ; #############################################
  if n_elements(convstar) ne 1 then convstar =  -3.69206
  if n_elements(convstar_rerr) ne 1 then convstar_rerr =  0.
  if n_elements(convgas) ne 1 then convgas =  2.30794
  if n_elements(convgas_rerr) ne 1 then convgas_rerr =  0.
  if n_elements(convstar3) ne 1 then convstar3 =  0.
  if n_elements(convstar3_rerr) ne 1 then convstar3_rerr =  0.
  if n_elements(lighttomass) ne 1 then lighttomass =  0.002
  if n_elements(momratetomass) ne 1 then momratetomass =  5.e-10
  if n_elements(star_tot_val) ne 1 then star_tot_val =  1.0
  if n_elements(star_tot_err) ne 1 then star_tot_err =  0.1
  if n_elements(gas_tot_val) ne 1 then gas_tot_val =  1.0e9
  if n_elements(gas_tot_err) ne 1 then gas_tot_err =  1.0e8

  convstar_str = strcompress(string(convstar))
  convstar_rerr_str = strcompress(string(convstar_rerr))
  convgas_str = strcompress(string(convgas))
  convgas_rerr_str = strcompress(string(convgas_rerr))
  convstar3_str = strcompress(string(convstar3))
  convstar3_rerr_str = strcompress(string(convstar3_rerr))
  lighttomass_str = strcompress(string(lighttomass))
  momratetomass_str = strcompress(string(momratetomass))
  star_tot_val_str = strcompress(string(star_tot_val))
  star_tot_err_str = strcompress(string(star_tot_err))
  gas_tot_val_str = strcompress(string(gas_tot_val))
  gas_tot_err_str = strcompress(string(gas_tot_err))

  max_string_len = max([strlen(convstar_str), strlen(convstar_rerr_str), strlen(convgas_str), strlen(convgas_rerr_str), strlen(convstar3_str), strlen(convstar3_rerr_str), strlen(lighttomass_str), strlen(momratetomass_str), strlen(star_tot_val_str), strlen(star_tot_err_str), strlen(gas_tot_val_str), strlen(gas_tot_err_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(convstar_str) lt max_string_len) do convstar_str = convstar_str + pad_string
  while(strlen(convstar_rerr_str) lt max_string_len) do convstar_rerr_str = convstar_rerr_str + pad_string
  while(strlen(convgas_str) lt max_string_len) do convgas_str = convgas_str + pad_string
  while(strlen(convgas_rerr_str) lt max_string_len) do convgas_rerr_str = convgas_rerr_str + pad_string
  while(strlen(convstar3_str) lt max_string_len) do convstar3_str = convstar3_str + pad_string
  while(strlen(convstar3_rerr_str) lt max_string_len) do convstar3_rerr_str = convstar3_rerr_str + pad_string
  while(strlen(lighttomass_str) lt max_string_len) do lighttomass_str = lighttomass_str + pad_string
  while(strlen(momratetomass_str) lt max_string_len) do momratetomass_str = momratetomass_str + pad_string
  while(strlen(star_tot_val_str) lt max_string_len) do star_tot_val_str = star_tot_val_str + pad_string
  while(strlen(star_tot_err_str) lt max_string_len) do star_tot_err_str = star_tot_err_str + pad_string
  while(strlen(gas_tot_val_str) lt max_string_len) do gas_tot_val_str = gas_tot_val_str + pad_string
  while(strlen(gas_tot_err_str) lt max_string_len) do gas_tot_err_str = gas_tot_err_str + pad_string



  ; #############################################
  ; # INPUT PARAMETERS 8 (sensitivity)
  ; #############################################


  if n_elements(use_stds) ne 1 then use_stds =  0
  if n_elements(std_star1) ne 1 then std_star1 =  0.1
  if n_elements(std_star3) ne 1 then std_star3 =  0.1
  if n_elements(std_gas) ne 1 then std_gas =  0.1

  use_stds_str = strcompress(string(use_stds))
  std_star_str = strcompress(string(std_star1))
  std_star3_str = strcompress(string(std_star3))
  std_gas_str = strcompress(string(std_gas))

  max_string_len = max([strlen(use_stds_str), strlen(std_star_str), strlen(std_star3_str), strlen(std_gas_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(use_stds_str) lt max_string_len) do use_stds_str = use_stds_str + pad_string
  while(strlen(std_star_str) lt max_string_len) do std_star_str = std_star_str + pad_string
  while(strlen(std_star3_str) lt max_string_len) do std_star3_str = std_star3_str + pad_string
  while(strlen(std_gas_str) lt max_string_len) do std_gas_str = std_gas_str + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 9 (noise threshold)
  ; #############################################


  if n_elements(use_noisecut) ne 1 then use_noisecut =  0
  if n_elements(noisethresh_s) ne 1 then noisethresh_s =  0.1
  if n_elements(noisethresh_g) ne 1 then noisethresh_g =  0.1

  use_noisecut_str = strcompress(string(use_noisecut))
  noisethresh_s_str = strcompress(string(noisethresh_s))
  noisethresh_g_str = strcompress(string(noisethresh_g))

  max_string_len = max([strlen(use_noisecut_str), strlen(noisethresh_s_str), strlen(noisethresh_g_str)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(use_noisecut_str) lt max_string_len) do use_noisecut_str = use_noisecut_str + pad_string
  while(strlen(noisethresh_s_str) lt max_string_len) do noisethresh_s_str = noisethresh_s_str + pad_string
  while(strlen(noisethresh_g_str) lt max_string_len) do noisethresh_g_str = noisethresh_g_str + pad_string



  ; ##########################################################################################
  ; # INPUT PARAMETERS 9 (Fourier diffuse removal iteration) -- Optional
  ; ##########################################################################################
  iter_input_switch = 0
  if n_elements(use_guess) eq 1 && n_elements(initial_guess) eq 1 && n_elements(iter_criterion) eq 1 && n_elements(iter_crit_len) eq 1 && n_elements(iter_nmax) eq 1 && n_elements(iter_filter) eq 1 && n_elements(iter_bwo) eq 1 && n_elements(iter_len_conv) eq 1 && n_elements(iter_rpeak_mode) eq 1 && n_elements(iter_tot_mode_s) eq 1 && n_elements(iter_tot_mode_g) eq 1  && n_elements(iter_autoexit) eq 1 && n_elements(use_nice) eq 1 && n_elements(nice_value) eq 1 then begin
    iter_input_switch = 1 ; use so don't need to repeat test


    use_guess_str = strcompress(string(use_guess))
    initial_guess_str = strcompress(string(initial_guess))
    iter_criterion_str = strcompress(string(iter_criterion))
    iter_crit_len_str = strcompress(string(iter_crit_len))
    iter_nmax_str = strcompress(string(iter_nmax))
    iter_filter_str = strcompress(string(iter_filter))
    iter_bwo_str = strcompress(string(iter_bwo))
    iter_len_conv_str = strcompress(string(iter_len_conv))
    iter_rpeak_mode_str = strcompress(string(iter_rpeak_mode))
    iter_tot_mode_s_str = strcompress(string(iter_tot_mode_s))
    iter_tot_mode_g_str = strcompress(string(iter_tot_mode_g))
    iter_autoexit_str = strcompress(string(iter_autoexit))
    use_nice_str = strcompress(string(use_nice))
    nice_value_str = strcompress(string(nice_value))



    max_string_len = max([strlen(use_guess_str), strlen(initial_guess_str), strlen(iter_criterion_str), strlen(iter_crit_len_str), strlen(iter_nmax_str), strlen(iter_filter_str), strlen(iter_bwo_str), strlen(iter_len_conv_str), strlen(iter_tot_mode_s_str), strlen(iter_tot_mode_g_str), strlen(iter_rpeak_mode_str), strlen(iter_autoexit_str), strlen(max_string_len), strlen(use_nice_str), strlen(nice_value_str)])
    max_string_len = max([max_string_len,22])


    while(strlen(use_guess_str) lt max_string_len) do use_guess_str = use_guess_str + pad_string
    while(strlen(initial_guess_str) lt max_string_len) do initial_guess_str = initial_guess_str + pad_string
    while(strlen(iter_criterion_str) lt max_string_len) do iter_criterion_str = iter_criterion_str + pad_string
    while(strlen(iter_crit_len_str) lt max_string_len) do iter_crit_len_str = iter_crit_len_str + pad_string
    while(strlen(iter_nmax_str) lt max_string_len) do iter_nmax_str = iter_nmax_str + pad_string
    while(strlen(iter_filter_str) lt max_string_len) do iter_filter_str = iter_filter_str + pad_string
    while(strlen(iter_bwo_str) lt max_string_len) do iter_bwo_str = iter_bwo_str + pad_string
    while(strlen(iter_len_conv_str) lt max_string_len) do iter_len_conv_str = iter_len_conv_str + pad_string
    while(strlen(iter_rpeak_mode_str) lt max_string_len) do iter_rpeak_mode_str = iter_rpeak_mode_str + pad_string
    while(strlen(iter_tot_mode_s_str) lt max_string_len) do iter_tot_mode_s_str = iter_tot_mode_s_str + pad_string
    while(strlen(iter_tot_mode_g_str) lt max_string_len) do iter_tot_mode_g_str = iter_tot_mode_g_str + pad_string
    while(strlen(iter_autoexit_str) lt max_string_len) do iter_autoexit_str = iter_autoexit_str + pad_string
    while(strlen(use_nice_str) lt max_string_len) do use_nice_str = use_nice_str + pad_string
    while(strlen(nice_value_str) lt max_string_len) do nice_value_str = nice_value_str + pad_string



  endif


  ; #############################################
  ; # write input file
  ; #############################################


  printf, inp_lun, '################################################################################'
  printf, inp_lun, '#                                                                              #'
  printf, inp_lun, '#                  FIT KL14 PRINCIPLE TO OBSERVED GALAXY MAPS                  #'
  printf, inp_lun, '#  start environment with >> idl kl14 -arg [full/absolute path of input file]  #'
  printf, inp_lun, '#                                                                              #'
  if iter_input_switch eq 1 then begin ; print instructions for iterative running
    printf, inp_lun, '#  for iterative diffuse filtering, start environment with:                    #'
    printf, inp_lun, '#   >> idl iterate_kl14 -arg [full/absolute path of input file]                #'
    printf, inp_lun, '#                                                                              #'
  endif
  printf, inp_lun, '#                             BEGIN PARAMETER FILE                             #'
  printf, inp_lun, '#                                                                              #'
  printf, inp_lun, '################################################################################'
  printf, inp_lun, ''
  printf, inp_lun, ''

  printf, inp_lun, '# FILE NAMES'
  printf, inp_lun, 'datadir         ', datadir_str                                                                , '# Full path of data directory'
  printf, inp_lun, 'galaxy          ', galaxy_str                                                             , '# Name of data set'
  printf, inp_lun, 'starfile        ', starfile_str                                                                , '# Name of primary stellar map'
  printf, inp_lun, 'starfile2       ', starfile2_str                                                               , '# Name of secondary stellar map (specifically for peak identification - only used if use_star2=1)'
  printf, inp_lun, 'gasfile         ', gasfile_str                                                                 , '# Name of primary gas map'
  printf, inp_lun, 'gasfile2        ', gasfile2_str                                                                       , '# Name of secondary gas map (specifically for peak identification - only used if use_gas2=1)'
  printf, inp_lun, 'starfile3       ', starfile3_str                                                                       , '# Name of additional stellar map (only used if use_star3=1)'

  printf, inp_lun, '# MASK FILE NAMES'
  printf, inp_lun, 'maskdir            ', maskdir_str,           '# Full path of mask directory'
  printf, inp_lun, 'star_ext_mask      ', star_ext_mask_str,     '# Region file for masking areas in starfile exterior to closed shapes (only used if mask_star=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask      ', star_int_mask_str,     '# Region file for masking areas in starfile interior to closed shapes (only used if mask_star=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_ext_mask       ', gas_ext_mask_str,      '# Region file for masking areas in gasfile exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_int_mask       ', gas_int_mask_str,      '# Region file for masking areas in gasfile interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_ext_mask2     ', star_ext_mask2_str,    '# Region file for masking areas in starfile2 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask2     ', star_int_mask2_str,    '# Region file for masking areas in starfile2 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_ext_mask2      ', gas_ext_mask2_str,     '# Region file for masking areas in gasfile2 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_int_mask2      ', gas_int_mask2_str,     '# Region file for masking areas in gasfile2 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_ext_mask3     ', star_ext_mask3_str,    '# Region file for masking areas in starfile3 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask3     ', star_int_mask3_str,    '# Region file for masking areas in starfile3 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'


  printf, inp_lun, '# UNFILTERED FILE NAMES'
  printf, inp_lun, 'unfiltdir          ', unfiltdir_str,         '# Full path of unfiltered image directory'
  printf, inp_lun, 'star_unfilt_file   ', star_unfilt_file_str
  printf, inp_lun, 'gas_unfilt_file    ', gas_unfilt_file_str


  printf, inp_lun, '# FLAGS 1 (module switches)'
  printf, inp_lun, 'mask_images     ', mask_images_str,    '# Mask images (on/off)'
  printf, inp_lun, 'regrid          ', regrid_str,         '# Read and regrid original files and synchronise masks across all used images, i.e. pixels masked in one image will be masked in all images (on/off)'
  printf, inp_lun, 'smoothen        ', smoothen_str,       '# Create smoothened maps for each aperture size (on/off)'
  printf, inp_lun, 'sensitivity     ', sensitivity_str,    '# Fit Gaussians to pixel intensity PDFs to determine sensitivity limits (on/off)'
  printf, inp_lun, 'id_peaks        ', id_peaks_str,       '# Identify peaks and save their locations (on/off)'
  printf, inp_lun, 'calc_ap_flux    ', calc_ap_flux_str,   '# Calculate the enclosed flux for each peak and aperture size (on/off)'
  printf, inp_lun, 'generate_plot   ', generate_plot_str,  '# Generate output plots and save to PostScript files (on/off)'
  printf, inp_lun, 'get_distances   ', get_distances_str,  '# Calculate distances between all peak pairs (on/off)'
  printf, inp_lun, 'calc_obs        ', calc_obs_str,       '# Monte-Carlo sample peaks, get observed tuning fork (on/off)'
  printf, inp_lun, 'calc_fit        ', calc_fit_str,       '# Fit KL14 principle model (on/off)'
  printf, inp_lun, 'diffuse_frac    ', diffuse_frac_str,   '# Calculate diffuse fraction in images (on/off)'
  printf, inp_lun, 'derive_phys     ', derive_phys_str,    '# Calculate derived physical quantities (on/off)'
  printf, inp_lun, 'write_output    ', write_output_str,   '# Write results to output files (on/off)'
  printf, inp_lun, 'cleanup         ', cleanup_str,        '# [0] Keep auxiliary files to allow any of the above flags to be set to 0 for future runs; [1, DEFAULT] Delete auxiliary files after confirmation prompt; [2] WARNING: Automatically delete auxiliary files'
  printf, inp_lun, 'autoexit        ', autoexit_str,       '# Exit IDL automatically upon successful run completion (on/off)'

  printf, inp_lun, '# FLAGS 2 (use ancillary files)'
  printf, inp_lun, 'use_star2       ', use_star2_str,      '# Use different map for identifying stellar peaks and use the original for performing the flux calculation (on/off)'
  printf, inp_lun, 'use_gas2        ', use_gas2_str,       '# Use different map for identifying gas peaks and use the original for performing the flux calculation (on/off)'
  printf, inp_lun, 'use_star3       ', use_star3_str,      '# Use an additional stellar tracer map to mask pixels based on their flux ratios (on/off; e.g. Hα/FUV)'

  printf, inp_lun, '# FLAGS 3 (masking-related options)'
  printf, inp_lun, 'mstar_ext       ', mstar_ext_str,      '# Mask areas in starfile exterior to regions described in star_ext_mask (on/off)'
  printf, inp_lun, 'mstar_int       ', mstar_int_str,      '# Mask areas in starfile interior to regions described in star_int_mask (on/off)'
  printf, inp_lun, 'mgas_ext        ', mgas_ext_str,       '# Mask areas in gasfile exterior to regions described in gas_ext_mask (on/off)'
  printf, inp_lun, 'mgas_int        ', mgas_int_str,       '# Mask areas in gasfile interior to regions described in gas_int_mask (on/off)'
  printf, inp_lun, 'mstar_ext2      ', mstar_ext2_str,     '# Mask areas in starfile2 exterior to regions described in star_ext_mask2 (on/off)'
  printf, inp_lun, 'mstar_int2      ', mstar_int2_str,     '# Mask areas in starfile2 interior to regions described in star_int_mask2 (on/off)'
  printf, inp_lun, 'mgas_ext2       ', mgas_ext2_str,      '# Mask areas in gasfile2 exterior to regions described in gas_ext_mask2 (on/off)'
  printf, inp_lun, 'mgas_int2       ', mgas_int2_str,      '# Mask areas in gasfile2 interior to regions described in gas_int_mask2 (on/off)'
  printf, inp_lun, 'mstar_ext3      ', mstar_ext3_str,     '# Mask areas in starfile3 exterior to regions described in star_ext_mask3 (on/off)'
  printf, inp_lun, 'mstar_int3      ', mstar_int3_str,     '# Mask areas in starfile3 interior to regions described in star_int_mask3 (on/off)'
  printf, inp_lun, 'convert_masks   ', convert_masks_str,  '# Convert masks from other ds9 region coordinate systems to image coordinates (requires command line version of ds9) (on/off)'
  printf, inp_lun, 'cut_radius      ', cut_radius_str,     '# Mask maps outside a specified radial interval within the galaxy (on/off)'

  printf, inp_lun, '# FLAGS 4 (choose analysis options)'
  printf, inp_lun, 'set_centre      ', set_centre_str,     '# [0] Pixel coordinates of galaxy centre are set to the central pixel of starfile; [1, DEFAULT] Specify pixel coordinates of galaxy centre under INPUT PARAMETERS 1 (basic map analysis)'
  printf, inp_lun, 'tophat          ', tophat_str,         '# [0] Use Gaussian kernel to smoothen maps to larger aperture sizes; [1, DEFAULT] Use tophat kernel to smoothen maps to larger aperture sizes'
  printf, inp_lun, 'loglevels       ', loglevels_str,      '# [0] Linear equal-spacing of peak identification contour levels; [1, DEFAULT] Logarithmic equal-spacing of peak identification contour levels'
  printf, inp_lun, 'peak_find_tui   ', peak_find_tui_str,  '# [0, DEFAULT] select peaks only using INPUT PARAMETERS 3 (peak identification) without user interactivity; [1] use INPUT PARAMETERS 3 (peak identification) for initial peak finding and then enter interactive peak finding Text User Interface (TUI)'
  printf, inp_lun, 'flux_weight     ', flux_weight_str,    '# [0, DEFAULT] Peak positions correspond to brightest pixel in peak area; [1] Peak positions correspond to flux-weighted mean position of peak area'
  printf, inp_lun, 'calc_ap_area    ', calc_ap_area_str,   '# Calculate the area of apertures used based on number of unmasked pixels within the target aperture defined under INPUT PARAMETERS 2 (apertures) (on/off)'
  printf, inp_lun, 'tstar_incl      ', tstar_incl_str,     '# [0, DEFAULT] Reference time-scale tstariso does not include the overlap phase (tstar=tstariso+tover); [1] tstariso includes the overlap phase (tstar=tstariso)'
  printf, inp_lun, 'peak_prof       ', peak_prof_str,      '# [0] Model independent regions as points; [1] Model independent regions as constant-surface density discs; [2, DEFAULT] Model independent regions as two-dimensional Gaussians'
  printf, inp_lun, 'map_units       ', map_units_str,      '# [0] Unknown, do not calculate derived quantities; [1, DEFAULT] star=SFR and gas=gas; [2] star=gas1 and gas=gas2; [3] star=SFR1 and gas=SFR2'
  printf, inp_lun, 'star_tot_mode   ', star_tot_mode_str,  '# [0, DEFAULT] Use sfr/gas mass total for the stellar map calculated on the input map [1] use values of sfr/gas mass specified with the parameter star_tot_val'
  printf, inp_lun, 'gas_tot_mode    ', gas_tot_mode_str,   '# [0, DEFAULT] Use sfr/gas mass total for the gas map calculated on the input map [1] use values of sfr/gas mass specified with the parameter gas_tot_val'
  printf, inp_lun, 'use_X11         ', use_X11_str,        '# [0] Not allowed to create X11 windows; [1] Allowed to create X11 windows'
  printf, inp_lun, 'log10_output    ', log10_output_str,   '# [0] Values of output parameters in output file: galaxy_output.dat are written a double precision floating point numbers [1] Values of output parameters in output file: galaxy_output.dat are written with as log10(value)'






  printf, inp_lun, '# INPUT PARAMETERS 1 (basic map analysis)'
  printf, inp_lun, 'distance        ', distance_str,       '# Distance to galaxy in pc'
  printf, inp_lun, 'inclination     ', inclination_str,    '# Inclination angle in degrees'
  printf, inp_lun, 'posangle        ', posangle_str,       '# Position angle in degrees'
  printf, inp_lun, 'centrex         ', centrex_str,        '# Index of pixel x-axis coordinate of galaxy centre in starfile (important: start counting at zero on the left-hand side of the image)'
  printf, inp_lun, 'centrey         ', centrey_str,        '# Index of pixel y-axis coordinate of galaxy centre in starfile (important: start counting at zero on the bottom side of the image)'
  printf, inp_lun, 'minradius       ', minradius_str,      '# Minimum radius for analysis in pc'
  printf, inp_lun, 'maxradius       ', maxradius_str,      '# Maximum radius for analysis in pc'
  printf, inp_lun, 'Fs1_Fs2_min     ', Fs1_Fs2_min_str,    '# Minimum primary-to-secondary star formation tracer flux ratio for including pixels (only if use_star3=1)'
  printf, inp_lun, 'max_sample      ', max_sample_str,     '# Maximum number of pixels per map resolution FWHM'
  printf, inp_lun, 'astr_tolerance  ', astr_tolerance_str, '# allowable tolerance in decimal degrees between astrometric position of image pixels of different images (e.g. 0.9" = 2.5e-4 decimal degrees)'
  printf, inp_lun, 'nbins           ', nbins_str,          '# Number of bins used during sensitivity limit calculation to sample the intensity PDFs and fit Gaussians'

  printf, inp_lun, '# INPUT PARAMETERS 2 (apertures)'
  printf, inp_lun, 'lapmin          ', lapmin_str,         '# Minimum aperture size (i.e. diameter) in pc to create smoothened maps for'
  printf, inp_lun, 'lapmax          ', lapmax_str,         '# Maximum aperture size (i.e. diameter) in pc to create smoothened maps for'
  printf, inp_lun, 'naperture       ', naperture_str,      '# Number of aperture sizes'
  printf, inp_lun, 'peak_res        ', peak_res_str,       '# Minimum aperture size used in fitting the KL14 principle model, also index of aperture size at which peaks are identified (start counting at 0) - is set to map resolution if originally chosen to be smaller'
  printf, inp_lun, 'max_res         ', max_res_str,        '# Maximum aperture size used in fitting the KL14 principle model, also index of aperture size at which fluxes best reflect galactic averages (start counting at 0)'


  printf, inp_lun, '# INPUT PARAMETERS 3 (peak identification)'
  printf, inp_lun, 'npixmin         ', npixmin_str,        '# Minimum number of pixels for a valid peak (use npixmin=1 to allow points)'
  printf, inp_lun, 'nsigma          ', nsigma_str,         '# Multiple of the sensitivity limit needed for a valid peak'
  printf, inp_lun, 'logrange_s      ', logrange_s_str,     '# Logarithmic range covered by flux contour levels for stellar peak identification, counting down from stellar flux maximum'
  printf, inp_lun, 'logspacing_s    ', logspacing_s_str,   '# Logarithmic interval between flux contour levels for stellar peak identification'
  printf, inp_lun, 'logrange_g      ', logrange_g_str,     '# Logarithmic range covered by flux contour levels for gas peak identification, counting down from gas flux maximum'
  printf, inp_lun, 'logspacing_g    ', logspacing_g_str,   '# Logarithmic interval between flux contour levels for gas peak identification'
  printf, inp_lun, 'nlinlevel_s     ', nlinlevel_s_str,    '# Number of flux contour levels for stellar peak identification between zero and stellar flux maximum (only if loglevels=0)'
  printf, inp_lun, 'nlinlevel_g     ', nlinlevel_g_str,    '# Number of flux contour levels for gas peak identification between zero and gas flux maximum (only if loglevels=0)'



  printf, inp_lun, '# INPUT PARAMETERS 4 (timeline)'
  printf, inp_lun, 'tstariso        ', tstariso_str,            '# Reference time-scale spanned by star formation tracer in Myr'
  printf, inp_lun, 'tstariso_errmin ', tstariso_errmin_str,     '# Downward standard error of the reference time-scale in Myr'
  printf, inp_lun, 'tstariso_errmax ', tstariso_errmax_str,     '# Upward standard error of the reference time-scale in Myr'
  printf, inp_lun, 'tgasmini        ', tgasmini_str,            '# Minimum value of tgas considered during fitting process in Myr'
  printf, inp_lun, 'tgasmaxi        ', tgasmaxi_str,            '# Maximum value of tgas considered during fitting process in Myr'
  printf, inp_lun, 'tovermini       ', tovermini_str,           '# Minimum value of tover considered during fitting process in Myr'



  printf, inp_lun, '# INPUT PARAMETERS 5 (fitting)'
  printf, inp_lun, 'nmc             ', nmc_str,            '# Number of Monte-Carlo peak drawing experiments to be executed and averaged over during fitting process'
  printf, inp_lun, 'ndepth          ', ndepth_str,         '# Maximum number of free parameter array refinement loops for obtaining best-fitting value'
  printf, inp_lun, 'ntry            ', ntry_str,           '# Size of each free parameter array to obtain the best-fitting value'
  printf, inp_lun, 'nphysmc         ', nphysmc_str,        '# Number of Monte-Carlo experiments for error propagation of derived physical quantities'



  printf, inp_lun, '# INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)'
  printf, inp_lun, 'use_unfilt_ims  ', use_unfilt_ims_str, '# Calculate the diffuse fraction of [0] starfile and gasfile, the images used for the fitting process [1] star_unfilt_file and gas_unfilt_file, pass-through images (in the case of iterative diffuse filtering to calculate diffuse fraction of the unfiltered images)'
  printf, inp_lun, 'diffuse_quant   ', diffuse_quant_str,  '# Calculate the diffuse fraction of [0] flux, [1] power'
  printf, inp_lun, 'f_filter_type   ', f_filter_type_str,  '# Filter type choice for Fourier filtering: [0] butterworth filter [1] gaussian filter [2] ideal filter (Only relevant if diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
  printf, inp_lun, 'bw_order        ', bw_order_str,       '# The order of the butterworth filter (Must be an integer. Only relevant if f_filter_type = 0 [butterworth] and diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
  printf, inp_lun, 'filter_len_conv ', filter_len_conv_str,'# conversion factor for Fourier filtering cut_length (cut_length = lambda * filter_len_conv)'
  printf, inp_lun, 'emfrac_cor_mode ', emfrac_cor_mode_str,'# Apply corrections to measured fgmc and fcl [0] apply no correction [1] apply correction for flux loss from signal regions [2] apply correction for the filling factor, zeta, only [3] apply both corrections (flux loss from signal regions and filling factor, zeta)'
  printf, inp_lun, 'rpeak_cor_mode  ', rpeak_cor_mode_str,   '# [0] use measured r_peaks [1] use supplied values of rpeaks_cor_val and rpeakg_cor_val'
  printf, inp_lun, 'rpeaks_cor_val  ', rpeaks_cor_val_str,   '# value of r_peak_star to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'
  printf, inp_lun, 'rpeaks_cor_emin ', rpeaks_cor_emin_str, '# downwards error of r_peak_star to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'
  printf, inp_lun, 'rpeaks_cor_emax ', rpeaks_cor_emax_str, '# upwards error of r_peak_star to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'
  printf, inp_lun, 'rpeakg_cor_val  ', rpeakg_cor_val_str,   '# value of r_peak_gas to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'
  printf, inp_lun, 'rpeakg_cor_emin ', rpeakg_cor_emin_str, '# downwards error of r_peak_gas to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'
  printf, inp_lun, 'rpeakg_cor_emax ', rpeakg_cor_emax_str, '# upwards error of r_peak_gas to use when correcting for flux_loss. Only used if emfrac_cor_mode = 1 or 3 and rpeak_cor_mode = 2'

  printf, inp_lun, '# INPUT PARAMETERS 7 (conversions and constants to calculate derived quantities)'
  printf, inp_lun, 'convstar        ', convstar_str,       '# Log of multiplication factor to convert the pixel value in starfile to a SFR in Msun/yr (if maptypes=1 or maptypes=3) or to a gas mass in Msun (if maptypes=2)'
  printf, inp_lun, 'convstar_rerr   ', convstar_rerr_str,  '# Relative error (sigma_x/x) of convstar'
  printf, inp_lun, 'convgas         ', convgas_str,        '# Log of multiplication factor to convert the pixel value in gasfile to a gas mass in Msun (if maptypes=1 or maptypes=2) or to a SFR in Msun/yr (if maptypes=3)'
  printf, inp_lun, 'convgas_rerr    ', convgas_rerr_str,   '# Relative error (sigma_x/x) of convgas'
  printf, inp_lun, 'convstar3       ', convstar3_str,      '# Log of multiplication factor to convert the pixel value in starfile3 to a SFR in Msun/yr'
  printf, inp_lun, 'convstar3_rerr  ', convstar3_rerr_str, '# Relative error (sigma_x/x) of convstar3'
  printf, inp_lun, 'lighttomass     ', lighttomass_str,    '# Light-to-mass ratio of a desired feedback mechanism in m^2 s^-3 (default 0.002 is SNe)'
  printf, inp_lun, 'momratetomass   ', momratetomass_str,  '# Momentum output rate per unit mass of a desired feedback mechanism in m s^-2 (default 5.e-10 is SNe+winds)'
  printf, inp_lun, 'star_tot_val    ', star_tot_val_str,   '# input value of the total sfr/gas mass of the star map (only used if star_tot_mode = 1)'
  printf, inp_lun, 'star_tot_err    ', star_tot_err_str,   '# Error on the input value of the total sfr/gas mass of the star map (only used if star_tot_mode = 1)'
  printf, inp_lun, 'gas_tot_val     ', gas_tot_val_str,    '# input value of the total sfr/gas mass of the gas map (only used if gas_tot_mode = 1)'
  printf, inp_lun, 'gas_tot_err     ', gas_tot_err_str,    '# Error on the input value of the total sfr/gas mass of the gas map (only used if star_tot_mode = 1)'

  printf, inp_lun, '# INPUT PARAMETERS 8 (sensitivity)'
  printf, inp_lun, 'use_stds        ', use_stds_str,       '# [0] calculate standard deviations of images [1] Use supplied standard deviation measurements'
  printf, inp_lun, 'std_star        ', std_star_str,       '# The measured standard deviation of starfile'
  printf, inp_lun, 'std_star3       ', std_star3_str,      '# The measured standard deviation of starfile3 (only used of use_star3 = 1)'
  printf, inp_lun, 'std_gas         ', std_gas_str,        '# The measured standard deviation of gasfile'


  printf, inp_lun, '# INPUT PARAMETERS 9 (noise threshold)'
  printf, inp_lun, 'use_noisecut    ', use_noisecut_str,       '# [0] calculate standard deviations of images [1] Use supplied standard deviation measurements'
  printf, inp_lun, 'noisethresh_s   ', noisethresh_s_str,      '# The measured standard deviation of starfile'
  printf, inp_lun, 'noisethresh_g   ', noisethresh_g_str,      '# The measured standard deviation of starfile3 (only used of use_star3 = 1)'

  if iter_input_switch eq 1 then begin ; create input file for iteration
    printf, inp_lun, '# INPUT PARAMETERS 10 (Fourier diffuse removal iteration)'
    printf, inp_lun, 'use_guess       ', use_guess_str,      '#  [0] Do not Filter images before first KL14 run and ignore parameter: "initial_guess" [1] Filter images with initial guess of lambda before first KL14 run'
    printf, inp_lun, 'initial_guess   ', initial_guess_str,  '#  Length in pc of the initial estimate of lambda with which to filter the images'
    printf, inp_lun, 'iter_criterion  ', iter_criterion_str, '#  Fractional difference between lambda and previous lambda value(s) that triggers an end to the iteration process.'
    printf, inp_lun, 'iter_crit_len   ', iter_crit_len_str,  '#  Number of previous runs of KL14 over which iter_criterion must be true, i.e. if iter_crit_len = 1 then the current value of lambda is compared to the previous value of lambda, if iter_crit_len = 2 then then the current value of lambda is compared to the previous value of lambda and the value of lambda prior to the previous value'
    printf, inp_lun, 'iter_nmax       ', iter_nmax_str,      '#  Maximum number of KL14 runs allowed in the iteration process. If iter_criterion is not reached before this number of runs the programme will exit anyway.'
    printf, inp_lun, 'iter_filter     ', iter_filter_str,    '#  Filter type choice for iterative Fourier filtering: [0] butterworth filter [1] gaussian filter [2] ideal filter (Only relevant if diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
    printf, inp_lun, 'iter_bwo        ', iter_bwo_str,       '#  The order of the butterworth filter for iterative Fourier filtering (Must be an integer. Only relevant if iter_bwo = 0 [butterworth], but a dummy value should be supplied anyway)'
    printf, inp_lun, 'iter_len_conv   ', iter_len_conv_str,  '#  conversion factor for iterative Fourier filtering cut_length (cut_length = lambda * filter_len_conv)'
    printf, inp_lun, 'iter_rpeak_mode ', iter_rpeak_mode_str,'# [0] use mode specified by rpeak_cor_mode [1] use rpeak values from iter0'
    printf, inp_lun, 'iter_tot_mode_s', iter_tot_mode_s_str, '# [0, DEFAULT] Use the sfr/gas mass total for the stellar map calculated on the input map of each iteration [1] Use the sfr/gas mass total for the stellar map calculated on the initial iteration for all subsequent iterations [2] Use the sfr/gas mass total for the stellar map specified by the parameter star_tot_val for all iterations'
    printf, inp_lun, 'iter_tot_mode_g', iter_tot_mode_g_str, '# [0, DEFAULT] Use the sfr/gas mass total for the gas map calculated on the input map of each iteration [1] Use the sfr/gas mass total for the gas map calculated on the initial iteration for all subsequent iterations [2] Use the sfr/gas mass total for the gas map specified by the parameter star_tot_val for all iterations'
    printf, inp_lun, 'iter_autoexit   ', iter_autoexit_str,  '# Exit IDL automatically upon successful iteration run completion (on/off)'
    printf, inp_lun, 'use_nice        ', use_nice_str,       '# (on/off) Use the "nice" command when spawning KL14 runs to the terminal to adjust the process priority. Note that this may not work depending on the system configuration'
    printf, inp_lun, 'nice_value      ', nice_value_str,     '# Value of the priority adjustment must be an integer between -20 (the highest) to 19 (the lowest). (Note this varies according to the implementation of nice)'


  endif

  printf, inp_lun, '################################################'
  printf, inp_lun, '#                                              #'
  printf, inp_lun, '#              END PARAMETER FILE              #'
  printf, inp_lun, '#                                              #'
  printf, inp_lun, '################################################'



  free_lun, inp_lun


end
