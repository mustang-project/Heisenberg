pro make_input_file, input_file_filepath   $ ; general variables
                   , datadir, galaxy, starfile, gasfile   $  ; File names compulsory variables

                   , distance, inclination, posangle, centrex, centrey, minradius, maxradius, Fs1_Fs2_min, max_sample, astr_tolerance, nbins $ ; # INPUT PARAMETERS 1 (basic map analysis)
                   , lapmin, lapmax, naperture, peak_res, max_res $  ; # INPUT PARAMETERS 2 (apertures)
                  ; ################################################################################################################################################################################################################################################################

                   , starfile2 = starfile2, gasfile2 = gasfile2, starfile3 = starfile3  $  ; File names keywords
                   , maskdir = maskdir, star_ext_mask1 = star_ext_mask, star_int_mask1 = star_int_mask, gas_ext_mask1 = gas_ext_mask, gas_int_mask1 = gas_int_mask, star_ext_mask2 = star_ext_mask2, star_int_mask2 = star_int_mask2, gas_ext_mask2 = gas_ext_mask2, gas_int_mask2 = gas_int_mask2, star_ext_mask3 = star_ext_mask3, star_int_mask3 = star_int_mask3  $  ; Mask file names keywords ; note star_ext_mask1 = star_ext_mask etc. prevents ambigious keyword error
                   , mask_images = mask_images, regrid = regrid, smoothen = smoothen, sensitivity = sensitivity, id_peaks = id_peaks, calc_ap_flux = calc_ap_flux, generate_plot = generate_plot, get_distances = get_distances, calc_obs = calc_obs, calc_fit = calc_fit, diffuse_frac = diffuse_frac, derive_phys = derive_phys, write_output = write_output, cleanup = cleanup, autoexit = autoexit $ ; FLAGS 1 (keywords)'
                   , use_star2 = use_star2, use_gas2 = use_gas2, use_star3 = use_star3 $ ;  FLAGS 2 keywords
                   , mstar_ext1 = mstar_ext, mstar_int1 = mstar_int, mgas_ext1 = mgas_ext, mgas_int1 = mgas_int, mstar_ext2 = mstar_ext2, mstar_int2 = mstar_int2, mgas_ext2 = mgas_ext2, mgas_int2 = mgas_int2, mstar_ext3 = mstar_ext3, mstar_int3 = mstar_int3, convert_masks = convert_masks, cut_radius = cut_radius $ ; # FLAGS 3 (masking-related options) keywords ; note mstar_ext1 = mstar_ext etc. prevents ambigious keyword error
                   , set_centre = set_centre, tophat = tophat, loglevels = loglevels, flux_weight = flux_weight, calc_ap_area = calc_ap_area, tstar_incl = tstar_incl, peak_prof = peak_prof, map_units = map_units, use_X11 = use_X11 $ ; # FLAGS 4 (choose analysis options)

                   , npixmin = npixmin, nsigma = nsigma, logrange_s = logrange_s, logspacing_s = logspacing_s, logrange_g = logrange_g, logspacing_g = logspacing_g, nlinlevel_s = nlinlevel_s, nlinlevel_g = nlinlevel_g $ ; # INPUT PARAMETERS 3 (peak identification)
                   , tstariso_val = tstariso, tstariso_errmin = tstariso_errmin, tstariso_errmax = tstariso_errmax, tgasmini = tgasmini, tgasmaxi = tgasmaxi, tovermini = tovermini $ ; # INPUT PARAMETERS 4 (timeline) ; note tstariso_val = tstariso prevents ambigious keyword error
                   , nmc = nmc, ndepth = ndepth, ntry = ntry, nphysmc = nphysmc $ ; # INPUT PARAMETERS 5 (fitting)
                   , diffuse_quant, f_filter_type, bw_order, filter_len_conv $ ; # INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)
                   , convstar_val = convstar, convstar_rerr = convstar_rerr, convgas_val = convgas, convgas_rerr = convgas_rerr, convstar3_val = convstar3, convstar3_rerr = convstar3_rerr, lighttomass = lighttomass, momratetomass = momratetomass $; # INPUT PARAMETERS 6 (conversions and constants to calculate derived quantities) ; note convgas_val = convgas avoids % Ambiguous keyword abbreviation: CONVGAS.
                   , use_stds = use_stds, std_star1 = std_star, std_star3 = std_star3, std_gas = std_gas $ ; # INPUT PARAMETERS 8 (sensitivity) ; star_std1 prevents ambigious keywords

                   ; ################################################################################################################################################################################################################################################################
                   , use_guess = use_guess, initial_guess = initial_guess, iter_criterion = iter_criterion, iter_crit_len = iter_crit_len, iter_nmax = iter_nmax, iter_filter = iter_filter, iter_bwo = iter_bwo, iter_len_conv = iter_len_conv ; # INPUT PARAMETERS 9 (Fourier diffuse removal iteration) # note that these parameters are ignored if the call sequence idl kl14 -arg [full/absolute path of input file] is used. They are only used if the call sequence idl iterate_kl14 -arg [full/absolute path of input file] is used.



  openw, inp_lun, input_file_filepath, width = 350, /get_lun


  pad_string = ' '
  comment_sep_length = 9



  ; #############################################
  ; # FILE NAMES
  ; #############################################
  if n_elements(starfile2) ne 1 then starfile2 = '-'
  if n_elements(gasfile2) ne 1 then gasfile2 = '-'
  if n_elements(starfile3) ne 1 then starfile3 = '-'

  datadir = string(datadir)
  galaxy = string(galaxy)
  starfile = string(starfile)
  starfile2 = string(starfile2)
  gasfile = string(gasfile)
  gasfile2 = string(gasfile2)
  starfile3 = string(starfile3)

  max_string_len = max([strlen(datadir), strlen(galaxy), strlen(starfile), strlen(gasfile), strlen(starfile2), strlen(gasfile2), strlen(starfile3)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(datadir) lt max_string_len) do datadir = datadir + pad_string
  while(strlen(galaxy) lt max_string_len) do galaxy = galaxy + pad_string
  while(strlen(starfile) lt max_string_len) do starfile = starfile + pad_string
  while(strlen(starfile2) lt max_string_len) do starfile2 = starfile2 + pad_string
  while(strlen(gasfile) lt max_string_len) do gasfile = gasfile + pad_string
  while(strlen(gasfile2) lt max_string_len) do gasfile2 = gasfile2 + pad_string
  while(strlen(starfile3) lt max_string_len) do starfile3 = starfile3 + pad_string





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

  maskdir = string(maskdir)
  star_ext_mask = string(star_ext_mask)
  star_int_mask = string(star_int_mask)
  gas_ext_mask = string(gas_ext_mask)
  gas_int_mask = string(gas_int_mask)
  star_ext_mask2 = string(star_ext_mask2)
  star_int_mask2 = string(star_int_mask2)
  gas_ext_mask2 = string(gas_ext_mask2)
  gas_int_mask2 = string(gas_int_mask2)
  star_ext_mask3 = string(star_ext_mask3)
  star_int_mask3 = string(star_int_mask3)

  max_string_len = max([strlen(maskdir), strlen(star_ext_mask), strlen(star_int_mask), strlen(gas_ext_mask), strlen(gas_int_mask), strlen(star_ext_mask2), strlen(star_int_mask2), strlen(gas_ext_mask2), strlen(gas_int_mask2), strlen(star_ext_mask3), strlen(star_int_mask3)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while( strlen(maskdir) lt max_string_len) do maskdir = maskdir + pad_string
  while( strlen(star_ext_mask) lt max_string_len) do star_ext_mask = star_ext_mask + pad_string
  while( strlen(star_int_mask) lt max_string_len) do star_int_mask = star_int_mask + pad_string
  while( strlen(gas_ext_mask) lt max_string_len) do gas_ext_mask = gas_ext_mask + pad_string
  while( strlen(gas_int_mask) lt max_string_len) do gas_int_mask = gas_int_mask + pad_string
  while( strlen(star_ext_mask2) lt max_string_len) do star_ext_mask2 = star_ext_mask2 + pad_string
  while( strlen(star_int_mask2) lt max_string_len) do star_int_mask2 = star_int_mask2 + pad_string
  while( strlen(gas_ext_mask2) lt max_string_len) do gas_ext_mask2 = gas_ext_mask2 + pad_string
  while( strlen(gas_int_mask2) lt max_string_len) do gas_int_mask2 = gas_int_mask2 + pad_string
  while( strlen(star_ext_mask3) lt max_string_len) do star_ext_mask3 = star_ext_mask3 + pad_string
  while( strlen(star_int_mask3) lt max_string_len) do star_int_mask3 = star_int_mask3 + pad_string



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

  mask_images = strcompress(string(mask_images))
  regrid = strcompress(string(regrid))
  smoothen = strcompress(string(smoothen))
  sensitivity = strcompress(string(sensitivity))
  id_peaks = strcompress(string(id_peaks))
  calc_ap_flux = strcompress(string(calc_ap_flux))
  generate_plot = strcompress(string(generate_plot))
  get_distances = strcompress(string(get_distances))
  calc_obs = strcompress(string(calc_obs))
  calc_fit = strcompress(string(calc_fit))
  diffuse_frac = strcompress(string(diffuse_frac))
  derive_phys = strcompress(string(derive_phys))
  write_output = strcompress(string(write_output))
  cleanup = strcompress(string(cleanup))
  autoexit = strcompress(string(autoexit))


  max_string_len = max([strlen(mask_images), strlen(regrid), strlen(smoothen), strlen(sensitivity), strlen(id_peaks), strlen(calc_ap_flux), strlen(generate_plot), strlen(get_distances), strlen(calc_obs), strlen(calc_fit), strlen(derive_phys), strlen(write_output), strlen(cleanup), strlen(autoexit)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while (strlen(mask_images) lt max_string_len) do mask_images = mask_images + pad_string
  while (strlen(regrid) lt max_string_len) do regrid = regrid + pad_string
  while (strlen(smoothen) lt max_string_len) do smoothen = smoothen + pad_string
  while (strlen(sensitivity) lt max_string_len) do sensitivity = sensitivity + pad_string
  while (strlen(id_peaks) lt max_string_len) do id_peaks = id_peaks + pad_string
  while (strlen(calc_ap_flux) lt max_string_len) do calc_ap_flux = calc_ap_flux + pad_string
  while (strlen(generate_plot) lt max_string_len) do generate_plot = generate_plot + pad_string
  while (strlen(get_distances) lt max_string_len) do get_distances = get_distances + pad_string
  while (strlen(calc_obs) lt max_string_len) do calc_obs = calc_obs + pad_string
  while (strlen(calc_fit) lt max_string_len) do calc_fit = calc_fit + pad_string
  while (strlen(diffuse_frac) lt max_string_len) do diffuse_frac = diffuse_frac + pad_string
  while (strlen(derive_phys) lt max_string_len) do derive_phys = derive_phys + pad_string
  while (strlen(write_output) lt max_string_len) do write_output = write_output + pad_string
  while (strlen(cleanup) lt max_string_len) do cleanup = cleanup + pad_string
  while (strlen(autoexit) lt max_string_len) do autoexit = autoexit + pad_string


  ; #############################################
  ;'# FLAGS 2 (use ancillary files)''
  ; #############################################

  if n_elements(use_star2) ne 1 then use_star2 = 0
  if n_elements(use_gas2) ne 1 then use_gas2 = 0
  if n_elements(use_star3) ne 1 then use_star3 = 0

  use_star2 = strcompress(string(use_star2))
  use_gas2 = strcompress(string(use_gas2))
  use_star3 = strcompress(string(use_star3))

  max_string_len = max([strlen(use_star2), strlen(use_gas2), strlen(use_star3)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(use_star2) lt max_string_len) do use_star2 = use_star2 + pad_string
  while(strlen(use_gas2) lt max_string_len) do use_gas2 = use_gas2 + pad_string
  while(strlen(use_star3) lt max_string_len) do use_star3 = use_star3 + pad_string




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

  mstar_ext = strcompress(string(mstar_ext))
  mstar_int = strcompress(string(mstar_int))
  mgas_ext = strcompress(string(mgas_ext))
  mgas_int = strcompress(string(mgas_int))
  mstar_ext2 = strcompress(string(mstar_ext2))
  mstar_int2 = strcompress(string(mstar_int2))
  mgas_ext2 = strcompress(string(mgas_ext2))
  mgas_int2 = strcompress(string(mgas_int2))
  mstar_ext3 = strcompress(string(mstar_ext3))
  mstar_int3 = strcompress(string(mstar_int3))
  convert_masks = strcompress(string(convert_masks))
  cut_radius = strcompress(string(cut_radius))

  max_string_len = max([strlen(mstar_ext), strlen(mstar_int), strlen(mgas_ext), strlen(mgas_int), strlen(mstar_ext2), strlen(mstar_int2), strlen(mgas_ext2), strlen(mgas_int2), strlen(mstar_ext3), strlen(mstar_int3), strlen(convert_masks), strlen(cut_radius)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(mstar_ext) lt max_string_len) do mstar_ext = mstar_ext + pad_string
  while(strlen(mstar_int) lt max_string_len) do mstar_int = mstar_int + pad_string
  while(strlen(mgas_ext) lt max_string_len) do mgas_ext = mgas_ext + pad_string
  while(strlen(mgas_int) lt max_string_len) do mgas_int = mgas_int + pad_string
  while(strlen(mstar_ext2) lt max_string_len) do mstar_ext2 = mstar_ext2 + pad_string
  while(strlen(mstar_int2) lt max_string_len) do mstar_int2 = mstar_int2 + pad_string
  while(strlen(mgas_ext2) lt max_string_len) do mgas_ext2 = mgas_ext2 + pad_string
  while(strlen(mgas_int2) lt max_string_len) do mgas_int2 = mgas_int2 + pad_string
  while(strlen(mstar_ext3) lt max_string_len) do mstar_ext3 = mstar_ext3 + pad_string
  while(strlen(mstar_int3) lt max_string_len) do mstar_int3 = mstar_int3 + pad_string
  while(strlen(convert_masks) lt max_string_len) do convert_masks = convert_masks + pad_string
  while(strlen(cut_radius) lt max_string_len) do cut_radius = cut_radius + pad_string




  ; #############################################
  ; # FLAGS 4 (choose analysis options)
  ; #############################################

  if n_elements(set_centre) ne 1 then set_centre = 1
  if n_elements(tophat) ne 1 then tophat = 1
  if n_elements(loglevels) ne 1 then loglevels = 1
  if n_elements(flux_weight) ne 1 then flux_weight = 0
  if n_elements(calc_ap_area) ne 1 then calc_ap_area = 1
  if n_elements(tstar_incl) ne 1 then tstar_incl = 0
  if n_elements(peak_prof) ne 1 then peak_prof = 2
  if n_elements(map_units) ne 1 then map_units = 1
  if n_elements(use_X11) ne 1 then use_X11 = 1

  set_centre = strcompress(string(set_centre))
  tophat = strcompress(string(tophat))
  loglevels = strcompress(string(loglevels))
  flux_weight = strcompress(string(flux_weight))
  calc_ap_area = strcompress(string(calc_ap_area))
  tstar_incl = strcompress(string(tstar_incl))
  peak_prof = strcompress(string(peak_prof))
  map_units = strcompress(string(map_units))
  use_X11 = strcompress(string(use_X11))

  max_string_len = max([strlen(set_centre), strlen(tophat), strlen(loglevels), strlen(flux_weight), strlen(calc_ap_area), strlen(tstar_incl), strlen(peak_prof), strlen(map_units), strlen(use_X11)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(set_centre) lt max_string_len) do set_centre = set_centre + pad_string
  while(strlen(tophat) lt max_string_len) do tophat = tophat + pad_string
  while(strlen(loglevels) lt max_string_len) do loglevels = loglevels + pad_string
  while(strlen(flux_weight) lt max_string_len) do flux_weight = flux_weight + pad_string
  while(strlen(calc_ap_area) lt max_string_len) do calc_ap_area = calc_ap_area + pad_string
  while(strlen(tstar_incl) lt max_string_len) do tstar_incl = tstar_incl + pad_string
  while(strlen(peak_prof) lt max_string_len) do peak_prof = peak_prof + pad_string
  while(strlen(map_units) lt max_string_len) do map_units = map_units + pad_string
  while(strlen(use_X11) lt max_string_len) do use_X11 = use_X11 + pad_string










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


  distance = strcompress(string(distance))
  inclination = strcompress(string(inclination))
  posangle = strcompress(string(posangle))
  centrex = strcompress(string(centrex))
  centrey = strcompress(string(centrey))
  minradius = strcompress(string(minradius))
  maxradius = strcompress(string(maxradius))
  Fs1_Fs2_min = strcompress(string(Fs1_Fs2_min))
  max_sample = strcompress(string(max_sample))
  astr_tolerance = strcompress(string(astr_tolerance))
  nbins = strcompress(string(nbins))

  max_string_len = max([strlen(distance), strlen(inclination), strlen(posangle), strlen(centrex), strlen(centrey), strlen(minradius), strlen(maxradius), strlen(Fs1_Fs2_min), strlen(max_sample), strlen(astr_tolerance), strlen(nbins)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(distance) lt max_string_len) do distance = distance + pad_string
  while(strlen(inclination) lt max_string_len) do inclination = inclination + pad_string
  while(strlen(posangle) lt max_string_len) do posangle = posangle + pad_string
  while(strlen(centrex) lt max_string_len) do centrex = centrex + pad_string
  while(strlen(centrey) lt max_string_len) do centrey = centrey + pad_string
  while(strlen(minradius) lt max_string_len) do minradius = minradius + pad_string
  while(strlen(maxradius) lt max_string_len) do maxradius = maxradius + pad_string
  while(strlen(Fs1_Fs2_min) lt max_string_len) do Fs1_Fs2_min = Fs1_Fs2_min + pad_string
  while(strlen(max_sample) lt max_string_len) do max_sample = max_sample + pad_string
  while(strlen(astr_tolerance) lt max_string_len) do astr_tolerance = astr_tolerance + pad_string
  while(strlen(nbins) lt max_string_len) do nbins = nbins + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 2 (apertures)
  ; #############################################

  if n_elements(lapmin) ne 1 then stop
  if n_elements(lapmax) ne 1 then stop
  if n_elements(naperture) ne 1 then stop
  if n_elements(peak_res) ne 1 then stop
  if n_elements(max_res) ne 1 then stop

  lapmin = strcompress(string(lapmin))
  lapmax = strcompress(string(lapmax))
  naperture = strcompress(string(naperture))
  peak_res = strcompress(string(peak_res))
  max_res = strcompress(string(max_res))

  max_string_len = max([strlen(lapmin), strlen(lapmax), strlen(naperture), strlen(peak_res), strlen(max_res)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(lapmin) lt max_string_len) do lapmin = lapmin + pad_string
  while(strlen(lapmax) lt max_string_len) do lapmax = lapmax + pad_string
  while(strlen(naperture) lt max_string_len) do naperture = naperture + pad_string
  while(strlen(peak_res) lt max_string_len) do peak_res = peak_res + pad_string
  while(strlen(max_res) lt max_string_len) do max_res = max_res + pad_string

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

  npixmin = strcompress(string(npixmin))
  nsigma = strcompress(string(nsigma))
  logrange_s = strcompress(string(logrange_s))
  logspacing_s = strcompress(string(logspacing_s))
  logrange_g = strcompress(string(logrange_g))
  logspacing_g = strcompress(string(logspacing_g))
  nlinlevel_s = strcompress(string(nlinlevel_s))
  nlinlevel_g = strcompress(string(nlinlevel_g))

  max_string_len = max([strlen(npixmin), strlen(nsigma), strlen(logrange_s), strlen(logspacing_s), strlen(logrange_g), strlen(logspacing_g), strlen(nlinlevel_s), strlen(nlinlevel_g)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(npixmin) lt max_string_len) do npixmin = npixmin + pad_string
  while(strlen(nsigma) lt max_string_len) do nsigma = nsigma + pad_string
  while(strlen(logrange_s) lt max_string_len) do logrange_s = logrange_s + pad_string
  while(strlen(logspacing_s) lt max_string_len) do logspacing_s = logspacing_s + pad_string
  while(strlen(logrange_g) lt max_string_len) do logrange_g = logrange_g + pad_string
  while(strlen(logspacing_g) lt max_string_len) do logspacing_g = logspacing_g + pad_string
  while(strlen(nlinlevel_s) lt max_string_len) do nlinlevel_s = nlinlevel_s + pad_string
  while(strlen(nlinlevel_g) lt max_string_len) do nlinlevel_g = nlinlevel_g + pad_string

  ; #############################################
  ; # INPUT PARAMETERS 4 (timeline)
  ; #############################################

  if  n_elements(tstariso) ne 1 then tstariso = 1.
  if  n_elements(tstariso_errmin) ne 1 then tstariso_errmin = 0.
  if  n_elements(tstariso_errmax) ne 1 then tstariso_errmax = 0.
  if  n_elements(tgasmini) ne 1 then tgasmini = 0.1
  if  n_elements(tgasmaxi) ne 1 then tgasmaxi = 1000.
  if  n_elements(tovermini) ne 1 then tovermini = 0.01

  tstariso = strcompress(string(tstariso))
  tstariso_errmin = strcompress(string(tstariso_errmin))
  tstariso_errmax = strcompress(string(tstariso_errmax))
  tgasmini = strcompress(string(tgasmini))
  tgasmaxi = strcompress(string(tgasmaxi))
  tovermini = strcompress(string(tovermini))

  max_string_len = max([strlen(tstariso), strlen(tstariso_errmin), strlen(tstariso_errmax), strlen(tgasmini), strlen(tgasmaxi), strlen(tovermini)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(tstariso) lt max_string_len) do tstariso = tstariso + pad_string
  while(strlen(tstariso_errmin) lt max_string_len) do tstariso_errmin = tstariso_errmin + pad_string
  while(strlen(tstariso_errmax) lt max_string_len) do tstariso_errmax = tstariso_errmax + pad_string
  while(strlen(tgasmini) lt max_string_len) do tgasmini = tgasmini + pad_string
  while(strlen(tgasmaxi) lt max_string_len) do tgasmaxi = tgasmaxi + pad_string
  while(strlen(tovermini) lt max_string_len) do tovermini = tovermini + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 5 (fitting)
  ; #############################################
  if n_elements(nmc) ne 1 then nmc = 1000
  if n_elements(ndepth) ne 1 then ndepth = 4
  if n_elements(ntry) ne 1 then ntry = 101
  if n_elements(nphysmc) ne 1 then nphysmc = 1000000

  nmc = strcompress(string(nmc))
  ndepth = strcompress(string(ndepth))
  ntry = strcompress(string(ntry))
  nphysmc = strcompress(string(nphysmc))

  max_string_len = max( [strlen(nmc), strlen(ndepth), strlen(ntry), strlen(nphysmc)]) + comment_sep_length
  max_string_len = max([max_string_len,22])


  while(strlen(nmc) lt max_string_len) do nmc = nmc + pad_string
  while(strlen(ndepth) lt max_string_len) do ndepth = ndepth + pad_string
  while(strlen(ntry) lt max_string_len) do ntry = ntry + pad_string
  while(strlen(nphysmc) lt max_string_len) do nphysmc = nphysmc + pad_string


  ; #############################################
  ; # INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)
  ; #############################################
  if n_elements(diffuse_quant) ne 1 then diffuse_quant = 1
  if n_elements(f_filter_type) ne 1 then f_filter_type = 2
  if n_elements(bw_order) ne 1 then bw_order = 2
  if n_elements(filter_len_conv) ne 1 then filter_len_conv = 1.0

  diffuse_quant = strcompress(string(diffuse_quant))
  f_filter_type = strcompress(string(f_filter_type))
  bw_order = strcompress(string(bw_order))
  filter_len_conv = strcompress(string(filter_len_conv))

  max_string_len = max([strlen(diffuse_quant), strlen(f_filter_type), strlen(bw_order), strlen(filter_len_conv)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(diffuse_quant) lt max_string_len) do diffuse_quant = diffuse_quant + pad_string
  while(strlen(f_filter_type) lt max_string_len) do f_filter_type = f_filter_type + pad_string
  while(strlen(bw_order) lt max_string_len) do bw_order = bw_order + pad_string
  while(strlen(filter_len_conv) lt max_string_len) do filter_len_conv = filter_len_conv + pad_string

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

  convstar = strcompress(string(convstar))
  convstar_rerr = strcompress(string(convstar_rerr))
  convgas = strcompress(string(convgas))
  convgas_rerr = strcompress(string(convgas_rerr))
  convstar3 = strcompress(string(convstar3))
  convstar3_rerr = strcompress(string(convstar3_rerr))
  lighttomass = strcompress(string(lighttomass))
  momratetomass = strcompress(string(momratetomass))

  max_string_len = max([strlen(convstar), strlen(convstar_rerr), strlen(convgas), strlen(convgas_rerr), strlen(convstar3), strlen(convstar3_rerr), strlen(lighttomass), strlen(momratetomass)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(convstar) lt max_string_len) do convstar = convstar + pad_string
  while(strlen(convstar_rerr) lt max_string_len) do convstar_rerr = convstar_rerr + pad_string
  while(strlen(convgas) lt max_string_len) do convgas = convgas + pad_string
  while(strlen(convgas_rerr) lt max_string_len) do convgas_rerr = convgas_rerr + pad_string
  while(strlen(convstar3) lt max_string_len) do convstar3 = convstar3 + pad_string
  while(strlen(convstar3_rerr) lt max_string_len) do convstar3_rerr = convstar3_rerr + pad_string
  while(strlen(lighttomass) lt max_string_len) do lighttomass = lighttomass + pad_string
  while(strlen(momratetomass) lt max_string_len) do momratetomass = momratetomass + pad_string



  ; #############################################
  ; # INPUT PARAMETERS 8 (sensitivity)
  ; #############################################


  if n_elements(use_stds) ne 1 then use_stds =  0
  if n_elements(std_star) ne 1 then std_star =  0.1
  if n_elements(std_star3) ne 1 then std_star3 =  0.1
  if n_elements(std_gas) ne 1 then std_gas =  0.1

  use_stds = strcompress(string(use_stds))
  std_star = strcompress(string(std_star))
  std_star3 = strcompress(string(std_star3))
  std_gas = strcompress(string(std_gas))

  max_string_len = max([strlen(use_stds), strlen(std_star), strlen(std_star3), strlen(std_gas)]) + comment_sep_length
  max_string_len = max([max_string_len,22])

  while(strlen(use_stds) lt max_string_len) do use_stds = use_stds + pad_string
  while(strlen(std_star) lt max_string_len) do std_star = std_star + pad_string
  while(strlen(std_star3) lt max_string_len) do std_star3 = std_star3 + pad_string
  while(strlen(std_gas) lt max_string_len) do std_gas = std_gas + pad_string





  ; ##########################################################################################
  ; # INPUT PARAMETERS 9 (Fourier diffuse removal iteration) -- Optional
  ; ##########################################################################################
  iter_input_switch = 0
  if n_elements(use_guess) eq 1 && n_elements(initial_guess) eq 1 && n_elements(iter_criterion) eq 1 && n_elements(iter_crit_len) eq 1 && n_elements(iter_nmax) eq 1 && n_elements(iter_filter) eq 1 && n_elements(iter_bwo) eq 1 && n_elements(iter_len_conv) eq 1 then begin
    iter_input_switch = 1 ; use so don't need to repeat test


    use_guess = strcompress(string(use_guess))
    initial_guess = strcompress(string(initial_guess))
    iter_criterion = strcompress(string(iter_criterion))
    iter_crit_len = strcompress(string(iter_crit_len))
    iter_nmax = strcompress(string(iter_nmax))
    iter_filter = strcompress(string(iter_filter))
    iter_bwo = strcompress(string(iter_bwo))
    iter_len_conv = strcompress(string(iter_len_conv))

    max_string_len = max([strlen(use_guess), strlen(initial_guess), strlen(iter_criterion), strlen(iter_crit_len), strlen(iter_nmax), strlen(iter_filter), strlen(iter_bwo), strlen(iter_len_conv), strlen(max_string_len)])


    while(strlen(use_guess) lt max_string_len) do use_guess = use_guess + pad_string
    while(strlen(initial_guess) lt max_string_len) do initial_guess = initial_guess + pad_string
    while(strlen(iter_criterion) lt max_string_len) do iter_criterion = iter_criterion + pad_string
    while(strlen(iter_crit_len) lt max_string_len) do iter_crit_len = iter_crit_len + pad_string
    while(strlen(iter_nmax) lt max_string_len) do iter_nmax = iter_nmax + pad_string
    while(strlen(iter_filter) lt max_string_len) do iter_filter = iter_filter + pad_string
    while(strlen(iter_bwo) lt max_string_len) do iter_bwo = iter_bwo + pad_string
    while(strlen(iter_len_conv) lt max_string_len) do iter_len_conv = iter_len_conv + pad_string



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
  printf, inp_lun, 'datadir         ', datadir                                                                , '# Full path of data directory'
  printf, inp_lun, 'galaxy          ', galaxy                                                             , '# Name of data set'
  printf, inp_lun, 'starfile        ', starfile                                                                , '# Name of primary stellar map'
  printf, inp_lun, 'starfile2       ', starfile2                                                               , '# Name of secondary stellar map (specifically for peak identification - only used if use_star2=1)'
  printf, inp_lun, 'gasfile         ', gasfile                                                                 , '# Name of primary gas map'
  printf, inp_lun, 'gasfile2        ', gasfile2                                                                       , '# Name of secondary gas map (specifically for peak identification - only used if use_gas2=1)'
  printf, inp_lun, 'starfile3       ', starfile3                                                                       , '# Name of additional stellar map (only used if use_star3=1)'

  printf, inp_lun, '# MASK FILE NAMES'
  printf, inp_lun, 'maskdir            ', maskdir,           '# Full path of mask directory'
  printf, inp_lun, 'star_ext_mask      ', star_ext_mask,     '# Region file for masking areas in starfile exterior to closed shapes (only used if mask_star=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask      ', star_int_mask,     '# Region file for masking areas in starfile interior to closed shapes (only used if mask_star=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_ext_mask       ', gas_ext_mask,      '# Region file for masking areas in gasfile exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_int_mask       ', gas_int_mask,      '# Region file for masking areas in gasfile interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_ext_mask2     ', star_ext_mask2,    '# Region file for masking areas in starfile2 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask2     ', star_int_mask2,    '# Region file for masking areas in starfile2 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_ext_mask2      ', gas_ext_mask2,     '# Region file for masking areas in gasfile2 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'gas_int_mask2      ', gas_int_mask2,     '# Region file for masking areas in gasfile2 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_ext_mask3     ', star_ext_mask3,    '# Region file for masking areas in starfile3 exterior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'
  printf, inp_lun, 'star_int_mask3     ', star_int_mask3,    '# Region file for masking areas in starfile3 interior to closed shapes (only used if mask_gas=1); should be a ds9 region file using the image coordinate system (other coordinate systems can be used if convert_masks=1, which calls ds9)'

  printf, inp_lun, '# FLAGS 1 (module switches)'
  printf, inp_lun, 'mask_images     ', mask_images,    '# Mask images (on/off)'
  printf, inp_lun, 'regrid          ', regrid,         '# Read and regrid original files and synchronise masks across all used images, i.e. pixels masked in one image will be masked in all images (on/off)'
  printf, inp_lun, 'smoothen        ', smoothen,       '# Create smoothened maps for each aperture size (on/off)'
  printf, inp_lun, 'sensitivity     ', sensitivity,    '# Fit Gaussians to pixel intensity PDFs to determine sensitivity limits (on/off)'
  printf, inp_lun, 'id_peaks        ', id_peaks,       '# Identify peaks and save their locations (on/off)'
  printf, inp_lun, 'calc_ap_flux    ', calc_ap_flux,   '# Calculate the enclosed flux for each peak and aperture size (on/off)'
  printf, inp_lun, 'generate_plot   ', generate_plot,  '# Generate output plots and save to PostScript files (on/off)'
  printf, inp_lun, 'get_distances   ', get_distances,  '# Calculate distances between all peak pairs (on/off)'
  printf, inp_lun, 'calc_obs        ', calc_obs,       '# Monte-Carlo sample peaks, get observed tuning fork (on/off)'
  printf, inp_lun, 'calc_fit        ', calc_fit,       '# Fit KL14 principle model (on/off)'
  printf, inp_lun, 'diffuse_frac    ', diffuse_frac,   '# Calculate diffuse fraction in images (on/off)'
  printf, inp_lun, 'derive_phys     ', derive_phys,    '# Calculate derived physical quantities (on/off)'
  printf, inp_lun, 'write_output    ', write_output,   '# Write results to output files (on/off)'
  printf, inp_lun, 'cleanup         ', cleanup,        '# [0] Keep auxiliary files to allow any of the above flags to be set to 0 for future runs; [1, DEFAULT] Delete auxiliary files after confirmation prompt; [2] WARNING: Automatically delete auxiliary files'
  printf, inp_lun, 'autoexit        ', autoexit,       '# Exit IDL automatically upon successful run completion (on/off)'

  printf, inp_lun, '# FLAGS 2 (use ancillary files)'
  printf, inp_lun, 'use_star2       ', use_star2,      '# Use different map for identifying stellar peaks and use the original for performing the flux calculation (on/off)'
  printf, inp_lun, 'use_gas2        ', use_gas2,       '# Use different map for identifying gas peaks and use the original for performing the flux calculation (on/off)'
  printf, inp_lun, 'use_star3       ', use_star3,      '# Use an additional stellar tracer map to mask pixels based on their flux ratios (on/off; e.g. Hα/FUV)'

  printf, inp_lun, '# FLAGS 3 (masking-related options)'
  printf, inp_lun, 'mstar_ext       ', mstar_ext,      '# Mask areas in starfile exterior to regions described in star_ext_mask (on/off)'
  printf, inp_lun, 'mstar_int       ', mstar_int,      '# Mask areas in starfile interior to regions described in star_int_mask (on/off)'
  printf, inp_lun, 'mgas_ext        ', mgas_ext,       '# Mask areas in gasfile exterior to regions described in gas_ext_mask (on/off)'
  printf, inp_lun, 'mgas_int        ', mgas_int,       '# Mask areas in gasfile interior to regions described in gas_int_mask (on/off)'
  printf, inp_lun, 'mstar_ext2      ', mstar_ext2,     '# Mask areas in starfile2 exterior to regions described in star_ext_mask2 (on/off)'
  printf, inp_lun, 'mstar_int2      ', mstar_int2,     '# Mask areas in starfile2 interior to regions described in star_int_mask2 (on/off)'
  printf, inp_lun, 'mgas_ext2       ', mgas_ext2,      '# Mask areas in gasfile2 exterior to regions described in gas_ext_mask2 (on/off)'
  printf, inp_lun, 'mgas_int2       ', mgas_int2,      '# Mask areas in gasfile2 interior to regions described in gas_int_mask2 (on/off)'
  printf, inp_lun, 'mstar_ext3      ', mstar_ext3,     '# Mask areas in starfile3 exterior to regions described in star_ext_mask3 (on/off)'
  printf, inp_lun, 'mstar_int3      ', mstar_int3,     '# Mask areas in starfile3 interior to regions described in star_int_mask3 (on/off)'
  printf, inp_lun, 'convert_masks   ', convert_masks,  '# Convert masks from other ds9 region coordinate systems to image coordinates (requires command line version of ds9) (on/off)'
  printf, inp_lun, 'cut_radius      ', cut_radius,     '# Mask maps outside a specified radial interval within the galaxy (on/off)'

  printf, inp_lun, '# FLAGS 4 (choose analysis options)'
  printf, inp_lun, 'set_centre      ', set_centre,     '# [0] Pixel coordinates of galaxy centre are set to the central pixel of starfile; [1, DEFAULT] Specify pixel coordinates of galaxy centre under INPUT PARAMETERS 1 (basic map analysis)'
  printf, inp_lun, 'tophat          ', tophat,         '# [0] Use Gaussian kernel to smoothen maps to larger aperture sizes; [1, DEFAULT] Use tophat kernel to smoothen maps to larger aperture sizes'
  printf, inp_lun, 'loglevels       ', loglevels,      '# [0] Linear equal-spacing of peak identification contour levels; [1, DEFAULT] Logarithmic equal-spacing of peak identification contour levels'
  printf, inp_lun, 'flux_weight     ', flux_weight,    '# [0, DEFAULT] Peak positions correspond to brightest pixel in peak area; [1] Peak positions correspond to flux-weighted mean position of peak area'
  printf, inp_lun, 'calc_ap_area    ', calc_ap_area,   '# Calculate the area of apertures used based on number of unmasked pixels within the target aperture defined under INPUT PARAMETERS 2 (apertures) (on/off)'
  printf, inp_lun, 'tstar_incl      ', tstar_incl,     '# [0, DEFAULT] Reference time-scale tstariso does not include the overlap phase (tstar=tstariso+tover); [1] tstariso includes the overlap phase (tstar=tstariso)'
  printf, inp_lun, 'peak_prof       ', peak_prof,      '# [0] Model independent regions as points; [1] Model independent regions as constant-surface density discs; [2, DEFAULT] Model independent regions as two-dimensional Gaussians'
  printf, inp_lun, 'map_units       ', map_units,      '# [0] Unknown, do not calculate derived quantities; [1, DEFAULT] star=SFR and gas=gas; [2] star=gas1 and gas=gas2; [3] star=SFR1 and gas=SFR2'
  printf, inp_lun, 'use_X11         ', use_X11,        '# [0] Not allowed to create X11 windows; [1] Allowed to create X11 windows'







  printf, inp_lun, '# INPUT PARAMETERS 1 (basic map analysis)'
  printf, inp_lun, 'distance        ', distance,       '# Distance to galaxy in pc'
  printf, inp_lun, 'inclination     ', inclination,    '# Inclination angle in degrees'
  printf, inp_lun, 'posangle        ', posangle,       '# Position angle in degrees'
  printf, inp_lun, 'centrex         ', centrex,        '# Index of pixel x-axis coordinate of galaxy centre in starfile (important: start counting at zero on the left-hand side of the image)'
  printf, inp_lun, 'centrey         ', centrey,        '# Index of pixel y-axis coordinate of galaxy centre in starfile (important: start counting at zero on the bottom side of the image)'
  printf, inp_lun, 'minradius       ', minradius,      '# Minimum radius for analysis in pc'
  printf, inp_lun, 'maxradius       ', maxradius,      '# Maximum radius for analysis in pc'
  printf, inp_lun, 'Fs1_Fs2_min     ', Fs1_Fs2_min,    '# Minimum primary-to-secondary star formation tracer flux ratio for including pixels (only if use_star3=1)'
  printf, inp_lun, 'max_sample      ', max_sample,     '# Maximum number of pixels per map resolution FWHM'
  printf, inp_lun, 'astr_tolerance  ', astr_tolerance, '# allowable tolerance in decimal degrees between astrometric position of image pixels of different images (e.g. 0.9" = 2.5e-4 decimal degrees)'
  printf, inp_lun, 'nbins           ', nbins,          '# Number of bins used during sensitivity limit calculation to sample the intensity PDFs and fit Gaussians'

  printf, inp_lun, '# INPUT PARAMETERS 2 (apertures)'
  printf, inp_lun, 'lapmin          ', lapmin,         '# Minimum aperture size (i.e. diameter) in pc to create smoothened maps for'
  printf, inp_lun, 'lapmax          ', lapmax,         '# Maximum aperture size (i.e. diameter) in pc to create smoothened maps for'
  printf, inp_lun, 'naperture       ', naperture,      '# Number of aperture sizes'
  printf, inp_lun, 'peak_res        ', peak_res,       '# Minimum aperture size used in fitting the KL14 principle model, also index of aperture size at which peaks are identified (start counting at 0) - is set to map resolution if originally chosen to be smaller'
  printf, inp_lun, 'max_res         ', max_res,        '# Maximum aperture size used in fitting the KL14 principle model, also index of aperture size at which fluxes best reflect galactic averages (start counting at 0)'


  printf, inp_lun, '# INPUT PARAMETERS 3 (peak identification)'
  printf, inp_lun, 'npixmin         ', npixmin,        '# Minimum number of pixels for a valid peak (use npixmin=1 to allow points)'
  printf, inp_lun, 'nsigma          ', nsigma,         '# Multiple of the sensitivity limit needed for a valid peak'
  printf, inp_lun, 'logrange_s      ', logrange_s,     '# Logarithmic range covered by flux contour levels for stellar peak identification, counting down from stellar flux maximum'
  printf, inp_lun, 'logspacing_s    ', logspacing_s,   '# Logarithmic interval between flux contour levels for stellar peak identification'
  printf, inp_lun, 'logrange_g      ', logrange_g,     '# Logarithmic range covered by flux contour levels for gas peak identification, counting down from gas flux maximum'
  printf, inp_lun, 'logspacing_g    ', logspacing_g,   '# Logarithmic interval between flux contour levels for gas peak identification'
  printf, inp_lun, 'nlinlevel_s     ', nlinlevel_s,    '# Number of flux contour levels for stellar peak identification between zero and stellar flux maximum (only if loglevels=0)'
  printf, inp_lun, 'nlinlevel_g     ', nlinlevel_g,    '# Number of flux contour levels for gas peak identification between zero and gas flux maximum (only if loglevels=0)'



  printf, inp_lun, '# INPUT PARAMETERS 4 (timeline)'
  printf, inp_lun, 'tstariso        ', tstariso,            '# Reference time-scale spanned by star formation tracer in Myr'
  printf, inp_lun, 'tstariso_errmin ', tstariso_errmin,     '# Downward standard error of the reference time-scale in Myr'
  printf, inp_lun, 'tstariso_errmax ', tstariso_errmax,     '# Upward standard error of the reference time-scale in Myr'
  printf, inp_lun, 'tgasmini        ', tgasmini,            '# Minimum value of tgas considered during fitting process in Myr'
  printf, inp_lun, 'tgasmaxi        ', tgasmaxi,            '# Maximum value of tgas considered during fitting process in Myr'
  printf, inp_lun, 'tovermini       ', tovermini,           '# Minimum value of tover considered during fitting process in Myr'



  printf, inp_lun, '# INPUT PARAMETERS 5 (fitting)'
  printf, inp_lun, 'nmc             ', nmc,            '# Number of Monte-Carlo peak drawing experiments to be executed and averaged over during fitting process'
  printf, inp_lun, 'ndepth          ', ndepth,         '# Maximum number of free parameter array refinement loops for obtaining best-fitting value'
  printf, inp_lun, 'ntry            ', ntry,           '# Size of each free parameter array to obtain the best-fitting value'
  printf, inp_lun, 'nphysmc         ', nphysmc,        '# Number of Monte-Carlo experiments for error propagation of derived physical quantities'



  printf, inp_lun, '# INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)''
  printf, inp_lun, 'diffuse_quant   ', diffuse_quant,  '# Calculate the diffuse fraction of [0] flux, [1] power'
  printf, inp_lun, 'f_filter_type   ', f_filter_type,  '# Filter type choice for Fourier filtering: [0] butterworth filter [1] gaussian filter [2] ideal filter (Only relevant if diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
  printf, inp_lun, 'bw_order        ', bw_order,       '# The order of the butterworth filter (Must be an integer. Only relevant if f_filter_type = 0 [butterworth] and diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
  printf, inp_lun, 'filter_len_conv ', filter_len_conv,'# conversion factor for Fourier filtering cut_length (cut_length = lambda * filter_len_conv)'



  printf, inp_lun, '# INPUT PARAMETERS 7 (conversions and constants to calculate derived quantities)'
  printf, inp_lun, 'convstar        ', convstar,       '# Log of multiplication factor to convert the pixel value in starfile to a SFR in Msun/yr (if maptypes=1 or maptypes=3) or to a gas mass in Msun (if maptypes=2)'
  printf, inp_lun, 'convstar_rerr   ', convstar_rerr,  '# Relative error (sigma_x/x) of convstar'
  printf, inp_lun, 'convgas         ', convgas,        '# Log of multiplication factor to convert the pixel value in gasfile to a gas mass in Msun (if maptypes=1 or maptypes=2) or to a SFR in Msun/yr (if maptypes=3)'
  printf, inp_lun, 'convgas_rerr    ', convgas_rerr,   '# Relative error (sigma_x/x) of convgas'
  printf, inp_lun, 'convstar3       ', convstar3,      '# Log of multiplication factor to convert the pixel value in starfile3 to a SFR in Msun/yr'
  printf, inp_lun, 'convstar3_rerr  ', convstar3_rerr, '# Relative error (sigma_x/x) of convstar3'
  printf, inp_lun, 'lighttomass     ', lighttomass,    '# Light-to-mass ratio of a desired feedback mechanism in m^2 s^-3 (default 0.002 is SNe)'
  printf, inp_lun, 'momratetomass   ', momratetomass,  '# Momentum output rate per unit mass of a desired feedback mechanism in m s^-2 (default 5.e-10 is SNe+winds)'


  printf, inp_lun, '# INPUT PARAMETERS 8 (sensitivity)'
  printf, inp_lun, 'use_stds        ', use_stds,       '# [0] calculate standard deviations of images [1] Use supplied standard deviation measurements'
  printf, inp_lun, 'std_star        ', std_star,       '# The measured standard deviation of starfile'
  printf, inp_lun, 'std_star3       ', std_star3,      '# The measured standard deviation of gasfile'
  printf, inp_lun, 'std_gas         ', std_gas,        '# The measured standard deviation of starfile3 (only used of use_star3 = 1)'

  if iter_input_switch eq 1 then begin ; create input file for iteration
    printf, inp_lun, '# INPUT PARAMETERS 9 (Fourier diffuse removal iteration)''
    printf, inp_lun, 'use_guess       ', use_guess,      '#  [0] Do not Filter images before first KL14 run and ignore parameter: "initial_guess" [1] Filter images with initial guess of lambda before first KL14 run'
    printf, inp_lun, 'initial_guess   ', initial_guess,  '#  Length in pc of the initial estimate of lambda with which to filter the images'
    printf, inp_lun, 'iter_criterion  ', iter_criterion, '#  Fractional difference between lambda and previous lambda value(s) that triggers an end to the iteration process.'
    printf, inp_lun, 'iter_crit_len   ', iter_crit_len,  '#  Number of previous runs of KL14 over which iter_criterion must be true, i.e. if iter_crit_len = 1 then the current value of lambda is compared to the previous value of lambda, if iter_crit_len = 2 then then the current value of lambda is compared to the previous value of lambda and the value of lambda prior to the previous value'
    printf, inp_lun, 'iter_nmax       ', iter_nmax,      '#  Maximum number of KL14 runs allowed in the iteration process. If iter_criterion is not reached before this number of runs the programme will exit anyway.'
    printf, inp_lun, 'iter_filter     ', iter_filter,    '#  Filter type choice for iterative Fourier filtering: [0] butterworth filter [1] gaussian filter [2] ideal filter (Only relevant if diffuse_quant = 1 [flux], but a dummy value should be supplied anyway)'
    printf, inp_lun, 'iter_bwo        ', iter_bwo,       '#  The order of the butterworth filter for iterative Fourier filtering (Must be an integer. Only relevant if iter_bwo = 0 [butterworth], but a dummy value should be supplied anyway)'
    printf, inp_lun, 'iter_len_conv   ', iter_len_conv,  '#  conversion factor for iterative Fourier filtering cut_length (cut_length = lambda * filter_len_conv)'
  endif

  printf, inp_lun, '################################################'
  printf, inp_lun, '#                                              #'
  printf, inp_lun, '#              END PARAMETER FILE              #'
  printf, inp_lun, '#                                              #'
  printf, inp_lun, '################################################'



  free_lun, inp_lun


end
