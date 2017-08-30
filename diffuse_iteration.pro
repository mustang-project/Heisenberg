pro diffuse_iteration, master_inputfile
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range

  master_inputfile = master_inputfile[0]

  input_split = strsplit(master_inputfile, path_sep(), /extract)
  ifile_shortname = input_split[-1]

  ; #############################
  ; COPY from tuningfork.pro -> needs modularisation
  ; #############################


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;READ INPUT FILE AND VERIFY;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  openr,lun,master_inputfile,/get_lun ;open the input file
  line='' ;initialise line
  print,' ==> reading '+master_inputfile+' with the following variables:'
  while ~eof(lun) do begin ;read lines until end string is reached
      readf,lun,line ;read line containing variable
      linearr=strsplit(line,' ',/extract) ;split line into space-separated parts
      if n_elements(linearr) gt 1 && strmid(linearr[0],0,1) ne '#' then begin ;if line is not empty or comment
          value=linearr[1]
          is_num=strnumber(value,value)
          (scope_varfetch(linearr[0],/enter))=value ;assign variable name and value
          spacelen=20-strlen(linearr[0])
          space=''
          for i=0,spacelen-1 do space+=' '
          print,'     '+linearr[0]+space+linearr[1]
      endif
  endwhile
  close,lun ;close the logical unit
  free_lun,lun ;make the logical unit available again

  expected_flags1=['mask_images','regrid','smoothen','sensitivity','id_peaks','calc_ap_flux','generate_plot','get_distances','calc_obs','calc_fit','diffuse_frac','derive_phys','write_output','cleanup','autoexit'] ;variable names of expected flags (1)
  expected_flags2=['use_star2','use_gas2','use_star3'] ;variable names of expected flags (2)
  expected_flags3=['mstar_ext','mstar_int','mgas_ext','mgas_int','mstar_ext2','mstar_int2','mgas_ext2','mgas_int2','mstar_ext3','mstar_int3','convert_masks','cut_radius'] ;variable names of expected flags (3)
  expected_flags4=['set_centre','tophat','loglevels','flux_weight','calc_ap_area','tstar_incl','peak_prof','map_units','use_X11'] ;variable names of expected flags (4)
  expected_flags=[expected_flags1,expected_flags2,expected_flags3,expected_flags4] ;variable names of expected flags (all)
  expected_filenames=['datadir','galaxy','starfile','starfile2','gasfile','gasfile2','starfile3'] ;variable names of expected filenames
  expected_masknames=['maskdir','star_ext_mask','star_int_mask','gas_ext_mask','gas_int_mask','star_ext_mask2','star_int_mask2','gas_ext_mask2','gas_int_mask2','star_ext_mask3','star_int_mask3'] ;variable names of expected mask filenames
  expected_params1=['distance','inclination','posangle','centrex','centrey','minradius','maxradius','Fs1_Fs2_min','max_sample','astr_tolerance','nbins'] ;variable names of expected input parameters (1)
  expected_params2=['lapmin','lapmax','naperture','peak_res','max_res'] ;variable names of expected input parameters (2)
  expected_params3=['npixmin','nsigma','logrange_s','logspacing_s','logrange_g','logspacing_g'] ;variable names of expected input parameters (3)
  expected_params4=['tstariso','tstariso_errmin','tstariso_errmax','tgasmini','tgasmaxi','tovermini'] ;variable names of expected input parameters (4)
  expected_params5=['nmc','ndepth','ntry','nphysmc'] ;variable names of expected input parameters (5)
  expected_params6=['diffuse_quant','filter_len_conv','f_filter_type','bw_order'] ;variable names of expected input parameters (6)
  expected_params7=['convstar','convstar_rerr','convgas','convgas_rerr','convstar3','convstar3_rerr','lighttomass','momratetomass'] ;variable names of expected input parameters (7)
  expected_params8=['use_stds','std_star','std_star3','std_gas'] ;variable names of expected input parameters (8)
  expected_params9=['use_guess','initial_guess','iter_criterion','iter_crit_len','iter_nmax','iter_filter','iter_bwo','iter_len_conv'] ;variable names of expected input parameters (8)
  expected_params=[expected_params1,expected_params2,expected_params3,expected_params4,expected_params5,expected_params6,expected_params7,expected_params8,expected_params9] ;variable names of expected input parameters (all)
  expected_vars=[expected_flags,expected_filenames,expected_masknames,expected_params] ;names of all expected variables
  nvars=n_elements(expected_vars) ;number of expected variables
  for i=0,nvars-1 do if n_elements(scope_varfetch(expected_vars[i])) eq 0 then f_error,'variable '+expected_vars[i]+' not present in input file' ;verify input reading
  if peak_prof eq 0 and npixmin gt 1 then f_error,['peak_prof=0 (points), but minimum number of pixels npixmin>1','Please set npixmin=1 or peak_prof>0']


  ; #############################
  ; COPY from tuningfork.pro -> needs modularisation
  ; #############################


  ; ******************************************************
  ; create master directory for iterations
  ; ******************************************************

  master_inputdir=file_dirname(master_inputfile)+path_sep() ;input directory
  master_runname=file_basename(master_inputfile)+'_iterrun' ;main run directory name
  master_rundir=strcompress(master_inputdir+master_runname+path_sep(),/remove_all) ;main run directory full path

  dummy=file_search(master_rundir,count=ct) ;check if rundir exists
  if ct eq 0 then spawn,'mkdir '+master_rundir ;if not, create it
  if ct ne 0 then print,' WARNING: iteration run directory already exists, some or all files may be overwritten'



  ; ******************************************************
  ; get key information about original images
  ; ******************************************************

  master_datadir = strcompress(datadir,/remove_all)
  if strmid(master_datadir, 0, /reverse) ne path_sep() then master_datadir = master_datadir + path_sep() ; directory check

  datadir = strcompress(master_rundir + 'Fourier_datadir' + path_sep(), /remove_all) ; new datadir for Fourier_filtered data
  dummy=file_search(datadir,count=ct) ;check if rundir exists
  if ct eq 0 then spawn,'mkdir '+datadir ;if not, create it
  if ct ne 0 then print,' WARNING: iteration data directory already exists, some or all files may be overwritten'


  master_starfiletot = master_datadir + starfile
  master_gasfiletot = master_datadir + gasfile

  master_starhdr = headfits(master_starfiletot)
  master_gashdr = headfits(master_gasfiletot)
  ; ********************************
  ; CONVOLVE TO COMMON BEAM!!!!
  ; ********************************


  getrot, master_starhdr, rotation, cdelt_var ; get cdelt value ; only single precision
  if n_elements(cdelt_var) eq 2 && abs(abs(cdelt_var[0]) - abs(cdelt_var[1])) le astr_tolerance then begin
    master_star_cdelt = mean(abs(cdelt_var))
  endif else begin
    print, "star image not square"
    stop
  endelse

  master_star_pix_to_pc = distance*tan(!dtor*master_star_cdelt)/sqrt(cos(inclination))



  getrot, master_gashdr, rotation, cdelt_var ; get cdelt value ; only single precision
  if n_elements(cdelt_var) eq 2 && abs(abs(cdelt_var[0]) - abs(cdelt_var[1])) le astr_tolerance then begin
    master_gas_cdelt = mean(abs(cdelt_var))
  endif else begin
    print, "gas image not square"
    stop
  endelse

  master_gas_pix_to_pc = distance*tan(!dtor*master_gas_cdelt)/sqrt(cos(inclination))




  ; ******************************************************
  ; apply initial guess if requested
  ; ******************************************************

  if use_guess eq 1 then begin ; initiial fourier cut
    print, "functionality not implemented"
    stop
  endif else begin ; copy original files

    starfile = "starfile_iter0.fits"
    gasfile = "gasfile_iter0.fits"

    base_starfile = datadir + starfile ; the base file for fourier iterations
    base_gasfile = datadir + gasfile ; the base file for fourier iterations


    file_copy, master_starfiletot, base_starfile, /overwrite
    file_copy, master_gasfiletot, base_gasfile, /overwrite

  endelse


  ; ************************************************
  ; *** Cut runs
  ; ************************************************
  iter_tgas_vec = dblarr(iter_nmax) ; initialise vectors for results             [] ; initialise blank vecotr
  iter_tgas_errmin_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_tgas_errmax_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_tover_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_tover_errmin_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_tover_errmax_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_lambda_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_lambda_errmin_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr
  iter_lambda_errmax_vec = dblarr(iter_nmax) ;                                    [] ; initialise blank vecotr


  ; ******************************************************
  ; iteration loop
  ; ******************************************************

  iter_report_file = master_rundir + ifile_shortname + '_iteration_report.dat'

  openw, iter_lun, iter_report_file, width = 450, /get_lun
  printf, iter_lun, '  tgas [Myr]       ',  '  tgas_errmin [Myr]       ',  '  tgas_errmax [Myr]       ',  '  tover [Myr]       ',  '  tover_errmin [Myr]       ',  '  tover_errmax [Myr]       ',  '  lambda [pc]        ',  '  lambda_errmin [pc]        ',  '  lambda_errmax [pc]        '

  iter_break = 0
  iter_num = 0
  while(iter_break eq 0 && iter_num lt iter_nmax) do begin



    input_file_name = strcompress(master_inputfile + ifile_shortname + '_iter' + string(iter_num), /remove_all)
    input_file_filepath = strcompress(master_rundir + input_file_name, /remove_all)
    ; **********************
    ; * generate input file
    ; **********************
    make_input_file, input_file_filepath   $ ; general variables
     , datadir, galaxy, starfile, gasfile   $  ; File names compulsory variables
     , distance, inclination, posangle, centrex, centrey, minradius, maxradius, Fs1_Fs2_min, max_sample, astr_tolerance, nbins $ ; # INPUT PARAMETERS 1 (basic map analysis)
     , lapmin, lapmax, naperture, peak_res, max_res $  ; # INPUT PARAMETERS 2 (apertures)
     , starfile2 = starfile2, gasfile2 = gasfile2, starfile3 = starfile3  $  ; File names keywords
     , maskdir = maskdir, star_ext_mask1 = star_ext_mask, star_int_mask1 = star_int_mask, gas_ext_mask1 = gas_ext_mask, gas_int_mask1 = gas_int_mask, star_ext_mask2 = star_ext_mask2, star_int_mask2 = star_int_mask2, gas_ext_mask2 = gas_ext_mask2, gas_int_mask2 = gas_int_mask2, star_ext_mask3 = star_ext_mask3, star_int_mask3 = star_int_mask3  $  ; Mask file names keywords
     , mask_images = mask_images, regrid = regrid, smoothen = smoothen, sensitivity = sensitivity, id_peaks = id_peaks, calc_ap_flux = calc_ap_flux, generate_plot = generate_plot, get_distances = get_distances, calc_obs = calc_obs, calc_fit = calc_fit, diffuse_frac = diffuse_frac, derive_phys = derive_phys, write_output = write_output, cleanup = cleanup, autoexit = autoexit $ ; FLAGS 1 (keywords)'
     , use_star2 = use_star2, use_gas2 = use_gas2, use_star3 = use_star3 $ ;  FLAGS 2 keywords
     , mstar_ext1 = mstar_ext, mstar_int1 = mstar_int, mgas_ext1 = mgas_ext, mgas_int1 = mgas_int, mstar_ext2 = mstar_ext2, mstar_int2 = mstar_int2, mgas_ext2 = mgas_ext2, mgas_int2 = mgas_int2, mstar_ext3 = mstar_ext3, mstar_int3 = mstar_int3, convert_masks = convert_masks, cut_radius = cut_radius $ ; # FLAGS 3 (masking-related options) keywords ; note mstar_ext1 = mstar_ext etc. prevents ambigious keyword error
     , set_centre = set_centre, tophat = tophat, loglevels = loglevels, flux_weight = flux_weight, calc_ap_area = calc_ap_area, tstar_incl = tstar_incl, peak_prof = peak_prof, map_units = map_units, use_X11 = use_X11 $ ; # FLAGS 4 (choose analysis options)
     , npixmin = npixmin, nsigma = nsigma, logrange_s = logrange_s, logspacing_s = logspacing_s, logrange_g = logrange_g, logspacing_g = logspacing_g, nlinlevel_s = nlinlevel_s, nlinlevel_g = nlinlevel_g $ ; # INPUT PARAMETERS 3 (peak identification)
     , tstariso_val = tstariso, tstariso_errmin = tstariso_errmin, tstariso_errmax = tstariso_errmax, tgasmini = tgasmini, tgasmaxi = tgasmaxi, tovermini = tovermini $ ; # INPUT PARAMETERS 4 (timeline) ; note tstariso_val = tstariso prevents ambigious keyword error
     , nmc = nmc, ndepth = ndepth, ntry = ntry, nphysmc = nphysmc $ ; # INPUT PARAMETERS 5 (fitting)
     , diffuse_quant, f_filter_type, bw_order, filter_len_conv $ ; # INPUT PARAMETERS 6 (Fourier filtering for diffuse gas calculation)
     , convstar_val = convstar, convstar_rerr = convstar_rerr, convgas_val = convgas, convgas_rerr = convgas_rerr, convstar3_val = convstar3, convstar3_rerr = convstar3_rerr, lighttomass = lighttomass, momratetomass = momratetomass $ ; # INPUT PARAMETERS 6 (conversions and constants to calculate derived quantities) ; note convgas_val = convgas avoids % Ambiguous keyword abbreviation: CONVGAS.
     , use_stds = use_stds, std_star1 = std_star, std_star3 = std_star3, std_gas = std_gas ; # INPUT PARAMETERS 8 (sensitivity)

     ; **********************
     ; * call KL14
     ; **********************
     spawn, 'idl kl14.pro -arg ' + input_file_filepath
     ; **********************
     ; * get variables
     ; **********************
     output_filename = strcompress(input_file_filepath+ '_run/output/' + galaxy + '_output.dat', /remove_all)
     read_kl14_tablerow, output_filename, name_vec, value_vec, /de_log, /compress_names ; get values, convert from stored log value to the real value and remove whitespace from names

     var_loc = where(name_vec eq 'tgas', var_count)
     if var_count eq 1 then iter_tgas_vec[iter_num] = value_vec[var_loc] else begin
     ;if var_loc ne !null then iter_tgas_vec = [iter_tgas_vec, value_vec[var_loc]] else begin
       print, 'variable tgas not found'
       stop
     endelse

    var_loc = where(name_vec eq 'tgas_errmin', var_count)
    if var_count eq 1 then iter_tgas_errmin_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_tgas_errmin_vec = [iter_tgas_errmin_vec, value_vec[var_loc]] else begin
     print, 'variable tgas_errmin not found'
     stop
    endelse

    var_loc = where(name_vec eq 'tgas_errmax', var_count)
    if var_count eq 1 then iter_tgas_errmax_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_tgas_errmax_vec = [iter_tgas_errmax_vec, value_vec[var_loc]] else begin
     print, 'variable tgas_errmax not found'
     stop
    endelse


    var_loc = where(name_vec eq 'tover', var_count)
    if var_count eq 1 then iter_tover_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_tover_vec = [iter_tover_vec, value_vec[var_loc]] else begin
     print, 'variable tover not found'
     stop
    endelse

    var_loc = where(name_vec eq 'tover_errmin', var_count)
    if var_count eq 1 then iter_tover_errmin_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_tover_errmin_vec = [iter_tover_errmin_vec, value_vec[var_loc]] else begin
     print, 'variable tover_errmin not found'
     stop
    endelse

    var_loc = where(name_vec eq 'tover_errmax', var_count)
    if var_count eq 1 then iter_tover_errmax_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_tover_errmax_vec = [iter_tover_errmax_vec, value_vec[var_loc]] else begin
     print, 'variable tover_errmax not found'
     stop
    endelse


    var_loc = where(name_vec eq 'lambda', var_count)
    if var_count eq 1 then iter_lambda_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_lambda_vec = [iter_lambda_vec, value_vec[var_loc]] else begin
     print, 'variable lambda not found'
     stop
    endelse

    var_loc = where(name_vec eq 'lambda_errmin', var_count)
    if var_count eq 1 then iter_lambda_errmin_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_lambda_errmin_vec = [iter_lambda_errmin_vec, value_vec[var_loc]] else begin
     print, 'variable lambda_errmin not found'
     stop
    endelse

    var_loc = where(name_vec eq 'lambda_errmax', var_count)
    if var_count eq 1 then iter_lambda_errmax_vec[iter_num] = value_vec[var_loc] else begin
    ;if var_loc ne !null then iter_lambda_errmax_vec = [iter_lambda_errmax_vec, value_vec[var_loc]] else begin
     print, 'variable lambda_errmax not found'
     stop
    endelse


    ; ; ************************************************
    ; ; *** output
    ; ; ************************************************
    printf, iter_lun,  iter_tgas_vec[iter_num], iter_tgas_errmin_vec[iter_num], iter_tgas_errmax_vec[iter_num], iter_tover_vec[iter_num], iter_tover_errmin_vec[iter_num], iter_tover_errmax_vec[iter_num], iter_lambda_vec[iter_num], iter_lambda_errmin_vec[iter_num], iter_lambda_errmax_vec[iter_num]


    plot_filename = strcompress(master_rundir + ifile_shortname + '_tgas_iteration.ps', /remove_all)
    iteration_plot, plot_filename, iter_tgas_vec[0:iter_num], iter_tgas_errmin_vec[0:iter_num], iter_tgas_errmax_vec[0:iter_num], 'Tgas [Myr]'

    plot_filename = strcompress(master_rundir + ifile_shortname + '_tover_iteration.ps', /remove_all)
    iteration_plot, plot_filename, iter_tover_vec[0:iter_num], iter_tover_errmin_vec[0:iter_num], iter_tover_errmax_vec[0:iter_num], 'Tover [Myr]'

    plot_filename = strcompress(master_rundir + ifile_shortname + '_lambda_iteration.ps', /remove_all)
    iteration_plot, plot_filename, iter_lambda_vec[0:iter_num], iter_lambda_errmin_vec[0:iter_num], iter_lambda_errmax_vec[0:iter_num], 'Lambda [pc]'



    lambda_val = iter_lambda_vec[iter_num] ; lambda value for this iteration
    if iter_num ge iter_crit_len  then begin ; check iteration condition if nim number of iterations reached

      iter_break = 1
      for cc = 1,iter_crit_len, 1 do begin
        if ~(abs((lambda_val - iter_lambda_vec[iter_num-cc]) / lambda_val) lt iter_criterion) then begin
          iter_break = 0
        endif
      endfor

      if iter_break eq 1 then begin
        print, "Iteration condition reahced. Ending iteration"
        break ; end Fourier iteration
      endif
    endif


    ; ********************************************
    ; * Fourier filter master images images
    ; ********************************************
    iter_num ++ ; increase iteration counter


    if iter_filter eq 0 then filter_choice = 'butterworth' else $
        if iter_filter eq 1 then filter_choice = 'gaussian' else $
        if iter_filter eq 2 then filter_choice = 'ideal'

    star_cut_length = (lambda_val/master_star_pix_to_pc) * iter_len_conv ; in pixels
    starfile = strcompress("starfile_iter" + string(iter_num) + ".fits",/remove_all) ; change starfile for input_file
    filtered_starfile = strcompress(datadir + starfile,/remove_all)

    fourier_filter_tool, filter_choice, iter_bwo, 'high', star_cut_length, $ ; description of the filter
      base_starfile, $ ; input image by filename
      filtered_image_path = filtered_starfile ; output file for the filtered image

    star_arr = readfits(filtered_starfile, star_hdr)
    szero_list = where(star_arr lt 0.0e0, szero_count)
    if szero_count ne 0 then star_arr[szero_list] = 0.0
    writefits, filtered_starfile, star_arr, star_hdr



    gas_cut_length = (lambda_val/master_gas_pix_to_pc) * iter_len_conv ; in pixels
    gasfile = strcompress("gasfile_iter" + string(iter_num) + ".fits",/remove_all) ; change gasfile for input_file
    filtered_gasfile = strcompress(datadir + gasfile,/remove_all)

    fourier_filter_tool, filter_choice, iter_bwo, 'high', star_cut_length, $ ; description of the filter
      base_gasfile, $ ; input image by filename
      filtered_image_path = filtered_gasfile ; output file for the filtered image

    gas_arr = readfits(filtered_gasfile, gas_hdr)
    gzero_list = where(gas_arr lt 0.0e0, gzero_count)
    if gzero_count ne 0 then gas_arr[gzero_list] = 0.0
    writefits, filtered_gasfile, gas_arr, gas_hdr



  endwhile

  free_lun, iter_lun ; free output_file



end
