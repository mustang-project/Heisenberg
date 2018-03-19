pro read_kl14_input_file, input_file, $ ; input_file_filepath
  mask_images, regrid, smoothen, sensitivity,id_peaks, calc_ap_flux ,generate_plot, get_distances, calc_obs, calc_fit, diffuse_frac, derive_phys, write_output, cleanup , autoexit, $ ;variable names of expected flags (1)
  use_star2 ,use_gas2 ,use_star3, $  ;variable names of expected flags (2)
  mstar_ext, mstar_int, mgas_ext, mgas_int,mstar_ext2, mstar_int2 ,mgas_ext2, mgas_int2, mstar_ext3, mstar_int3, convert_masks, cut_radius, $ ;variable names of expected flags (3)
  set_centre, tophat, loglevels, peak_find_tui,flux_weight, calc_ap_area ,tstar_incl, peak_prof, map_units, use_X11, log10_output, $;variable names of expected flags (4)
  datadir, galaxy, starfile, starfile2,gasfile, gasfile2 ,starfile3, $ ;variable names of expected filenames
  maskdir, star_ext_mask, star_int_mask, gas_ext_mask,gas_int_mask, star_ext_mask2 ,star_int_mask2, gas_ext_mask2, gas_int_mask2, star_ext_mask3, star_int_mask3, $ ;variable names of expected mask filenames
  unfiltdir, star_unfilt_file, gas_unfilt_file, $
  ; parameters:
  distance, inclination, posangle, centrex,centrey, minradius ,maxradius, Fs1_Fs2_min, max_sample, astr_tolerance, nbins, $ ;variable names of expected input parameters (1)
  lapmin, lapmax, naperture, peak_res,max_res, $  ;variable names of expected input parameters (2)
  npixmin, nsigma, logrange_s, logspacing_s,logrange_g, logspacing_g, nlinlevel_s, nlinlevel_g, $ ;variable names of expected input parameters (3)
  tstariso, tstariso_errmin, tstariso_errmax, tgasmini,tgasmaxi, tovermini, $ ;variable names of expected input parameters (4)
  nmc, ndepth, ntry, nphysmc, $ ;variable names of expected input parameters (5)
  use_unfilt_ims, diffuse_quant, filter_len_conv, emfrac_cor_mode, f_filter_type, bw_order, $ ;variable names of expected input parameters (6)
  convstar, convstar_rerr, convgas, convgas_rerr,convstar3, convstar3_rerr ,lighttomass, momratetomass, $ ;variable names of expected input parameters (7)
  use_stds, std_star, std_star3, std_gas, $ ;variable names of expected input parameters (8)
  use_guess, initial_guess, iter_criterion, iter_crit_len,iter_nmax, iter_filter ,iter_bwo, iter_len_conv, iter_autoexit, use_nice, nice_value ;variable names of expected input parameters (9)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;READ INPUT FILE AND VERIFY;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  openr,lun,input_file,/get_lun ;open the input file
  line='' ;initialise line
  print,' ==> reading '+input_file+' with the following variables:'
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
  expected_flags4=['set_centre','tophat','loglevels','peak_find_tui','flux_weight','calc_ap_area','tstar_incl','peak_prof','map_units','use_X11','log10_output'] ;variable names of expected flags (4)
  expected_flags=[expected_flags1,expected_flags2,expected_flags3,expected_flags4] ;variable names of expected flags (all)

  expected_filenames=['datadir','galaxy','starfile','starfile2','gasfile','gasfile2','starfile3'] ;variable names of expected filenames
  expected_masknames=['maskdir','star_ext_mask','star_int_mask','gas_ext_mask','gas_int_mask','star_ext_mask2','star_int_mask2','gas_ext_mask2','gas_int_mask2','star_ext_mask3','star_int_mask3'] ;variable names of expected mask filenames
  expected_unfiltered=['unfiltdir','star_unfilt_file','gas_unfilt_file'] ;variable names of unfiltered filenames

  expected_params1=['distance','inclination','posangle','centrex','centrey','minradius','maxradius','Fs1_Fs2_min','max_sample','astr_tolerance','nbins'] ;variable names of expected input parameters (1)
  expected_params2=['lapmin','lapmax','naperture','peak_res','max_res'] ;variable names of expected input parameters (2)
  expected_params3=['npixmin','nsigma','logrange_s','logspacing_s','logrange_g','logspacing_g','nlinlevel_s','nlinlevel_g'] ;variable names of expected input parameters (3)
  expected_params4=['tstariso','tstariso_errmin','tstariso_errmax','tgasmini','tgasmaxi','tovermini'] ;variable names of expected input parameters (4)
  expected_params5=['nmc','ndepth','ntry','nphysmc'] ;variable names of expected input parameters (5)
  expected_params6=['use_unfilt_ims','diffuse_quant','filter_len_conv','emfrac_cor_mode','f_filter_type','bw_order'] ;variable names of expected input parameters (6)
  expected_params7=['convstar','convstar_rerr','convgas','convgas_rerr','convstar3','convstar3_rerr','lighttomass','momratetomass'] ;variable names of expected input parameters (7)
  expected_params8=['use_stds','std_star','std_star3','std_gas'] ;variable names of expected input parameters (8)
  expected_params=[expected_params1,expected_params2,expected_params3,expected_params4,expected_params5,expected_params6,expected_params7,expected_params8] ;variable names of expected input parameters (all)

  n_min_vars = n_elements([expected_flags,expected_filenames,expected_masknames,expected_unfiltered,expected_params]) + 1 ; plus 1 for input file
  ; check for optional input parameters:
  expected_params9=['use_guess','initial_guess','iter_criterion','iter_crit_len','iter_nmax','iter_filter','iter_bwo','iter_len_conv','iter_autoexit', 'use_nice', 'nice_value'] ;variable names of expected input parameters (9)
  if n_params() eq (n_min_vars + n_elements(expected_params9)) then begin
    expected_params = [expected_params,expected_params9]
  endif else if n_params() gt n_min_vars then begin
    f_error,['Incorrect number of parameters supplied to read_kl14_input_file']
  endif




  ;
  ; if n_params eq 117 then begin; check if user has requested iteration parameters
  ;   expected_params9=['use_guess','initial_guess','iter_criterion','iter_crit_len','iter_nmax','iter_filter','iter_bwo','iter_len_conv','iter_autoexit'] ;variable names of expected input parameters (9)
  ;   expected_params = [expected_params,expected_params9]
  ; endif else if n_params lt 117 && gt 108 then begin ; lt gt
  ;   sdljf
  ; endif


  expected_vars=[expected_flags,expected_filenames,expected_masknames,expected_unfiltered,expected_params] ;names of all expected variables
  nvars=n_elements(expected_vars) ;number of expected variables
  for i=0,nvars-1 do if n_elements(scope_varfetch(expected_vars[i])) eq 0 then f_error,'variable '+expected_vars[i]+' not present in input file' ;verify input reading
  if peak_prof eq 0 and npixmin gt 1 then f_error,['peak_prof=0 (points), but minimum number of pixels npixmin>1','Please set npixmin=1 or peak_prof>0']




end
