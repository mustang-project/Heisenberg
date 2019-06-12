pro interactive_peak_find, peak_find_tui, $
  nsigma, npixmin, $ ; global variables
  loglevels, flux_weight, use_gas2, use_star2, $ ; control switches
  peakdir, starfile2short, gasfile2short, outputdir, galaxy,  $ ; filepaths
  peakidstar,logspacing_s,logrange_s, nlevels_s, nlinlevel_s, $ ; star levels variables
  peakidgas,logspacing_g,logrange_g, nlevels_g, nlinlevel_g, $  ; gas levels variables
  sensstar, offstar, $ ; star sensitivity variables
  sensgas, offgas, $   ; gas sensitivity variables
  ; ########################################################################
  starpeaks, gaspeaks ; output variables with peaks
  ; ########################################################################
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range



  run_stars = 1  ; variables to control running of peak finding on gas and stars
  run_gas = 1    ;
  break_test = 0
  counter = -1
  while(break_test eq 0) do begin ; start peak_finding
    empty_stars = 0
    empty_gas = 0


    if run_stars eq 1 then begin
      starpeaks = peak_find(peak_find_tui, $
        loglevels, peakidstar,logspacing_s,logrange_s, nlevels_s, $
        peakdir, starfile2short, $
        flux_weight, npixmin, use_star2, $ ; inputs (x is gas or stars as appropriate)
        nsigma, sensstar, offstar)

      if (size(starpeaks,/type) eq 7) then begin ; catch errors. 7 = string
        if starpeaks eq 'TOOFEWLINES' then begin
          empty_stars = 1
        endif

      endif

    endif
    if run_gas eq 1 then begin
      gaspeaks = peak_find(peak_find_tui, $
        loglevels, peakidgas,logspacing_g,logrange_g, nlevels_g, $
        peakdir, gasfile2short, $
        flux_weight, npixmin, use_gas2, $ ; inputs (x is gas or stars as appropriate)
        nsigma, sensgas, offgas)

      if (size(gaspeaks,/type) eq 7) then begin ; catch errors. 7 = string
        if gaspeaks eq 'TOOFEWLINES' then begin
          empty_gas = 1
        endif
      endif

    endif



    if peak_find_tui ne 1 then break_test = 1 else begin ; if requested, begin peak-identification TUI

      run_test = 0 ; test for run_stars and run_gas = 0
      while (run_test eq 0) do begin

        if run_stars eq 1 || run_gas eq 1 then begin ; don't reshow peaks if looping to ensure changed variables
          counter ++
          print, 'PeakID # ' +  strcompress(string(counter), /remove_all) + '. Current settings:'
          print, 'npixmin          ', npixmin
          print, 'nsigma           ', nsigma
          print, 'logrange_s       ', logrange_s
          print, 'logspacing_s     ', logspacing_s
          print, 'logrange_g       ', logrange_g
          print, 'logspacing_g     ', logspacing_g
          print, 'nlinlevel_s      ', nlinlevel_s
          print, 'nlinlevel_g      ', nlinlevel_g


          star_window_title = 'starmap_peakID_#' + strcompress(string(counter), /remove_all)
          gas_window_title = 'gasmap_peakID_#' + strcompress(string(counter), /remove_all)
          if empty_stars eq 1 then ds9_display_peaks, peakdir, starfile2short, star_window_title else  ds9_display_peaks, peakdir, starfile2short, star_window_title, starpeaks[*,0], starpeaks[*,1]
          if empty_gas eq 1 then ds9_display_peaks, peakdir, gasfile2short, gas_window_title else ds9_display_peaks, peakdir, gasfile2short, gas_window_title, gaspeaks[*,0], gaspeaks[*,1]
        endif


        if empty_stars eq 1 || empty_gas eq 1 then begin
          test = 'n' ; must redo the peak selection
          if (empty_stars eq 1 && empty_gas eq 1) then empty_maps = "star and gas maps" else $
          if (empty_stars eq 1) then  empty_maps = "star map" else $
          if (empty_gas eq 1) then  empty_maps = "gas map"
          print, "Too few peaks were found in the " + empty_maps + " . You must change at least one peak identification parameter"


        endif else begin
          print, "Are you satisfied with the identified peaks? Type 'Yes' or 'No' "
          test =''
          beep
          read, test
          test = strlowcase(test)
        endelse

        if test eq 'y' || test eq 'yes' then begin
          break_test = 1 ; to break out of main loop
          break ; break out of run testing loop
        endif else begin  ; add check for no or mistaken input...
          run_stars = 0
          run_gas = 0

          ; 1) global variables
          ; npixmin
          npixmin_old = npixmin
          npixmin = peak_var_read('npixmin', npixmin_old, /global_var)

          if npixmin ne npixmin_old then begin
            run_stars = 1
            run_gas = 1
          endif

          ; nsigma
          nsigma_old = nsigma
          nsigma = peak_var_read('nsigma', nsigma_old, /global_var)

          if nsigma ne nsigma_old then begin
            run_stars = 1 ; global variable so must rerun both stars and gas
            run_gas = 1
          endif

          if loglevels eq 1 then begin ; change log spacing parameters
            logrange_s_old = logrange_s
            logrange_s = peak_var_read('logrange_s', logrange_s_old)
            logspacing_s_old = logspacing_s
            logspacing_s = peak_var_read('logspacing_s', logspacing_s_old)
            if logrange_s ne logrange_s_old || logspacing_s ne logspacing_s_old then run_stars = 1

            logrange_g_old = logrange_g
            logrange_g = peak_var_read('logrange_g', logrange_g_old)
            logspacing_g_old = logspacing_g
            logspacing_g = peak_var_read('logspacing_g', logspacing_g_old)
            if logrange_g ne logrange_g_old || logspacing_g ne logspacing_g_old then run_gas = 1

          endif else begin ; change linear spacing parameters
            nlinlevel_s_old = nlinlevel_s
            nlinlevel_s = peak_var_read('nlinlevel_s', nlinlevel_s_old)
            if nlinlevel_s ne nlinlevel_s_old then run_stars = 1

            nlinlevel_g_old = nlinlevel_g
            nlinlevel_g = peak_var_read('nlinlevel_g', nlinlevel_g_old)
            if nlinlevel_g ne nlinlevel_g_old then run_gas = 1

          endelse

          nlevels_s = calc_levels(loglevels, logrange_s, logspacing_s, nlinlevels_s) ; calculate new stellar levels
          nlevels_g = calc_levels(loglevels, logrange_g, logspacing_g, nlinlevels_g) ; calculate new gas levels

        endelse
        if run_stars eq 1 || run_gas eq 1 then run_test = 1 else print, "No peakfinding variables are different to the variables of the previous peakfinding run. Please change at least one variable if you want to rerun."
      endwhile

    endelse
  endwhile


  if peak_find_tui eq 1 then begin
    ; test to see if the user whishes to repeat interactive peak finding in future steps
    final_break_test = 0
    while(final_break_test eq 0) do begin ; start peak_finding
      print, "Do you want to finalise peak finding parameters for future iterations? Type 'Yes' or 'No' "
      print, " 'Yes' = do not adjust peak finding in future iterations, 'No' = continue to adjust peak finding in future iterations "
      finaltest =''
      beep
      read, finaltest
      finaltest = strlowcase(finaltest)
      if finaltest eq 'y' || finaltest eq 'yes' then begin
        finalised = 1
        print, 'interactive peakfinding will be disabled in future iterations'
        final_break_test = 1

      endif else  if finaltest eq 'n' || finaltest eq 'no' then begin
        finalised = 0
        print, 'interactive peakfinding will continue to be enabled in future iterations'
        final_break_test = 1

      endif else begin
        print, "Incorect string entered. Type 'Yes' or 'No' "
      endelse
    endwhile
  endif





  ; write out final region files without displaying in a ds9 window
  ds9_display_peaks, peakdir, starfile2short, star_window_title, starpeaks[*,0], starpeaks[*,1], /no_disp
  ds9_display_peaks, peakdir, gasfile2short, gas_window_title, gaspeaks[*,0], gaspeaks[*,1], /no_disp


  ; write out final parameters
  peak_report_file = outputdir+galaxy+ '_interactive_peak_ID_report.dat'
  openw, peak_lun, peak_report_file, width = 150, /get_lun
  printf, peak_lun, '# Interactive Peak Identification report. Final selected parameters'
  printf, peak_lun, 'npixmin          ', npixmin
  printf, peak_lun, 'nsigma           ', nsigma
  printf, peak_lun, 'logrange_s       ', logrange_s
  printf, peak_lun, 'logspacing_s     ', logspacing_s
  printf, peak_lun, 'logrange_g       ', logrange_g
  printf, peak_lun, 'logspacing_g     ', logspacing_g
  printf, peak_lun, 'nlinlevel_s      ', nlinlevel_s
  printf, peak_lun, 'nlinlevel_g      ', nlinlevel_g
  if n_elements(finalised) eq 0 then finalised = 0 ; if not already set then = 0
  printf, peak_lun, 'finalised        ', finalised
  free_lun, peak_lun




end
