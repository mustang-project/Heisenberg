;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                            ;
;              FIT KL14 PRINCIPLE (v0.2) TO OBSERVED GALAXY MAPS             ;
; start environment with >> idl kl14 -arg [full/absolute path of input file] ;
;                                                                            ;
;                      BEGIN MAIN ROUTINE tuningfork.pro                     ;
;                                                                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CHECK INPUT FILE VALIDITY;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

startrun=0
while startrun eq 0 && inputfile ne '0' do begin ;as long as the run is not greenlit and an input file is specified, check validity and request new file if needed
    filetest=file_search(inputfile,count=ctfile) ;test if input file exists
    if ctfile eq 0 then begin ;if it does not, then
        startrun1=0 ;prevent execution, request new input, and force another test of the new input
        print,' error: '+inputfile+' does not exist'
        read,' please specify the full/absolute path of the input file (enter 0 to stop autorun): ',inputfile
    endif else startrun1=1 ;if it does, then greenlight for file check
    dirtest=file_search(inputfile,count=ctdir,/test_directory) ;test if the input file is a directory
    if ctdir ne 0 then begin ;if it is, then
        startrun2=0 ;prevent execution, request new input, and force another test of the new input
        print,' error: '+inputfile+' is a directory'
        read,' please specify the full/absolute path of the input file (enter 0 to stop autorun): ',inputfile
    endif else startrun2=1 ;if it does, then greenlight for dir check
    startrun=startrun1*startrun2 ;if both tests are passed we can run
endwhile
if inputfile eq '0' then stop ;stop if no input file is provided


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE DIRECTORY STRUCTURE AND START RUN;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

inputdir=file_dirname(inputfile)+path_sep() ;input directory
runname=file_basename(inputfile)+'_run' ;main run directory name
rundir=inputdir+runname+path_sep() ;main run directory full path
dummy=file_search(rundir,count=ct) ;check if rundir exists
if ct eq 0 then spawn,'mkdir '+rundir ;if not, create it
if ct ne 0 then print,' WARNING: run directory already exists, some or all files may be overwritten'

griddir=rundir+'regrid'+path_sep() ;directory for saving regridded maps
dummy=file_search(griddir,count=ct) ;check if griddir exists
if ct eq 0 then spawn,'mkdir '+griddir ;if not, create it
peakdir=rundir+'peakid'+path_sep() ;directory for saving peak ID maps
dummy=file_search(peakdir,count=ct) ;check if peakdir exists
if ct eq 0 then spawn,'mkdir '+peakdir ;if not, create it
figdir=rundir+'figures'+path_sep() ;directory for saving figures
dummy=file_search(figdir,count=ct) ;check if figdir exists
if ct eq 0 then spawn,'mkdir '+figdir ;if not, create it
outputdir=rundir+'output'+path_sep() ;directory for saving output
dummy=file_search(outputdir,count=ct) ;check if outputdir exists
if ct eq 0 then spawn,'mkdir '+outputdir ;if not, create it
arrdir=rundir+'arrays'+path_sep() ;directory for saving temporary arrays
dummy=file_search(arrdir,count=ct) ;check if arrdir exists
if ct eq 0 then spawn,'mkdir '+arrdir ;if not, create it
maskeddir=rundir+'masked'+path_sep() ;directory for saving temporary arrays
dummy=file_search(maskeddir,count=ct) ;check if maskeddir exists
if ct eq 0 then spawn,'mkdir '+maskeddir ;if not, create it

starttime=systime(1)
logfile='logfile.txt'
journal,outputdir+logfile ;start keeping journal of output -- note that galaxy name is only specified at input file read-in, so journal file is renamed at end of run to include galaxy name
print,' ==> starting KL14 tuning fork analysis, date/time is ',systime(0)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;READ INPUT FILE AND VERIFY;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

openr,lun,inputfile,/get_lun ;open the input file
line='' ;initialise line
print,' ==> reading '+inputfile+' with the following variables:'
while ~eof(lun) do begin ;read lines until end string is reached
    readf,lun,line ;read line containing variable
    linearr=strsplit(line,' ',/extract) ;split line into space-separated parts
    if n_elements(linearr) gt 1 && strmid(linearr(0),0,1) ne '#' then begin ;if line is not empty or comment
        value=linearr(1)
        is_num=strnumber(value,value)
        (scope_varfetch(linearr(0),/enter))=value ;assign variable name and value
        spacelen=20-strlen(linearr(0))
        space=''
        for i=0,spacelen-1 do space+=' '
        print,'     '+linearr(0)+space+linearr(1)
    endif
endwhile
close,lun ;close the logical unit
free_lun,lun ;make the logical unit available again

expected_flags1=['mask_images','regrid','smoothen','sensitivity','id_peaks','calc_ap_flux','generate_plot','get_distances','calc_obs','calc_fit','diffuse_frac','derive_phys','write_output','cleanup','autoexit'] ;variable names of expected flags (1)
expected_flags2=['use_star2','use_gas2','use_star3'] ;variable names of expected flags (2)
expected_flags3=['mstar_ext','mstar_int','mgas_ext','mgas_int','mstar_ext2','mstar_int2','mgas_ext2','mgas_int2','mstar_ext3','mstar_int3','convert_masks','cut_radius'] ;variable names of expected flags (3)
expected_flags4=['set_centre','tophat','loglevels','peak_find_tui','flux_weight','calc_ap_area','tstar_incl','peak_prof','map_units','use_X11'] ;variable names of expected flags (4)
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
expected_params=[expected_params1,expected_params2,expected_params3,expected_params4,expected_params5,expected_params6,expected_params7,expected_params8] ;variable names of expected input parameters (all)
expected_vars=[expected_flags,expected_filenames,expected_masknames,expected_params] ;names of all expected variables
nvars=n_elements(expected_vars) ;number of expected variables
for i=0,nvars-1 do if n_elements(scope_varfetch(expected_vars(i))) eq 0 then f_error,'variable '+expected_vars(i)+' not present in input file' ;verify input reading
if peak_prof eq 0 and npixmin gt 1 then f_error,['peak_prof=0 (points), but minimum number of pixels npixmin>1','Please set npixmin=1 or peak_prof>0']

if use_X11 then begin ;set window properties
    window_plot='x'
    device,retain=2
    ; set color
    device,true_color=24
    device,decomposed=0
    window,xsize=600,ysize=840,title='IDL graphics',colors=100
    loadct,12
endif else window_plot='z' ;prevents window creation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DERIVE AND SET ADDITIONAL PROGRAM PARAMETERS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

inclination=inclination*!dtor ;inclination in radians
posangle=posangle*!dtor-!pi/2. ;position angle in radians and rotated by pi/2
apertures=lapmin*(lapmax/lapmin)^(dindgen(naperture)/(naperture-1)) ;aperture size (diameter) array in pc
res=f_string(apertures,0) ;aperture size string array
if loglevels then begin ;set number of logarithmic contour levels for peak identification
    nlevels_s=logrange_s/logspacing_s+1 ;number of stellar contour levels
    nlevels_g=logrange_g/logspacing_g+1 ;number of gas contour levels
endif else begin ;set number of linear contour levels for peak identification
    nlevels_s=nlinlevels_s ;number of stellar contour levels
    nlevels_g=nlinlevels_s ;number of gas contour levels
endelse
if tstar_incl eq 0 then tovermaxi=tgasmaxi else tovermaxi=min([tstariso,tgasmaxi])

mainseed=systime(1) ;source seed for random number generation
nseed=100*nmc*naperture ;number of random seeds needed
seedarr=floor(randomu(mainseed,nseed)*nseed) ;random seed array
iseed=long(0) ;index of next seed to be used


;;;;;;;;;;;
;READ DATA;
;;;;;;;;;;;

beamtest=0.
if regrid then resdir=datadir else resdir=griddir
print,' ==> reading maps from '+resdir
starfiletot=resdir+starfile ;Halpha map (total)
if mask_images then begin
    if mstar_ext then mstar_ext_path = maskdir + path_sep() + star_ext_mask
    if mstar_int then mstar_int_path = maskdir + path_sep() + star_int_mask
    masked_path_star = maskeddir + starfile
endif
mask_tool, starfiletot $ ;image variables
    , ds9_positive_path = mstar_ext_path, ds9_negative_path = mstar_int_path $ ;ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
    , masked_image_path = masked_path_star $ ;path to write masked image to
    , image_masked = starmap $ ;overwrite original image array with masked one
    , header_output = starmaphdr $
    , convert = convert_masks, /run_without_masks

cdeltstar=abs(sxpar(starmaphdr,'CDELT1',count=ct))
if ct eq 0 then cdeltstar = abs(sxpar(starmaphdr,'CD1_1'))  ;try alternative CDELT
beamfits=sqrt(sxpar(starmaphdr, 'BMAJ', count=ct)*sxpar(starmaphdr, 'BMIN', count=ct))
beamtest=max([beamtest,beamfits])
databeam=beamtest
if use_star2 then begin
    starfiletot2=resdir+starfile2 ;Halpha peak map
    if mask_images then begin
        if mstar_ext2 then mstar_ext_path2 = maskdir + path_sep() + star_ext_mask2
        if mstar_int2 then mstar_int_path2 = maskdir + path_sep() + star_int_mask2
        masked_path_star2 = maskeddir + starfile2
    endif

    mask_tool, starfiletot2 $ ;image variables
        , ds9_positive_path = mstar_ext_path2, ds9_negative_path = mstar_int_path2 $ ;ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
        , masked_image_path = masked_path_star2 $ ;path to write masked image to
        , image_masked = starmap2 $ ;overwrite original image array with masked one
        , header_output = starmaphdr2 $
        , convert = convert_masks, /run_without_masks

        cdeltstar2=abs(sxpar(starmaphdr2,'CDELT1',count=ct))
        if ct eq 0 then cdeltstar2 = abs(sxpar(starmaphdr2,'CD1_1'))  ;try alternative CDELT
        beamfits=sqrt(sxpar(starmaphdr2, 'BMAJ', count=ct)*sxpar(starmaphdr2, 'BMIN', count=ct))
        beamtest=max([beamtest,beamfits])
endif else begin
    starfile2=starfile
    starmap2=starmap
endelse

if use_star3 then begin
    starfiletot3=resdir+starfile3 ;FUV map
    if mask_images then begin
        if mstar_ext3 then mstar_ext_path3 = maskdir + path_sep() + star_ext_mask3
        if mstar_int3 then mstar_int_path3 = maskdir + path_sep() + star_int_mask3
        masked_path_star3 = maskeddir + starfile3
    endif

    mask_tool, starfiletot3 $ ;image variables
        , ds9_positive_path = mstar_ext_path3, ds9_negative_path = mstar_int_path3 $ ;ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
        , masked_image_path = masked_path_star3 $ ;path to write masked image to
        , image_masked = starmap3 $ ;overwrite original image array with masked one
        , header_output = starmaphdr3 $
        , convert = convert_masks, /run_without_masks

        cdeltstar3=abs(sxpar(starmaphdr3,'CDELT1',count=ct))
        if ct eq 0 then cdeltstar3 = abs(sxpar(starmaphdr3,'CD1_1'))  ;try alternative CDELT
        beamfits=sqrt(sxpar(starmaphdr3, 'BMAJ', count=ct)*sxpar(starmaphdr3, 'BMIN', count=ct))
        beamtest=max([beamtest,beamfits])
endif

gasfiletot=resdir+gasfile ;moment zero CO map (velocity mask for better flux estimate)
if mask_images then begin
    if mgas_ext then mgas_ext_path = maskdir + path_sep() + gas_ext_mask
    if mgas_int then mgas_int_path = maskdir + path_sep() + gas_int_mask
    masked_path_gas = maskeddir + gasfile
endif

mask_tool, gasfiletot $ ;image variables
    , ds9_positive_path = mgas_ext_path, ds9_negative_path = mgas_int_path $ ;ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
    , masked_image_path = masked_path_gas $ ;path to write masked image to
    , image_masked = gasmap $ ;overwrite original image array with masked one
    , header_output = gasmaphdr $
    , convert = convert_masks, /run_without_masks

cdeltgas=abs(sxpar(gasmaphdr,'CDELT1',count=ct))
if ct eq 0 then cdeltgas = abs(sxpar(gasmaphdr,'CD1_1'))  ;try alternative CDELT
beamfits=sqrt(sxpar(gasmaphdr, 'BMAJ', count=ct)*sxpar(gasmaphdr, 'BMIN', count=ct))
beamtest=max([beamtest,beamfits])
databeam=max([databeam,beamfits])
if use_gas2 then begin
    gasfiletot2=resdir+gasfile2 ;moment zero CO map (sigma mask for better peak ID)
    if mask_images then begin
        if mgas_ext2 then mgas_ext_path2 = maskdir + path_sep() + gas_ext_mask2
        if mgas_int2 then mgas_int_path2 = maskdir + path_sep() + gas_int_mask2
        masked_path_gas2 = maskeddir + gasfile2
    endif

    mask_tool, gasfiletot2 $ ;image variables
        , ds9_positive_path = mgas_ext_path2, ds9_negative_path = mgas_int_path2 $ ;ds9 region file keywords ; positive = allowed regions, negative = not allowed regions
        , masked_image_path = masked_path_gas2 $ ;path to write masked image to
        , image_masked = gasmap2 $ ;overwrite original image array with masked one
        , header_output = gasmaphdr2 $
        , convert = convert_masks, /run_without_masks

        cdeltgas2=abs(sxpar(gasmaphdr2,'CDELT1',count=ct))
        if ct eq 0 then cdeltgas2 = abs(sxpar(gasmaphdr2,'CD1_1'))  ;try alternative CDELT
        beamfits=sqrt(sxpar(gasmaphdr2, 'BMAJ', count=ct)*sxpar(gasmaphdr2, 'BMIN', count=ct))
        beamtest=max([beamtest,beamfits])
endif else begin
    gasfile2=gasfile
    gasmap2=gasmap
endelse

if ~cut_radius then begin
    minradius=0.
    maxradius=huge
endif
mapdim_orig=size(starmap)
nx_orig=mapdim_orig(1) ;number of x pixels in original stellar map
ny_orig=mapdim_orig(2) ; number of y pixels in original stellar map
if set_centre then begin
    centrefracx=centrex/long(nx_orig-1) ;central pixel x-coordinate measured from map edge
    centrefracy=centrey/long(ny_orig-1) ;central pixel y-coordinate measured from map edge
endif else begin
    centrefracx=0.5
    centrefracy=0.5
endelse
centre_orig=[centrefracx*(nx_orig-1),centrefracy*(ny_orig-1)] ;pixel coordinates of the galaxy centre in the original stellar map
beammaxpc=distance*beamtest*!dtor/sqrt(cos(inclination)) ;largest beam size in pc -- assumes small angles, i.e. tan(x)~x
beamaperture=min(where(abs(alog10(apertures/beammaxpc)) eq min(abs(alog10(apertures/beammaxpc)))))
peak_res=max([peak_res,beamaperture]) ;ensure that the smallest aperture size is the aperture closest to the beam size
fitap=fix(peak_res)+indgen(max_res-peak_res+1) ;aperture sizes used in fitting the KL14 principle model
lap_min=apertures(peak_res) ;size of smallest aperture


;;;;;;;;;;;;;
;REGRID MAPS;
;;;;;;;;;;;;;

if regrid then begin
    print,' ==> regridding maps'
    if mask_images then regriddir = maskeddir else regriddir = resdir ;if images are masked, get them from corresponding directory
    if diffuse_frac && beamtest gt databeam then f_error,'Maximum beam of all images (inlcuding starmap2, starmap3 and gasmap2) is larger than the maximum beam of the two main images (starmap and gasmap). This will affect the resolution avaliable for Fourier filtering in the calculation of fcl and fgmc'
    if diffuse_frac then targetres=beammaxpc<apertures(peak_res) else targetres=apertures(peak_res) ; if calculating diffuse fraction then preserve at least worst beam resolution to maximise fourier resolution, else cut resolution to smallest aperture
    starmap=smoothfits(distance,inclination,targetres,starfile,regriddir,griddir,max_sample,tophat=0,regrid=regrid)
    starfileshort=strmid(starfile,0,strpos(starfile,'.fits')) ;short name without extension
    starfile2short=strmid(starfile2,0,strpos(starfile2,'.fits')) ;short name without extension
    starfiletot=griddir+starfile ;Halpha map (total)
    starmap=readfits(starfiletot,hdr,/silent) ;read Halpha map
    starmaphdr=hdr ;header of Halpha map
    if use_star2 then begin
        starmap2=smoothfits(distance,inclination,targetres,starfile2,regriddir,griddir,max_sample,tophat=0,regrid=regrid)
        starfiletot2=griddir+starfile2 ;Halpha peak map
        starmap2=readfits(starfiletot2,hdr,/silent) ;read Halpha peak map
        starmaphdr2=hdr ;header of Halpha peak map
    endif
    if use_star3 then begin
        starmap3=smoothfits(distance,inclination,targetres,starfile3,regriddir,griddir,max_sample,tophat=0,regrid=regrid)
        starfile3short=strmid(starfile3,0,strpos(starfile3,'.fits')) ;short name without extension
        starfiletot3=griddir+starfile3 ;FUV map
        starmap3=readfits(starfiletot3,hdr,/silent) ;read FUV map
        starmaphdr3=hdr ;header of FUV map
    endif
    gasmap=smoothfits(distance,inclination,targetres,gasfile,regriddir,griddir,max_sample,tophat=0,regrid=regrid)
    gasfileshort=strmid(gasfile,0,strpos(gasfile,'.fits')) ;short name without extension
    gasfile2short=strmid(gasfile2,0,strpos(gasfile2,'.fits')) ;short name without extension
    gasfiletot=griddir+gasfile ;moment zero CO map (velocity mask for better flux estimate)
    gasmap=readfits(gasfiletot,hdr,/silent) ;read moment zero CO map
    gasmaphdr=hdr ;moment zero CO header
    if use_gas2 then begin
        gasmap2=smoothfits(distance,inclination,targetres,gasfile2,regriddir,griddir,max_sample,tophat=0,regrid=regrid)
        gasfiletot2=griddir+gasfile2 ;moment zero CO map (signal mask for better peak ID)
        gasmap2=readfits(gasfiletot2,hdr,/silent) ;read moment zero CO peak map
        gasmaphdr2=hdr ;moment zero CO peak header
    endif else begin
        gasmap2=gasmap
    endelse
    mapdim=size(starmap)
    nx=mapdim(1) ;number of x pixels in regridded stellar map
    ny=mapdim(2) ; number of y pixels in regridded stellar map
    centre=[centrefracx*(nx-1),centrefracy*(ny-1)] ;pixel coordinates of the galaxy centre in the regridded stellar map
endif else centre=centre_orig

cdelt=abs(sxpar(starmaphdr,'CDELT1',count=ct))
if ct eq 0 then cdelt = abs(sxpar(starmaphdr,'CD1_1'))  ;try alternative CDELT
pixtopc=distance*tan(!dtor*cdelt)/sqrt(cos(inclination)) ;pixel size in pc -- assumes small angles, i.e. tan(x)~x
if sqrt(2.)*pixtopc gt apertures(naperture-2) then f_error,'less than 2 aperture sizes exceed inclination-corrected pixel diagonal'
if sqrt(2.)*pixtopc gt apertures(peak_res) then begin
    print, ' WARNING: minimum aperture size ('+f_string(lap_min,0)+' pc) is smaller than inclination-corrected pixel diagonal ('+f_string(sqrt(2.)*pixtopc,1)+' pc)'
    peak_res=min(where(apertures gt sqrt(2.)*pixtopc)) ;ensure that the smallest aperture size exceeds the pixel size
    fitap=fix(peak_res)+indgen(max_res-peak_res+1) ;aperture sizes used in fitting the KL14 principle model
    lap_min=apertures(peak_res) ;size of smallest aperture
    print, ' WARNING: setting minimum aperture size to '+f_string(lap_min,0)+' pc to exceed pixel diagonal'
endif
pctopix=1./pixtopc ;pc in number of pixels
convstar=10.^convstar*(cdelt/cdeltstar)^2. ;change pixel conversion factor to physical units to linear scale and account for regridding
convstar_err=convstar*convstar_rerr
convgas=10.^convgas*(cdelt/cdeltgas)^2. ;change pixel conversion factor to physical units to linear scale and account for regridding
convgas_err=convgas*convgas_rerr
if use_star3 then begin
    convstar3=10.^convstar3*(cdelt/cdeltstar3)^2. ;change pixel conversion factor to physical units to linear scale and account for regridding
    convstar3_err=convstar3*convstar3_rerr
endif

tstariso_rerrmin=tstariso_errmin/tstariso ;Relative downward standard error on tstariso
tstariso_rerrmax=tstariso_errmax/tstariso ;Relative upward standard error on tstariso


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SYNCHRONISE MASKS FOR ALL IMAGES AND CALCULATE TOTAL GAS MASSES AND/OR STAR FORMATION RATES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if regrid then begin
    print,' ==> syncronising masks of maps from '+griddir
    masksimple = 1 ; propogate masks using simple pixel by pixel comparison
    if astrometry_equal(starmap, starmaphdr, gasmap, gasmaphdr, astr_tolerance) then begin ;check astrometry is equal
        if use_star2 then if astrometry_equal(starmap, starmaphdr, starmap2, starmaphdr2, astr_tolerance) ne 1 then masksimple = 0
        if use_star3 then if astrometry_equal(starmap, starmaphdr, starmap3, starmaphdr3, astr_tolerance) ne 1 then masksimple = 0
        if use_gas2 then if astrometry_equal(starmap, starmaphdr, gasmap2, gasmaphdr2, astr_tolerance) ne 1 then masksimple = 0
    endif else masksimple = 0

    if masksimple then begin ;create array of total mask
        mask_arr = starmap ;create array to hold masked pixels
        mask_arr[*] = 1.0 ;1.0 = not masked (i.e. allowed through)

        nan_list = where(finite(starmap, /nan), nancount) ;find masked pixels in starmap
        if nancount gt 0 then mask_arr[nan_list] = 0.0
        nan_list = where(finite(gasmap, /nan), nancount)  ;find masked pixels in gasmap
        if nancount gt 0 then mask_arr[nan_list] = 0.0
        if use_star2 then begin
            nan_list = where(finite(starmap2, /nan), nancount)  ;find masked pixels in starmap2
            if nancount gt 0 then mask_arr[nan_list] = 0.0
        endif
        if use_star3 then begin
            nan_list = where(finite(starmap3, /nan), nancount)  ;find masked pixels in starmap3
            if nancount gt 0 then mask_arr[nan_list] = 0.0
        endif
        if use_gas2 then begin
            nan_list = where(finite(gasmap2, /nan), nancount)  ;find masked pixels in gasmap2
            if nancount gt 0 then mask_arr[nan_list] = 0.0
        endif

        ;make radial cut
        pixx=dblarr(nx,ny) ;array containing map x pixel indices
        pixy=dblarr(nx,ny) ;array containing map y pixel indices
        for i=0,nx-1 do pixx(i,*)=i
        for i=0,ny-1 do pixy(*,i)=i
        pixdistxpix=pixx-centre(0) ;x distance of each pixel from centre in number of pixels
        pixdistypix=pixy-centre(1) ;y distance of each pixel from centre in number of pixels
        pixdistxpc=(cos(-posangle)*pixdistxpix-sin(-posangle)*pixdistypix)*distance*tan(!dtor*cdelt) ;x distance of each pixel from centre in pc
        pixdistypc=(sin(-posangle)*pixdistxpix+cos(-posangle)*pixdistypix)*distance*tan(!dtor*cdelt)/cos(inclination) ;y distance of each pixel from centre in pc
        pixradius=sqrt(pixdistxpc^2.+pixdistypc^2.) ;distance of each pixel from centre in pc
        inclpix=where(pixradius ge minradius and pixradius le maxradius, ninclpix, complement=exclpix, ncomplement=nexclpix) ;find included pixels, excluded pixels and number of both for the radial cut
        pixmeanradius=mean(pixradius(inclpix)) ;mean radius of included area
        if nexclpix gt 0 then mask_arr[exclpix] = 0.0 ;mask excluded area
        totalarea=total(mask_arr)*pixtopc^2. ;total included area
        starfluxtotal=total(starmap*mask_arr,/nan)*pixtopc^2. ;total SF flux in included area
        gasfluxtotal=total(gasmap*mask_arr,/nan)*pixtopc^2. ;total gas flux in included area

        ;propogate masks
        nan_list = where(mask_arr eq 0.0, nancount) ;final mask list
        if nancount gt 0 then starmap[nan_list] = !values.f_nan ;propogate masks
        if nancount gt 0 then gasmap[nan_list] = !values.f_nan ;propogate masks
        writefits, starfiletot, starmap, starmaphdr ;write out masked starmap
        writefits, gasfiletot, gasmap, gasmaphdr ;write out masked gasmap
        if use_star2 then begin
            if nancount gt 0 then starmap2[nan_list] = !values.f_nan
            writefits, starfiletot2, starmap2, starmaphdr2
        endif
        if use_gas2 then begin
            if nancount gt 0 then gasmap2[nan_list] = !values.f_nan
            writefits, gasfiletot2, gasmap2, gasmaphdr2 ;write out masked gasmap
        endif
        if use_star3 then begin
            if nancount gt 0 then starmap3[nan_list] = !values.f_nan
            writefits, starfiletot3, starmap3, starmaphdr3 ;write out masked starmap
        endif

        ;write out mask file
        maskfile = 'totalmask.fits' ;filename for the mask
        maskfiletot = maskeddir + maskfile
        maskmaphdr = starmaphdr ;use the starmap header
        writefits, maskfiletot, mask_arr, maskmaphdr ;save a fits image of synchronised mask to maskeddir, which is needed to create smoothmask and work out the area of the mask

    endif else f_error,['for masks to be propogated between images, all images must share the same astrometry','set regrid=1 in the parameter file to regrid images']
endif
if map_units eq 1 then begin
    sfr_galaxy=total(starmap,/nan)*convstar ;total star formation rate in SF map
    sfr_galaxy_err=convstar_rerr*sfr_galaxy ;standard error on SFR
    mgas_galaxy=total(gasmap,/nan)*convgas ;total gas mass in gas map
    mgas_galaxy_err=convgas_rerr*mgas_galaxy ;standard error on gas mass
    tdeplmax=mgas_galaxy/sfr_galaxy*(1.+sqrt((mgas_galaxy_err/mgas_galaxy)^2.+(sfr_galaxy_err/sfr_galaxy)^2.))/1.d9 ;gas depletion time + 1sigma
    if tgasmaxi gt tdeplmax*1.d3 then tgasmaxi=tdeplmax*1.d3 ;tgas cannot exceed the gas depletion time
endif
if map_units eq 2 then begin
    mgas_galaxy1=total(starmap,/nan)*convstar ;total gas mass in "SF" map
    mgas_galaxy2=total(gasmap,/nan)*convgas ;total gas mass in gas map
    mgas_galaxy1_err=convstar_rerr*mgas_galaxy1 ;standard error on gas mass 1
    mgas_galaxy2_err=convgas_rerr*mgas_galaxy2 ;standard error on gas mass 2
endif
if map_units eq 3 then begin
    sfr_galaxy1=total(starmap,/nan)*convstar ;total star formation rate in SF map
    sfr_galaxy2=total(gasmap,/nan)*convgas ;total star formation rate in "gas" map
    sfr_galaxy1_err=convstar_rerr*sfr_galaxy1 ;standard error on SFR 1
    sfr_galaxy2_err=convgas_rerr*sfr_galaxy2 ;standard error on SFR 2
endif


;;;;;;;;;;;;;;;
;SMOOTHEN MAPS;
;;;;;;;;;;;;;;;

if smoothen then begin
    print,' ==> smoothing data'
    peakidstar=smoothfits(distance,inclination,apertures(peak_res),starfile2,griddir,peakdir,max_sample,tophat=0) ;generate Gaussian psf-smoothened Ha map for peak ID
    peakidgas=smoothfits(distance,inclination,apertures(peak_res),gasfile2,griddir,peakdir,max_sample,tophat=0) ;generate Gaussian psf-smoothened gas map for peak ID
    smoothstar=dblarr([naperture,size(peakidstar,/dimensions)]) ;create Ha map array for different resolutions
    if use_star3 then smoothstar3=dblarr([naperture,size(peakidstar,/dimensions)]) ;create FUV map array for different resolutions
    smoothgas=dblarr([naperture,size(peakidstar,/dimensions)]) ;create gas map array for different resolutions
    smoothmask=dblarr([naperture,size(peakidstar,/dimensions)]) ;create total mask array for different resolutions
    for i=0,naperture-1 do begin ;generate smoothened maps
        smoothstartemp=smoothfits(distance,inclination,apertures(i),starfile,griddir,rundir+'res_'+res(i)+'pc'+path_sep(),max_sample,tophat=tophat) ;get smoothened Ha map
        smoothgastemp=smoothfits(distance,inclination,apertures(i),gasfile,griddir,rundir+'res_'+res(i)+'pc'+path_sep(),max_sample,tophat=tophat) ;get smoothened gas map
        smoothmasktemp=smoothfits(distance,inclination,apertures(i),maskfile,maskeddir,rundir+'res_'+res(i)+'pc'+path_sep(),max_sample,tophat=tophat) ;get smoothened mask
        smoothstar(i,*,*)=smoothstartemp
        smoothgas(i,*,*)=smoothgastemp
        smoothmask(i,*,*)=smoothmasktemp
        if use_star3 then begin
            smoothstar3temp=smoothfits(distance,inclination,apertures(i),starfile3,griddir,rundir+'res_'+res(i)+'pc'+path_sep(),max_sample,tophat=tophat) ;get smoothened FUV map
            smoothstar3(i,*,*)=smoothstar3temp
        endif
    endfor
endif else begin
    print,' ==> reading smoothened maps'
    peakidstar=readfits(peakdir+starfile2,hdr,/silent)
    peakidgas=readfits(peakdir+gasfile2,hdr,/silent)
    smoothstar=dblarr([naperture,size(peakidstar,/dimensions)]) ;create Ha map array for different resolutions
    if use_star3 then smoothstar3=dblarr([naperture,size(peakidstar,/dimensions)]) ;create FUV map array for different resolutions
    smoothgas=dblarr([naperture,size(peakidstar,/dimensions)]) ;create gas map array for different resolutions
    for i=0,naperture-1 do begin ;read smoothened maps
        ;check if target directory exists
        dir=rundir+'res_'+res(i)+'pc'+path_sep()
        dummy=file_search(dir,count=ct)
        if ct eq 0 then f_error,['data directory not present for '+res(i)+' pc resolution','first run with smoothen keyword enabled']
        smoothstar(i,*,*)=readfits(rundir+'res_'+res(i)+'pc'+path_sep()+starfile,hdr,/silent)
        if use_star3 then smoothstar3(i,*,*)=readfits(rundir+'res_'+res(i)+'pc'+path_sep()+starfile3,hdr,/silent)
        smoothgas(i,*,*)=readfits(rundir+'res_'+res(i)+'pc'+path_sep()+gasfile,hdr,/silent)
    endfor
endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DETERMINE SENSITIVITY LIMITS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if sensitivity then begin
    if use_stds then begin
        sensstar=std_star
        offstar=0.0d0
        sensstar3=std_star3
        offstar3=0.0d0
        sensgas=std_gas
        offgas=0.0d0
    endif else begin
        print,' ==> determining sensitivity limits'
        sensstarlist=reform(smoothstar(peak_res,*,*),n_elements(smoothstar(peak_res,*,*)))
        nans=where(finite(sensstarlist,/nan),ct)
        if ct ne 0 then remove,nans,sensstarlist
        sensstarmin=min(sensstarlist)
        sensstarmed=median(sensstarlist)
        if sensstarmin lt 0. then use=where(abs(sensstarlist) le abs(sensstarmin)) else use=[-1]
        if sensstarmin ge 0. || n_elements(use) lt nbins then use=where(abs(sensstarlist) le sensstarmed)
        if sensstarmin ge 0. && sensstarmed eq 0. then use=where(abs(sensstarlist) le min(sensstarlist(where(sensstarlist gt 0.)))) ;in practice should only happen with near-perfect S/N
        disp=sqrt(mean(sensstarlist(use)^2.)-mean(sensstarlist(use))^2.)
        binwidth=disp/nbins
        set_plot,window_plot
        histoplot,sensstarlist(use),histdata=hist,locations=bins,binsize=binwidth
        set_plot,'x'
        binmid=bins+0.5*binwidth
        fit=gaussfit(binmid,hist,vars,nterms=3)
        offstar=vars(1)
        sensstar=vars(2)

        if use_star3 then begin
            sensstarlist3=reform(smoothstar3(peak_res,*,*),n_elements(smoothstar3(peak_res,*,*)))
            nans=where(finite(sensstarlist3,/nan),ct)
            if ct ne 0 then remove,nans,sensstarlist3
            sensstarmin3=min(sensstarlist3)
            sensstarmed3=median(sensstarlist3)
            if sensstarmin3 lt 0. then use=where(abs(sensstarlist3) le abs(sensstarmin3)) else use=[-1]
            if sensstarmin3 ge 0. || n_elements(use) lt nbins then use=where(abs(sensstarlist3) le sensstarmed3)
            if sensstarmin3 ge 0. && sensstarmed3 eq 0. then use=where(abs(sensstarlist3) le min(sensstarlist3(where(sensstarlist3 gt 0.)))) ;in practice should only happen with near-perfect S/N
            disp=sqrt(mean(sensstarlist3(use)^2.)-mean(sensstarlist3(use))^2.)
            binwidth=disp/nbins
            set_plot,window_plot
            histoplot,sensstarlist3(use),histdata=hist,locations=bins,binsize=binwidth
            set_plot,'x'
            binmid=bins+0.5*binwidth
            fit=gaussfit(binmid,hist,vars,nterms=3)
            offstar3=vars(1)
            sensstar3=vars(2)
        endif

        sensgaslist=reform(smoothgas(peak_res,*,*),n_elements(smoothgas(peak_res,*,*)))
        nans=where(finite(sensgaslist,/nan),ct)
        if ct ne 0 then remove,nans,sensgaslist
        sensgasmin=min(sensgaslist)
        sensgasmed=median(sensgaslist)
        if sensgasmin lt 0. then use=where(abs(sensgaslist) le abs(sensgasmin)) else use=[-1]
        if sensgasmin ge 0. || n_elements(use) lt nbins then use=where(abs(sensgaslist) le sensgasmed)
        if sensgasmin ge 0. && sensgasmed eq 0. then use=where(abs(sensgaslist) le min(sensgaslist(where(sensgaslist gt 0.)))) ;in practice should only happen with near-perfect S/N
        disp=sqrt(mean(sensgaslist(use)^2.)-mean(sensgaslist(use))^2.)
        binwidth=disp/nbins
        set_plot,window_plot
        histoplot,sensgaslist(use),histdata=hist,locations=bins,binsize=binwidth
        set_plot,'x'
        binmid=bins+0.5*binwidth
        fit=gaussfit(binmid,hist,vars,nterms=3)
        offgas=vars(1)
        sensgas=vars(2)
    endelse
endif


;;;;;;;;;;;;;;;;
;IDENTIFY PEAKS;
;;;;;;;;;;;;;;;;

if id_peaks then begin
    print,' ==> identifying peaks'
    interactive_peak_find, peak_find_tui, $
      nsigma, npixmin, $ ; global variables
      loglevels, flux_weight, use_gas2, use_star2, $ ; control switches
      peakdir, starfile2short, gasfile2short, outputdir, galaxy,  $ ; filepaths
      peakidstar,logspacing_s,logrange_s, nlevels_s, nlinlevel_s, $ ; star levels variables
      peakidgas,logspacing_g,logrange_g, nlevels_g, nlinlevel_g, $  ; gas levels variables
      sensstar, offstar, $ ; star sensitivity variables
      sensgas, offgas, $   ; gas sensitivity variables
      ; ########################################################################
      starpeaks, gaspeaks ; output variables with peaks

    peaks=[starpeaks,gaspeaks]
    sz=size(peaks)
    sz_s=size(starpeaks)
    sz_g=size(gaspeaks)
    npeaks=sz(1) ;total number of all peaks
    nstarpeaks=sz_s(1) ;total number of stellar peaks
    ngaspeaks=sz_g(1) ;total number of gas peaks

    includepeak=mask_arr(peaks(*,0),peaks(*,1)) eq 1 ;set to 1 when peak is in included pixel, otherwise set to 0
    inclpeaks=where(includepeak ne 0,nincludepeak) ;include peak? only if within specified radius interval
    inclstar=where(includepeak(0:nstarpeaks-1) ne 0,nincludepeak_star) ;included stellar peaks
    inclgas=where(includepeak(nstarpeaks:npeaks-1) ne 0,nincludepeak_gas)+nstarpeaks ;included gas peaks
    peakradius=sqrt(pixdistxpc(peaks(*,0),peaks(*,1))^2.+pixdistypc(peaks(*,0),peaks(*,1))^2.)
    peakmeanradius=mean(peakradius(inclpeaks))
    lambda_map=2.*sqrt(totalarea/nincludepeak/!pi) ;geometric lambda from map area assuming points are randomly distributed
    lambda_map_star=2.*sqrt(totalarea/nincludepeak_star/!pi) ;geometric lambda from map area assuming points are randomly distributed
    lambda_map_gas=2.*sqrt(totalarea/nincludepeak_gas/!pi) ;geometric lambda from map area assuming points are randomly distributed
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATE APERTURE FLUXES AND AREAS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if calc_ap_flux then begin
    print,' ==> calculating aperture fluxes'
    starflux=dblarr(naperture,npeaks)
    if use_star3 then starflux3=dblarr(naperture,npeaks)
    gasflux=dblarr(naperture,npeaks)
    aperturearea_frac=dblarr(naperture,npeaks)
    for i=0,naperture-1 do begin
        for j=0,npeaks-1 do begin
            kx=peaks(j,0) ;get x index of peak
            ky=peaks(j,1) ;get y index of peak
            starflux(i,j)=smoothstar(i,kx,ky) ;flux of SF tracer at given peak position and aperture size
            if use_star3 then starflux3(i,j)=smoothstar3(i,kx,ky) ;flux of second SF tracer at given peak position and aperture size
            gasflux(i,j)=smoothgas(i,kx,ky) ;flux of gas tracer at given peak position and aperture size
            if calc_ap_area then begin
                aperturearea_frac[i,j]=smoothmask[i,kx,ky] ;area of aperture at given peak position and aperture size
            endif
        endfor
    endfor
endif


;;;;;;;;;;;
;PLOT MAPS;
;;;;;;;;;;;

if generate_plot then begin
    print,' ==> plotting maps'

    set_plot,'ps'
    dim=size(gasmap)
    plot,[gaspeaks(0)],title='!6' ;set default font

    nphi=1001
    phi=dindgen(nphi)/(nphi-1)*2.*!pi
    ring1x=minradius/distance/!dtor/cdelt*(cos(phi)*cos(-posangle)+sin(-posangle)*sin(phi)*cos(inclination))+centre(0)
    ring1y=minradius/distance/!dtor/cdelt*(-cos(phi)*sin(-posangle)+cos(-posangle)*sin(phi)*cos(inclination))+centre(1)
    ring2x=maxradius/distance/!dtor/cdelt*(cos(phi)*cos(-posangle)+sin(-posangle)*sin(phi)*cos(inclination))+centre(0)
    ring2y=maxradius/distance/!dtor/cdelt*(-cos(phi)*sin(-posangle)+cos(-posangle)*sin(phi)*cos(inclination))+centre(1)

    device,filename=figdir+galaxy+'_map_star.eps',xsize=10,ysize=10*dim(2)/dim(1),/color,bits_per_pixel=8,/encapsulated
    rtar=logrange_s+1.
    rmin=min(alog10(smoothstar(peak_res,*,*)),/nan)
    rmax=max(alog10(smoothstar(peak_res,*,*)),/nan)
    plotfits,rundir+'res_'+res(peak_res)+'pc'+path_sep()+starfileshort,rmax-rmin-rtar,0.,1,log=1
    oplot,starpeaks(inclstar,0),starpeaks(inclstar,1),psym=7,color=fsc_color('red'),symsize=.8
    oplot,ring1x,ring1y
    oplot,ring2x,ring2y
    device,/close

    device,filename=figdir+galaxy+'_map_gas.eps',xsize=10,ysize=10*dim(2)/dim(1),/color,bits_per_pixel=8,/encapsulated
    rtar=logrange_g+1.
    rmin=min(alog10(smoothgas(peak_res,*,*)),/nan)
    rmax=max(alog10(smoothgas(peak_res,*,*)),/nan)
    plotfits,rundir+'res_'+res(peak_res)+'pc'+path_sep()+gasfileshort,rmax-rmin-rtar,0.,3,log=1
    oplot,gaspeaks(inclgas-nstarpeaks,0),gaspeaks(inclgas-nstarpeaks,1),psym=7,color=fsc_color('blue'),symsize=.8
    oplot,ring1x,ring1y
    oplot,ring2x,ring2y
    device,/close

    if use_star3 then begin
        device,filename=figdir+galaxy+'_map_star3.eps',xsize=10,ysize=10*dim(2)/dim(1),/color,bits_per_pixel=8
        rtar=1.25
        rmin=min(alog10(smoothstar3(peak_res,*,*)),/nan)
        rmax=max(alog10(smoothstar3(peak_res,*,*)),/nan)
        plotfits,rundir+'res_'+res(peak_res)+'pc'+path_sep()+galaxy+'_star3',rmax-rmin-rtar,0.,1,log=1
        oplot,starpeaks(inclstar,0),starpeaks(inclstar,1),psym=7,color=fsc_color('red'),symsize=.8
        oplot,ring1x,ring1y
        oplot,ring2x,ring2y
        device,/close
    endif
    set_plot,'x'
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATE DISTANCES BETWEEN PEAKS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if get_distances then begin
    distances_all=dblarr(nincludepeak,nincludepeak) ;peak-peak distances (half-filled to avoid double counting)
    for i=0,nincludepeak-1 do begin
        for j=i+1,nincludepeak-1 do begin ;avoid counting pairs twice
            ipeak1=inclpeaks(i)
            ipeak2=inclpeaks(j)
            distances_all(i,j)=norm(reform(peaks(ipeak1,0:1)-peaks(ipeak2,0:1)))*pixtopc
        endfor
        progress,'     ==> calculating distances between peaks',i,nincludepeak-1
    endfor
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DERIVE OBSERVED FLUX RATIO BIASES AND INPUT QUANTITIES FOR FITTING KL14 PRINCIPLE;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if calc_obs then begin
    print,' ==> deriving observed fluxratio biases and input parameters for fitting KL14 principle'
    fluxratio_galaxy=gasfluxtotal/starfluxtotal ;average g/s flux ratio in the included area

    ;start calculating observed fluxratio biases with Monte Carlo sampling of non-overlapping peaks
    totgas_star=dblarr(naperture,nmc) ;total gas flux in apertures centered on SF peaks
    totstar_star=dblarr(naperture,nmc) ;total SF flux in apertures centered on SF peaks
    totgas_gas=dblarr(naperture,nmc) ;total gas flux in apertures centered on gas peaks
    totstar_gas=dblarr(naperture,nmc) ;total SF flux in apertures centered on gas peaks
    apertures_star_mc=dblarr(naperture,nmc)
    apertures_gas_mc=dblarr(naperture,nmc)
    apertures_star=dblarr(naperture) ;mean aperture size centred on SF peaks for ith target aperture size
    apertures_gas=dblarr(naperture) ;mean aperture size centred on gas peaks for ith target aperture size
    meanapgas_star=dblarr(naperture) ;mean gas flux in apertures centered on SF peaks
    meanapstar_star=dblarr(naperture) ;mean SF flux in apertures centered on SF peaks
    meanapgas_gas=dblarr(naperture) ;mean gas flux in apertures centered on gas peaks
    meanapstar_gas=dblarr(naperture) ;mean SF flux in apertures centered on gas peaks
    meantotgas_star=dblarr(naperture) ;mean total gas flux in apertures centered on SF peaks
    meantotstar_star=dblarr(naperture) ;mean total SF flux in apertures centered on SF peaks
    meantotgas_gas=dblarr(naperture) ;mean total gas flux in apertures centered on gas peaks
    meantotstar_gas=dblarr(naperture) ;mean total SF flux in apertures centered on gas peaks
    fluxratio_star=dblarr(naperture)
    fluxratio_gas=dblarr(naperture)
    err_sens_star2=dblarr(naperture)
    err_apgas_star2=dblarr(naperture)
    err_apstar_star2=dblarr(naperture)
    err_apcov_star2=dblarr(naperture)
    err_sens_gas2=dblarr(naperture)
    err_apgas_gas2=dblarr(naperture)
    err_apstar_gas2=dblarr(naperture)
    err_apcov_gas2=dblarr(naperture)
    err_star=dblarr(naperture)
    err_gas=dblarr(naperture)
    err_star_log=dblarr(naperture)
    err_gas_log=dblarr(naperture)
    nexpstar=dblarr(naperture)
    nexpgas=dblarr(naperture)
    nfitstar=dblarr(naperture)
    nfitgas=dblarr(naperture)
    nusestar=dblarr(naperture,nmc) ;number of used stellar peaks in MC sample
    nusegas=dblarr(naperture,nmc) ;number of used gas peaks in MC sample
    nusestariso=dblarr(nmc,nincludepeak_star+1) ;number of used isolated stellar peaks in MC sample
    nusegasiso=dblarr(nmc,nincludepeak_gas+1) ;number of used isolated gas peaks in MC sample
    nstarmc=dblarr(naperture)
    ngasmc=dblarr(naperture)
    fstarover=dindgen(nincludepeak_star+1)/double(nincludepeak_star)
    fgasover=dindgen(nincludepeak_gas+1)/double(nincludepeak_gas)
    betastarmc=dblarr(nmc,nincludepeak_star+1) ;betastar for each MC realisation and each possible value of tover
    betagasmc=dblarr(nmc,nincludepeak_gas+1) ;betagas for each MC realisation and each possible value of tover and tgas
    beta_star=dblarr(nincludepeak_star+1) ;betastar for each possible value of tover, averaged over all MC realisations
    beta_gas=dblarr(nincludepeak_gas+1) ;betagas for each possible value of tover and tgas, averaged over all MC realisations
    for i=peak_res,max_res do begin ;add up all apertures at a given aperture size to obtain flux ratio bias
        usestar=intarr(nmc,nincludepeak_star) ;used stellar peaks in MC sample
        usegas=intarr(nmc,nincludepeak_gas) ;used gas peaks in MC sample
        for j=0,nmc-1 do begin ;start MC drawing of non-overlapping peak samples (separately for stars and gas)
            rnd_star=randomu(seedarr(iseed),nincludepeak_star) ;random numbers for drawing stellar peaks in a random order
            iseed+=1 ;seed was used, on to the next one
            possiblestar=dindgen(nincludepeak_star) ;stellar peak indices
            candidate=floor(rnd_star(0)*(nincludepeak_star)) ;first stellar peak
            usestar(j,0)=possiblestar(candidate) ;include candidate peak
            remove,candidate,possiblestar ;remove used peak from list of possible candidates
            nusestar(i,j)+=1 ;got one stellar peak
            for k=1,nincludepeak_star-1 do begin ;start looping over other stellar peaks
                candidate=floor(rnd_star(k)*(nincludepeak_star-k)) ;randomly-drawn, new candidate stellar peak
                usedstar=reform(usestar(j,0:nusestar(i,j)-1)) ;other stellar peaks used so far
                distances_candidate_usestar=[reform(distances_all(possiblestar(candidate),usedstar)),reform(distances_all(usedstar,possiblestar(candidate)))] ;distances between candidate and used peaks
                distances_candidate_usestar=distances_candidate_usestar(where(distances_candidate_usestar gt 0.)) ;remove zeroes
                mindist=min(distances_candidate_usestar) ;smallest distance between candidate and any of the used peaks
                if mindist gt apertures(i) then begin ;if at least one aperture size away from all used peaks, then:
                    usestar(j,nusestar(i,j))=possiblestar(candidate) ;include candidate peak
                    nusestar(i,j)+=1 ;got one more stellar peak
                endif
                if n_elements(possiblestar) gt 1 then remove,candidate,possiblestar ;remove used peak from list of possible candidates
            endfor
            rnd_gas=randomu(seedarr(iseed),nincludepeak_gas) ;random numbers for drawing gas peaks in a random order
            iseed+=1 ;seed was used, on to the next one
            possiblegas=dindgen(nincludepeak_gas) ;gas peak indices
            candidate=floor(rnd_gas(0)*(nincludepeak_gas)) ;first gas peak
            usegas(j,0)=possiblegas(candidate) ;include candidate peak
            remove,candidate,possiblegas ;remove used peak from list of possible candidates
            nusegas(i,j)+=1 ;got one gas peak
            for k=1,nincludepeak_gas-1 do begin ;start looping over other gas peaks
                candidate=floor(rnd_gas(k)*(nincludepeak_gas-k)) ;randomly-drawn, new candidate gas peak
                usedgas=reform(usegas(j,0:nusegas(i,j)-1)) ;other gas peaks used so far
                distances_candidate_usegas=[reform(distances_all(possiblegas(candidate)+nincludepeak_star,usedgas+nincludepeak_star)), $
                                            reform(distances_all(usedgas+nincludepeak_star,possiblegas(candidate)+nincludepeak_star))]
                                            ;distances between candidate and used peaks
                distances_candidate_usegas=distances_candidate_usegas(where(distances_candidate_usegas gt 0.)) ;remove zeroes
                mindist=min(distances_candidate_usegas) ;smallest distance between candidate and any of the used peaks
                if mindist gt apertures(i) then begin ;if at least one aperture size away from all used peaks, then:
                    usegas(j,nusegas(i,j))=possiblegas(candidate) ;include candidate peak
                    nusegas(i,j)+=1 ;got one more gas peak
                endif
                if n_elements(possiblegas) gt 1 then remove,candidate,possiblegas ;remove used peak from list of possible candidates
            endfor
            usedstar=reform(usestar(j,0:nusestar(i,j)-1)) ;IDs of SF peaks used in this MC sample
            nusedstar=n_elements(usedstar) ;number of SF peaks used in this MC sample
            usedstarcount=dindgen(nusedstar) ;array to count used SF peaks
            usedgas=reform(usegas(j,0:nusegas(i,j)-1)) ;IDs of gas peaks used in this MC sample
            nusedgas=n_elements(usedgas) ;number of gas peaks used in this MC sample
            usedgascount=dindgen(nusedgas) ;array to count used gas peaks
            totgas_star(i,j)=max([tiny,total(gasflux(i,inclstar(usedstar)),/nan)]) ;total gas flux in apertures centered on SF peaks (must be > 0)
            totstar_star(i,j)=max([tiny,total(starflux(i,inclstar(usedstar)),/nan)]) ;total SF flux in apertures centered on SF peaks (must be > 0)
            totgas_gas(i,j)=max([tiny,total(gasflux(i,inclgas(usedgas)),/nan)]) ;total gas flux in apertures centered on gas peaks (must be > 0)
            totstar_gas(i,j)=max([tiny,total(starflux(i,inclgas(usedgas)),/nan)]) ;total SF flux in apertures centered on gas peaks (must be > 0)
            apertures_star_mc(i,j)=apertures[i]*sqrt(mean(aperturearea_frac[i,inclstar(usedstar)])) ;mean aperture size centred on SF peaks for ith target aperture size (order of sqrt(mean) is intentional)
            apertures_gas_mc(i,j)=apertures[i]*sqrt(mean(aperturearea_frac[i,inclgas(usedgas)])) ;mean aperture size centred on SF peaks for ith target aperture size (order of sqrt(mean) is intentional)
            if i eq peak_res then begin ;determine Monte Carlo betastar and betagas based on fractional timescale coverage of each phase projected onto list of peaks sorted by gas fraction-to-star fraction ratio
                starratio=reform(gasflux(i,inclstar(usedstar))/starflux(i,inclstar(usedstar))) ;ratio of gas flux to SF flux around SF peak
                sortstarratio=reverse(sort(starratio)) ;reverse-sort the above, i.e. in order of decreasing ratio and thus of increasing SF peak prominence relative to associated gas peak
                nstarover=fstarover*nusedstar ;number of SF peaks in overlap
                kminstar=min(where(nstarover gt 1)) ;avoid case of no overlap
                kmaxstar=max(where(nstarover lt nusedstar-1)) ;avoid case of all overlap
                for k=kminstar,kmaxstar do begin
                    star_critcon=interpol(starratio(sortstarratio),usedstarcount,nstarover(k)) ;critical difference in contrast with respect to background between gas and SF emission to call it an isolated SF peak
                    i_stariso=where(starratio le star_critcon,nusestarisotemp) ;IDs of isolated SF peaks
                    nusestariso(j,k)=nusestarisotemp ;number of isolated SF peaks
                    if nusestariso(j,k) gt 0 && nusestariso(j,k) lt nusedstar then begin ;if isolated SF peaks exist and not all SF peaks are isolated (should always hold for loop interval)
                        i_starover=usedstarcount ;create array with overlapping SF peak IDs
                        remove,reform(i_stariso),i_starover ;remove isolated SF peaks from overlapping SF peak ID array
                        betastarmc(j,k)=mean(starflux(i,inclstar(usedstar(i_starover))))/mean(starflux(i,inclstar(usedstar(i_stariso)))) ;define betastar as ratio between mean overlapping and isolated SF peak fluxes
                    endif else betastarmc(j,k)=1. ;if no isolated SF peaks exist or all SF peaks are isolated, set betastar=1
                endfor
                if kminstar gt 0 then betastarmc(j,0:kminstar-1)=betastarmc(j,kminstar) ;extend beta for no overlap case(s)
                if kmaxstar lt nusedstar then betastarmc(j,kmaxstar+1:nusedstar)=betastarmc(j,kmaxstar) ;extend beta for all overlap case(s)
                gasratio=reform(starflux(i,inclgas(usedgas))/gasflux(i,inclgas(usedgas))) ;ratio of SF flux to excess gas flux around gas peak
                sortgasratio=reverse(sort(gasratio)) ;reverse-sort the above, i.e. in order of decreasing ratio and thus of increasing gas peak prominence relative to associated SF peak
                ngasover=fgasover*nusedgas ;number of gas peaks in overlap
                kmingas=min(where(ngasover gt 1)) ;avoid case of no overlap
                kmaxgas=max(where(ngasover lt nusedgas-1)) ;avoid case of all overlap
                for k=kmingas,kmaxgas do begin
                    gas_critcon=interpol(gasratio(sortgasratio),usedgascount,ngasover(k)) ;critical difference in contrast with respect to background between gas and SF emission to call it an isolated gas peak
                    i_gasiso=where(gasratio le gas_critcon,nusegasisotemp) ;IDs of isolated gas peaks
                    nusegasiso(j,k)=nusegasisotemp ;number of isolated gas peaks
                    if nusegasiso(j,k) gt 0 && nusegasiso(j,k) lt nusedgas then begin ;if isolated gas peaks exist and not all gas peaks are isolated (should always hold for loop interval)
                        i_gasover=usedgascount ;create array with overlapping gas peak IDs
                        remove,reform(i_gasiso),i_gasover ;remove isolated gas peaks from overlapping gas peak ID array
                        betagasmc(j,k)=mean(gasflux(i,inclgas(usedgas(i_gasover))))/mean(gasflux(i,inclgas(usedgas(i_gasiso)))) ;define betagas as ratio between mean overlapping and isolated gas peak fluxes
                    endif else betagasmc(j,k)=1. ;if no isolated gas peaks exist or all gas peaks are isolated, set betagas=1
                endfor
                if kmingas gt 0 then betagasmc(j,0:kmingas-1)=betagasmc(j,kmingas) ;extend beta for no overlap case(s)
                if kmaxgas lt nusedgas then betagasmc(j,kmaxgas+1:nusedgas)=betagasmc(j,kmaxgas) ;extend beta for all overlap case(s)
            endif
            progress,'     ==> Monte-Carlo sampling peak maps to get uncorrelated peak samples',j+i*nmc,nmc*naperture-1
        endfor

        ;APERTURE AREAS
        if calc_ap_area then begin ;calculate aperture area
            apertures_star[i]=mean(apertures_star_mc[i,*]) ;mean aperture size centred on SF peaks for ith target aperture size (order of sqrt(mean) is intentional)
            apertures_gas[i]=mean(apertures_gas_mc[i,*]) ;mean aperture size centred on gas peaks for ith target aperture size (order of sqrt(mean) is intentional)
        endif else begin ;use the area of user-defined aperture areas
            apertures_star[i] = apertures[i] ;replicate non-mask functionality for fitKL14
            apertures_gas[i] = apertures[i]  ;replicate non-mask functionality for fitKL14
        endelse

        ;APERTURES CENTERED ON SF PEAKS -- obtain data points and errors
        nstarmc(i)=mean(nusestar(i,*)) ;mean number of SF peaks
        for k=0,nincludepeak_star do begin
            beta_star_use=where(betastarmc(*,k) gt 0.,n_beta_star_use) ;avoid including MC realisations where no beta_star could be estimated
            if n_beta_star_use gt tiny then beta_star(k)=mean(betastarmc(beta_star_use,k)) ;betastar for each possible value of tover, averaged over all MC realisations
        endfor
        meanapgas_star(i)=mean(gasflux(i,inclstar)) ;mean gas flux density per aperture in apertures centered on SF peaks
        meanapstar_star(i)=mean(starflux(i,inclstar)) ;mean SF flux density per aperture in apertures centered on SF peaks
        meantotgas_star(i)=mean(totgas_star(i,*)) ;mean total gas flux density across all apertures centered on SF peaks
        meantotstar_star(i)=mean(totstar_star(i,*)) ;mean total SF flux density across all apertures centered on SF peaks
        fluxratio_star(i)=meantotgas_star(i)/meantotstar_star(i) ;gas-to-stellar flux ratio across all apertures centered on SF peaks, averaged over all MC realisations
        ;obtain individual error components
        err_sensgas_star=sensgas*(apertures_star(peak_res)/apertures_star(i))*sqrt(nstarmc(i)) ;total gas flux error due to sensitivity
        err_sensstar_star=sensstar*(apertures_star(peak_res)/apertures_star(i))*sqrt(nstarmc(i)) ;total SF flux error due to sensitivity
        err_apgas_star=stdev(gasflux(i,inclstar)) ;standard deviation of the gas flux in individual apertures centered on SF peaks
        err_apstar_star=stdev(starflux(i,inclstar)) ;standard deviation of the stellar flux in individual apertures centered on SF peaks
        err_totgas_star=stdev(totgas_star(i,*)) ;standard deviation of the total gas flux for apertures centered on SF peaks, across all MC realisations -- ancillary quantity
        err_totstar_star=stdev(totstar_star(i,*)) ;standard deviation of the total stellar flux for apertures centered on SF peaks, across all MC realisations -- ancillary quantity
        ;obtain relative error terms and covariance of numerator and denominator
        nexpstar(i)=nstarmc(i) ;the number of experiments that we are averaging over
        err_sens_star2(i)=((err_sensgas_star/meantotgas_star(i))^2.+(err_sensstar_star/meantotstar_star(i))^2.) ;relative error due to sensitivity
        err_apgas_star2(i)=(err_apgas_star/meantotgas_star(i))^2.*nexpstar(i) ;relative error due to standard deviation of the aperture population (gas)
        err_apstar_star2(i)=(err_apstar_star/meantotstar_star(i))^2.*nexpstar(i) ;relative error due to standard deviation of the aperture population (stars)
        err_apcov_star2(i)=-2.*correlate(gasflux(i,inclstar),starflux(i,inclstar),/covariance)/meantotgas_star(i)/meantotstar_star(i)*nexpstar(i) ;covariance
        ;get total error
        err_star(i)=sqrt(err_sens_star2(i)+err_apgas_star2(i)+err_apstar_star2(i)+err_apcov_star2(i))*fluxratio_star(i)
        err_star_log(i)=err_star(i)/fluxratio_star(i)/alog(10.) ;error in log space

        ;APERTURES CENTERED ON GAS PEAKS -- obtain data points and errors
        ngasmc(i)=mean(nusegas(i,*)) ;mean number of gas peaks
        for k=0,nincludepeak_gas do begin
            beta_gas_use=where(betagasmc(*,k) gt 0.,n_beta_gas_use) ;avoid including MC realisations where no beta_gas could be estimated
            if n_beta_gas_use gt tiny then beta_gas(k)=mean(betagasmc(beta_gas_use,k)) ;betagas for each possible value of tover and tgas, averaged over all MC realisations
        endfor
        meanapgas_gas(i)=mean(gasflux(i,inclgas)) ;mean gas flux per aperture in apertures centered on gas peaks
        meanapstar_gas(i)=mean(starflux(i,inclgas)) ;mean SF flux per aperture in apertures centered on gas peaks
        meantotgas_gas(i)=mean(totgas_gas(i,*)) ;mean total gas flux across all apertures centered on gas peaks
        meantotstar_gas(i)=mean(totstar_gas(i,*)) ;mean total SF flux across all apertures centered on gas peaks
        fluxratio_gas(i)=meantotgas_gas(i)/meantotstar_gas(i) ;gas-to-stellar flux ratio across all apertures centered on gas peaks, averaged over all MC realisations
        ;obtain individual error components
        err_sensgas_gas=sensgas*(apertures_gas(peak_res)/apertures_gas(i))*sqrt(ngasmc(i)) ;total SF flux error due to sensitivity
        err_sensstar_gas=sensstar*(apertures_gas(peak_res)/apertures_gas(i))*sqrt(ngasmc(i)) ;total gas flux error due to sensitivity
        err_apgas_gas=stdev(gasflux(i,inclgas)) ;standard deviation of the gas flux in individual apertures centered on gas peaks
        err_apstar_gas=stdev(starflux(i,inclgas)) ;standard deviation of the stellar flux in individual apertures centered on gas peaks
        err_totgas_gas=stdev(totgas_gas(i,*)) ;standard deviation of the total gas flux for apertures centered on gas peaks, across all MC realisations -- ancillary quantity
        err_totstar_gas=stdev(totstar_gas(i,*)) ;standard deviation of the total stellar flux for apertures centered on gas peaks, across all MC realisations -- ancillary quantity
        ;obtain relative error terms and covariance of numerator and denominator
        nexpgas(i)=ngasmc(i) ;the number of experiments that we are averaging over
        err_sens_gas2(i)=((err_sensgas_gas/meantotgas_gas(i))^2.+(err_sensstar_gas/meantotstar_gas(i))^2.) ;relative error due to sensitivity
        err_apgas_gas2(i)=(err_apgas_gas/meantotgas_gas(i))^2.*nexpgas(i) ;relative error due to standard deviation of the aperture population (gas)
        err_apstar_gas2(i)=(err_apstar_gas/meantotstar_gas(i))^2.*nexpgas(i) ;relative error due to standard deviation of the aperture population (stars)
        err_apcov_gas2(i)=-2.*correlate(gasflux(i,inclgas),starflux(i,inclgas),/covariance)/meantotgas_gas(i)/meantotstar_gas(i)*nexpgas(i) ;covariance
        ;get total error
        err_gas(i)=sqrt(err_sens_gas2(i)+err_apgas_gas2(i)+err_apstar_gas2(i)+err_apcov_gas2(i))*fluxratio_gas(i)
        err_gas_log(i)=err_gas(i)/fluxratio_gas(i)/alog(10.) ;error in log space
    endfor
    corrgas_gas=dblarr(naperture,naperture)
    corrstar_gas=dblarr(naperture,naperture)
    corrgas_star=dblarr(naperture,naperture)
    corrstar_star=dblarr(naperture,naperture)
    for i=peak_res,max_res do begin
        for j=i,max_res do begin
            corrgas_gas(i,j)=meantotgas_gas(i)/meantotgas_gas(j)*ngasmc(j)/ngasmc(i)*apertures_gas(i)^2./apertures_gas(j)^2.
            corrstar_gas(i,j)=meantotstar_gas(i)/meantotstar_gas(j)*ngasmc(j)/ngasmc(i)*apertures_gas(i)^2./apertures_gas(j)^2.
            corrgas_star(i,j)=meantotgas_star(i)/meantotgas_star(j)*nstarmc(j)/nstarmc(i)*apertures_star(i)^2./apertures_star(j)^2.
            corrstar_star(i,j)=meantotstar_star(i)/meantotstar_star(j)*nstarmc(j)/nstarmc(i)*apertures_star(i)^2./apertures_star(j)^2.
            corrgas_gas(j,i)=corrgas_gas(i,j)
            corrstar_gas(j,i)=corrstar_gas(i,j)
            corrgas_star(j,i)=corrgas_star(i,j)
            corrstar_star(j,i)=corrstar_star(i,j)
        endfor
    endfor
    for i=peak_res,max_res do begin ;total number of independent data points for gas and SF peaks
        nfitgas(i)=.5/total(corrgas_gas(i,*))+.5/total(corrstar_gas(i,*)) ;take the mean of the number of independent datapoints for the numerator and the denominator (for N=>inf 1/2+1/2N is right, but 2/(1+N) is wrong)
        nfitstar(i)=.5/total(corrgas_star(i,*))+.5/total(corrstar_star(i,*))
    endfor
    surfglobals=meantotstar_star/nstarmc/starfluxtotal*totalarea ;surface density of SF peak relative to entire map
    surfglobalg=meantotgas_gas/ngasmc/gasfluxtotal*totalarea ;surface density of SF peak relative to entire map
    surfcontrasts=dblarr(naperture)
    surfcontrastg=dblarr(naperture)
    for i=peak_res,max_res do begin
        surfcontrasts(i)=max([1.,meantotstar_star(peak_res)/meantotstar_star(i)*nstarmc(i)/nstarmc(peak_res)],/nan) ;surface density contrast of SF peak
        surfcontrastg(i)=max([1.,meantotgas_gas(peak_res)/meantotgas_gas(i)*ngasmc(i)/ngasmc(peak_res)],/nan) ;surface density contrast of gas peak
    endfor
    bias_star=fluxratio_star/fluxratio_galaxy
    bias_gas=fluxratio_gas/fluxratio_galaxy
    err_star=err_star/fluxratio_galaxy ;scale linear error -- log error is unchanged
    err_gas=err_gas/fluxratio_galaxy ;scale linear error -- log error is unchanged
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FIT KL14 PRINCIPLE TO DATA;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if calc_fit then begin
    print,' ==> fitting KL14 principle'

    fit=fitKL14(fluxratio_star[fitap]/fluxratio_galaxy,fluxratio_gas[fitap]/fluxratio_galaxy, $
                err_star_log[fitap],err_gas_log[fitap],tstariso,beta_star,beta_gas,fstarover,fgasover,apertures_star[fitap],apertures_gas[fitap], $
                surfcontrasts[fitap],surfcontrastg[fitap],peak_prof,tstar_incl,tgasmini,tgasmaxi,tovermini,tovermaxi,nfitstar[fitap],nfitgas[fitap], $
                ndepth,ntry,galaxy,figdir,generate_plot,outputdir,arrdir,window_plot)
    tgas=fit(1)
    tgas_errmin=sqrt(fit(2)^2.+(tgas*tstariso_rerrmin)^2.)
    tgas_errmax=sqrt(fit(3)^2.+(tgas*tstariso_rerrmax)^2.)
    tover=fit(4)
    tover_errmin=sqrt(fit(5)^2.+(tover*tstariso_rerrmin)^2.)
    tover_errmax=sqrt(fit(6)^2.+(tover*tstariso_rerrmax)^2.)
    lambda=fit(7)
    lambda_errmin=fit(8)
    lambda_errmax=fit(9)
    if tstar_incl eq 0 then tstar=tstariso+tover else tstar=tstariso
    ttotal=tgas+tstar-tover
    beta_star_fit=fit(10)
    beta_gas_fit=fit(11)
    surfcontrasts_fit=fit(12)
    surfcontrastg_fit=fit(13)
    redchi2=fit(0)
    fit=fit(0:9)

    zetas=f_zeta(surfcontrasts_fit,ttotal,tstar,peak_prof)
    zetag=f_zeta(surfcontrastg_fit,ttotal,tgas,peak_prof)
    rpeaks=f_rpeak(zetas,lambda)
    rpeakg=f_rpeak(zetag,lambda)

    if tstar_incl eq 0 then terrs=0. else terrs=0.
    if tstar_incl eq 0 then terrg=0. else terrg=0.

    if lambda/2. le rpeaks || lambda/2. le rpeakg then begin
        print,''
        print,' WARNING: derived lambda is smaller than peak dispersion, suggests inadequate map resolution'
    endif
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculate fgmc and fcl;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if diffuse_frac eq 1 then begin
  if f_filter_type eq 0 then filter_choice = 'butterworth' else $
      if f_filter_type eq 1 then filter_choice = 'gaussian' else $
          if f_filter_type eq 2 then filter_choice = 'ideal'

  if (diffuse_quant eq 0) then begin ; diffuse_quant =0 -> flux
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;1) method using derivephys and arrays ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      flux_fraction_calc, filter_choice, bw_order, 'high' $ ; description of the filter   $
       , gasfiletot, starfiletot $ ; input image (use masked image, but not regridded/smoothened)
       , distance, inclination, astr_tolerance $
       , arrdir $
       , filter_len_conv $  ; fourier length conversion factor
       , lambda $
       , gas_flux_frac, star_flux_frac $
       , /save_arrays

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ; 2) direct method using lambda errors ;
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     fourier_diffuse_fraction, filter_choice, bw_order, 'high' $ ; description of the filter   $
       , gasfiletot $ ; input image
       , distance, inclination, astr_tolerance $ ;
       , filter_len_conv $  ; fourier length conversion factor
       , lambda $ ; lambda
       , gas_flux_frac $ ; output flux and diffuse fraction
       , lambda_errmin = lambda_errmin, lambda_errmax = lambda_errmax $ ; lambda errors
       , flux_frac_errmax = gas_flux_frac_errmax, flux_frac_errmin = gas_flux_frac_errmin ; output flux fraction errors

     fourier_diffuse_fraction, filter_choice, bw_order, 'high' $ ; description of the filter   $
       , starfiletot $ ; input image
       , distance, inclination, astr_tolerance $ ;
       , filter_len_conv $  ; fourier length conversion factor
       , lambda $ ; lambda
       , star_flux_frac $ ; output flux and diffuse fraction
       , lambda_errmin = lambda_errmin, lambda_errmax = lambda_errmax $ ; lambda errors
       , flux_frac_errmax = star_flux_frac_errmax, flux_frac_errmin = star_flux_frac_errmin ; output flux fraction errors
     fcl=star_flux_frac
     fcl_errmin=star_flux_frac_errmin
     fcl_errmax=star_flux_frac_errmax

     fgmc=gas_flux_frac
     fgmc_errmin=gas_flux_frac_errmin
     fgmc_errmax=gas_flux_frac_errmax


   endif else if (diffuse_quant eq 1) then begin ; diffuse_quant =1 -> power
      power_fraction_calc, gasfiletot, starfiletot $ ; input image (use masked image, but not regridded/smoothened)
          , distance, inclination, astr_tolerance $
          , arrdir $ ; directory to restore lambda array from
          , filter_len_conv $  ; fourier length conversion factor
          , lambda $
          , gas_flux_frac, star_flux_frac $
          , /save_arrays

          fcl=star_flux_frac
          fgmc=gas_flux_frac


  ;     power_fraction_calc, masked_path_gas $ ; input image (use masked image, but not regridded/smoothened)
  ;         , distance, inclination, astr_tolerance $
  ;         , arrdir $ ; directory to restore lambda array from
  ;         , fourier_length_conv $  ; fourier length conversion factor
  ;         , lambda $
  ;         , fit_power_frac $
  ;         , power_frac_arr $
  ;         , save_array = save_array
   ;
  ;  print, fit_power_frac
  ;  print, power_frac_arr
  ;  stop

   endif else begin
      f_error, 'The current value of diffuse_quant (', diffuse_quant, ') is not valued select 0 (flux) or 1 (power)'
   endelse
endif else begin ; calcuate fgmc and fcl not using fourier filtering
    if map_units eq 0 then begin
        fcl=1. ;SF flux tracer emission from peaks -- will become functional after including Fourier filtering
        fcl_errmin=tiny ; no error calculated as 1.0 is assumed so use tiny as dummy error
        fcl_errmax=tiny ; no error calculated as 1.0 is assumed

        fgmc=1. ;gas flux tracer emission from peaks -- will become functional after including Fourier filtering
        fgmc_errmin=tiny ; no error calculated as 1.0 is assumed so use tiny as dummy error
        fgmc_errmax=tiny ; no error calculated as 1.0 is assumed so use tiny as dummy error
    endif
    if map_units eq 1 then begin
        fcl=total(starmap,/nan)*convstar/sfr_galaxy ;SF flux tracer emission from peaks -- will become functional after including Fourier filtering
        fcl_errmin=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations
        fcl_errmax=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations

        fgmc=total(gasmap,/nan)*convgas/mgas_galaxy ;gas flux tracer emission from peaks -- will become functional after including Fourier filtering
        fgmc_errmin=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
        fgmc_errmax=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
    endif
    if map_units eq 2 then begin
        fcl=total(starmap,/nan)*convstar/mgas_galaxy1 ;SF flux tracer emission from peaks -- will become functional after including Fourier filtering
        fcl_errmin=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations
        fcl_errmax=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations

        fgmc=total(gasmap,/nan)*convgas/mgas_galaxy2 ;gas flux tracer emission from peaks -- will become functional after including Fourier filtering
        fgmc_errmin=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
        fgmc_errmax=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
    endif
    if map_units eq 3 then begin
        fcl=total(starmap,/nan)*convstar/sfr_galaxy1 ;SF flux tracer emission from peaks -- will become functional after including Fourier filtering
        fcl_errmin=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations
        fcl_errmax=fcl*convstar_err/convstar        ; pre-fourier-filtering error calculations

        fgmc=total(gasmap,/nan)*convgas/sfr_galaxy2 ;gas flux tracer emission from peaks -- will become functional after including Fourier filtering
        fgmc_errmin=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
        fgmc_errmax=fgmc*convgas_err/convgas        ; pre-fourier-filtering error calculations
      endif

      fsingle_array_make, arrdir, fcl, fgmc ; save out arrays of fcl and fgmc with same number of elements of lambdaarr


endelse








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATE DERIVED PHYSICAL QUANTITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if derive_phys then begin
    print,' ==> calculate derived physical quantities'
    ; now set surfsfr and surfgas
    if map_units eq 0 then begin
        surfsfr=tiny
        surfsfr_err=tiny
        surfgas=tiny
        surfgas_err=tiny
    endif
    if map_units eq 1 then begin
        surfsfr=sfr_galaxy/totalarea
        surfsfr_err=sfr_galaxy_err/totalarea
        surfgas=mgas_galaxy/totalarea
        surfgas_err=mgas_galaxy_err/totalarea
    endif
    if map_units eq 2 then begin
        surfsfr=mgas_galaxy1/totalarea
        surfsfr_err=mgas_galaxy1_err/totalarea
        surfgas=mgas_galaxy2/totalarea
        surfgas_err=mgas_galaxy2_err/totalarea
    endif
    if map_units eq 3 then begin
        surfsfr=sfr_galaxy1/totalarea
        surfsfr_err=sfr_galaxy1_err/totalarea
        surfgas=sfr_galaxy2/totalarea
        surfgas_err=sfr_galaxy2_err/totalarea
    endif
    if map_units gt 0 then begin
        ext=[surfsfr*totalarea,surfsfr_err*totalarea,surfsfr_err*totalarea,surfgas*totalarea,surfgas_err*totalarea,surfgas_err*totalarea,surfsfr,surfsfr_err,surfsfr_err,surfgas,surfgas_err,surfgas_err]+tiny

        der=derivephys(surfsfr,surfsfr_err,surfgas,surfgas_err,totalarea,tgas,tover,lambda,beta_star,beta_gas,fstarover,fgasover,fcl,fgmc, $
                       tstariso,tstariso_rerrmin,tstariso_rerrmax,tstar_incl, $
                       surfglobals[fitap],surfglobalg[fitap],surfcontrasts[fitap],surfcontrastg[fitap],apertures_star[fitap],apertures_gas[fitap], $
                       lighttomass,momratetomass,peak_prof,ntry,nphysmc,galaxy,outputdir,arrdir,figdir,map_units)
    endif else der=derivephys(surfsfr,surfsfr_err,surfgas,surfgas_err,totalarea,tgas,tover,lambda,beta_star,beta_gas,fstarover,fgasover,fcl,fgmc, $
                   tstariso,tstariso_rerrmin,tstariso_rerrmax,tstar_incl, $
                   surfglobals[fitap],surfglobalg[fitap],surfcontrasts[fitap],surfcontrastg[fitap],apertures_star[fitap],apertures_gas[fitap], $
                   lighttomass,momratetomass,peak_prof,ntry,nphysmc,galaxy,outputdir,arrdir,figdir,map_units)
    aux=[nincludepeak_star,nincludepeak_gas,lap_min]
endif


;;;;;;;;;;;;;;;;;;;;
;WRITE OUTPUT FILES;
;;;;;;;;;;;;;;;;;;;;

if write_output then begin
    print,' ==> writing output files'
    ;write table output row and output file
    fitqty=['redchi2', $
            'tgas','tgas_errmin','tgas_errmax', $
            'tover','tover_errmin','tover_errmax', $
            'lambda','lambda_errmin','lambda_errmax']
    fitunit=['','Myr','Myr','pc']
    fitad='Best-fitting'
    fitstrings=[fitqty(0), $
                fitad+' '+fitqty(1)+', '+fitqty(2)+', '+fitqty(3)+' ['+fitunit(1)+']', $
                fitad+' '+fitqty(4)+', '+fitqty(5)+', '+fitqty(6)+' ['+fitunit(2)+']', $
                fitad+' '+fitqty(7)+', '+fitqty(8)+', '+fitqty(9)+' ['+fitunit(3)+']']

    extunit=['Msun yr^-1','Msun','Msun yr^-1 pc^-2','Msun pc^-2']
    extad='Derived'
    if map_units eq 1 then begin
        extqty=['sfr_galaxy','sfr_galaxy_errmin','sfr_galaxy_errmax', $
                'mgas_galaxy','mgas_galaxy_errmin','mgas_galaxy_errmax', $
                'surfsfr','surfsfr_errmin','surfsfr_errmax', $
                'surfgas','surfgas_errmin','surfgas_errmax']
    endif
    if map_units eq 2 then begin
        extunit(0)='Msun'
        extunit(2)='Msun pc^-2'
        extqty=['mgas_galaxy1','mgas_galaxy1_errmin','mgas_galaxy1_errmax', $
                'mgas_galaxy2','mgas_galaxy2_errmin','mgas_galaxy2_errmax', $
                'surfgas1','surfgas1_errmin','surfgas1_errmax', $
                'surfgas2','surfgas2_errmin','surfgas2_errmax']
    endif
    if map_units eq 3 then begin
        extunit(1)='Msun yr^-1'
        extunit(3)='Msun yr^-1 pc^-2'
        extqty=['sfr_galaxy1','sfr_galaxy1_errmin','sfr_galaxy1_errmax', $
                'sfr_galaxy2','sfr_galaxy2_errmin','sfr_galaxy2_errmax', $
                'surfsfr1','surfsfr1_errmin','surfsfr1_errmax', $
                'surfsfr2','surfsfr2_errmin','surfsfr2_errmax']
    endif
    if map_units gt 0 then extstrings=[extad+' '+extqty(0)+', '+extqty(1)+', '+extqty(2)+' ['+extunit(0)+']', $
                            extad+' '+extqty(3)+', '+extqty(4)+', '+extqty(5)+' ['+extunit(1)+']', $
                            extad+' '+extqty(6)+', '+extqty(7)+', '+extqty(8)+' ['+extunit(2)+']', $
                            extad+' '+extqty(9)+', '+extqty(10)+', '+extqty(11)+' ['+extunit(3)+']']

    if map_units gt 0 then begin
        derqty=['tstar','tstar_errmin','tstar_errmax', $
                'ttotal','ttotal_errmin','ttotal_errmax', $
                'betastar','betastar_errmin','betastar_errmax', $
                'betagas','betagas_errmin','betagas_errmax', $
                'surfglobalstar','surfglobalstar_errmin','surfglobalstar_errmax', $
                'surfglobalgas','surfglobalgas_errmin','surfglobalgas_errmax', $
                'surfconstar','surfconstar_errmin','surfconstar_errmax', $
                'surfcongas','surfcongas_errmin','surfcongas_errmax', $

                'fcl','fcl_errmin','fcl_errmax', $
                'fgmc','fgmc_errmin','fgmc_errmax', $
                ; 'fdiffstar','fdiffstar_errmin','fdiffstar_errmax', $
                ; 'fdiffgas','fdiffgas_errmin','fdiffgas_errmax', $
                ;add: fdiffstar
                ;add: fdiffgas

                'zetastar','zetastar_errmin','zetastar_errmax', $
                'zetagas','zetagas_errmin','zetagas_errmax', $
                'rpeakstar','rpeakstar_errmin','rpeakstar_errmax', $
                'rpeakgas','rpeakgas_errmin','rpeakgas_errmax', $
                'vfb','vfb_errmin','vfb_errmax', $
                'tdepl','tdepl_errmin','tdepl_errmax', $
                'esf','esf_errmin','esf_errmax', $
                'mdotsf','mdotsf_errmin','mdotsf_errmax', $
                'mdotfb','mdotfb_errmin','mdotfb_errmax', $
                'etainst','etainst_errmin','etainst_errmax', $
                'etaavg','etaavg_errmin','etaavg_errmax', $
                'chie','chie_errmin','chie_errmax', $
                'chip','chip_errmin','chip_errmax']
        derunit=['Myr','Myr','','','','','','','','','','','pc','pc','km s^-1','Gyr','','Msun yr^-1','Msun yr^-1','','','','']
        derad='Derived'
        derstrings=[derad+' '+derqty(0)+', '+derqty(1)+', '+derqty(2)+' ['+derunit(0)+']', $
                    derad+' '+derqty(3)+', '+derqty(4)+', '+derqty(5)+' ['+derunit(1)+']', $
                    derad+' '+derqty(6)+', '+derqty(7)+', '+derqty(8)+' ['+derunit(2)+']', $
                    derad+' '+derqty(9)+', '+derqty(10)+', '+derqty(11)+' ['+derunit(3)+']', $
                    derad+' '+derqty(12)+', '+derqty(13)+', '+derqty(14)+' ['+derunit(4)+']', $
                    derad+' '+derqty(15)+', '+derqty(16)+', '+derqty(17)+' ['+derunit(5)+']', $
                    derad+' '+derqty(18)+', '+derqty(19)+', '+derqty(20)+' ['+derunit(6)+']', $
                    derad+' '+derqty(21)+', '+derqty(22)+', '+derqty(23)+' ['+derunit(7)+']', $
                    derad+' '+derqty(24)+', '+derqty(25)+', '+derqty(26)+' ['+derunit(8)+']', $
                    derad+' '+derqty(27)+', '+derqty(28)+', '+derqty(29)+' ['+derunit(9)+']', $
                    derad+' '+derqty(30)+', '+derqty(31)+', '+derqty(32)+' ['+derunit(10)+']', $
                    derad+' '+derqty(33)+', '+derqty(34)+', '+derqty(35)+' ['+derunit(11)+']', $
                    derad+' '+derqty(36)+', '+derqty(37)+', '+derqty(38)+' ['+derunit(12)+']', $
                    derad+' '+derqty(39)+', '+derqty(40)+', '+derqty(41)+' ['+derunit(13)+']', $
                    derad+' '+derqty(42)+', '+derqty(43)+', '+derqty(44)+' ['+derunit(14)+']', $
                    derad+' '+derqty(45)+', '+derqty(46)+', '+derqty(47)+' ['+derunit(15)+']', $
                    derad+' '+derqty(48)+', '+derqty(49)+', '+derqty(50)+' ['+derunit(16)+']', $
                    derad+' '+derqty(51)+', '+derqty(52)+', '+derqty(53)+' ['+derunit(17)+']', $
                    derad+' '+derqty(54)+', '+derqty(55)+', '+derqty(56)+' ['+derunit(18)+']', $
                    derad+' '+derqty(57)+', '+derqty(58)+', '+derqty(59)+' ['+derunit(19)+']', $
                    derad+' '+derqty(60)+', '+derqty(61)+', '+derqty(62)+' ['+derunit(20)+']', $
                    derad+' '+derqty(63)+', '+derqty(64)+', '+derqty(65)+' ['+derunit(21)+']', $
                    derad+' '+derqty(66)+', '+derqty(67)+', '+derqty(68)+' ['+derunit(22)+']']
    endif else begin
        derqty=['tstar','tstar_errmin','tstar_errmax', $
                'ttotal','ttotal_errmin','ttotal_errmax', $
                'betastar','betastar_errmin','betastar_errmax', $
                'betagas','betagas_errmin','betagas_errmax', $
                'surfglobalstar','surfglobalstar_errmin','surfglobalstar_errmax', $
                'surfglobalgas','surfglobalgas_errmin','surfglobalgas_errmax', $
                'surfconstar','surfconstar_errmin','surfconstar_errmax', $
                'surfcongas','surfcongas_errmin','surfcongas_errmax', $

                'fcl','fcl_errmin','fcl_errmax', $
                'fgmc','fgmc_errmin','fgmc_errmax', $
                ;add: fdiffstar
                ;add: fdiffgas


                'zetastar','zetastar_errmin','zetastar_errmax', $
                'zetagas','zetagas_errmin','zetagas_errmax', $
                'rpeakstar','rpeakstar_errmin','rpeakstar_errmax', $
                'rpeakgas','rpeakgas_errmin','rpeakgas_errmax', $
                'vfb','vfb_errmin','vfb_errmax']
        derunit=['Myr','Myr','','','','','','','','','','','pc','pc','km s^-1']
        derad='Derived'
        derstrings=[derad+' '+derqty(0)+', '+derqty(1)+', '+derqty(2)+' ['+derunit(0)+']', $
                    derad+' '+derqty(3)+', '+derqty(4)+', '+derqty(5)+' ['+derunit(1)+']', $
                    derad+' '+derqty(6)+', '+derqty(7)+', '+derqty(8)+' ['+derunit(2)+']', $
                    derad+' '+derqty(9)+', '+derqty(10)+', '+derqty(11)+' ['+derunit(3)+']', $
                    derad+' '+derqty(12)+', '+derqty(13)+', '+derqty(14)+' ['+derunit(4)+']', $
                    derad+' '+derqty(15)+', '+derqty(16)+', '+derqty(17)+' ['+derunit(5)+']', $
                    derad+' '+derqty(18)+', '+derqty(19)+', '+derqty(20)+' ['+derunit(6)+']', $
                    derad+' '+derqty(21)+', '+derqty(22)+', '+derqty(23)+' ['+derunit(7)+']', $
                    derad+' '+derqty(24)+', '+derqty(25)+', '+derqty(26)+' ['+derunit(8)+']', $
                    derad+' '+derqty(27)+', '+derqty(28)+', '+derqty(29)+' ['+derunit(9)+']', $
                    derad+' '+derqty(30)+', '+derqty(31)+', '+derqty(32)+' ['+derunit(10)+']', $
                    derad+' '+derqty(33)+', '+derqty(34)+', '+derqty(35)+' ['+derunit(11)+']', $
                    derad+' '+derqty(36)+', '+derqty(37)+', '+derqty(38)+' ['+derunit(12)+']', $
                    derad+' '+derqty(39)+', '+derqty(40)+', '+derqty(41)+' ['+derunit(13)+']', $
                    derad+' '+derqty(42)+', '+derqty(43)+', '+derqty(44)+' ['+derunit(14)+']']
    endelse

    auxqty=['npeak_star','npeak_gas', $
            'lap_min']
    auxunit=['','pc']
    auxad='Derived'
    auxstrings=[auxad+' '+auxqty(0)+', '+auxqty(1)+' ['+auxunit(0)+']', $
                auxad+' '+auxqty(2)+' ['+auxunit(1)+']']
    nfit=n_elements(fitstrings)
    nder=n_elements(derstrings)
    if map_units gt 0 then begin
        next=n_elements(extstrings)
    endif else begin
        next=0
    endelse
    naux=n_elements(auxstrings)
    ntot=nfit+next+nder+naux
    nvarfit0=1
    nvarfit=3
    nvarext=3
    nvarder=3
    nvaraux=2
    nvaraux2=1
    nentries=nvarfit0+(nfit-1)*nvarfit+next*nvarext+nder*nvarder+naux*nvaraux-1
    colstrings=strarr(ntot)
    onecol='Column'
    multicol='Columns'
    for i=0,ntot-1 do begin
        if i eq 0 then begin
            nvar=nvarfit0
            nstart=i
        endif
        if i ge 1 && i lt nfit then begin
            nvar=nvarfit
            nstart=1+nvar*(i-1)
        endif
        if i ge nfit && i lt nfit+next && map_units gt 0 then begin
            nvar=nvarext
            nstart=1+(nfit-1)*nvarfit+nvar*(i-nfit)
        endif
        if i ge nfit+next && i lt nfit+next+nder then begin
            nvar=nvarder
            nstart=1+(nfit-1)*nvarfit+next*nvarext+nvar*(i-nfit-next)
        endif
        if i ge nfit+next+nder && i lt nfit+next+nder+naux-1 then begin
            nvar=nvaraux
            nstart=1+(nfit-1)*nvarfit+next*nvarext+nder*nvarder+nvar*(i-nfit-next-nder)
        endif
        if i eq nfit+next+nder+naux-1 then begin
            nvar=nvaraux2
            nstart=1+(nfit-1)*nvarfit+next*nvarext+nder*nvarder+(naux-1)*nvaraux
        endif
        if nvar gt 1 then colstrings(i)=multicol+' '+f_string(nstart+1,0)+'-'+f_string(nstart+nvar,0)+': ' else colstrings(i)=onecol+' '+f_string(nstart+1,0)+': '
    endfor

    openw,lun,outputdir+galaxy+'_tablerow.dat',/get_lun
    printf,lun,'# Best-fitting values and derived quantities for run '+galaxy+', generated with the Kruijssen & Longmore (2014) uncertainty principle code'
    printf,lun,'# IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how these numbers were calculated'
    printf,lun,'# IMPORTANT: all values represent log10(listed quantity)'
    for i=0,nfit-1 do printf,lun,'# '+colstrings(i)+fitstrings(i)
    if map_units gt 0 then begin
        for i=0,next-1 do printf,lun,'# '+colstrings(i+nfit)+extstrings(i)
    endif
    for i=0,nder-1 do printf,lun,'# '+colstrings(i+nfit+next)+derstrings(i)
    for i=0,naux-1 do printf,lun,'# '+colstrings(i+nfit+next+nder)+auxstrings(i)
    printf,lun,'# '+onecol+' '+f_string(nstart+2,0)+': Run ID (galaxy)'
    if map_units gt 0 then printf,lun,format='(f12.5,'+f_string(nentries-1,0)+'(3x,f12.5),3x,a24)',alog10([fit,ext,der,aux]+tiny),galaxy $
                        else printf,lun,format='(f12.5,'+f_string(nentries-1,0)+'(1x,f12.5),3x,a24)',alog10([fit,der,aux]+tiny),galaxy
    close,lun
    free_lun,lun

    varlen=20
    unitlen=40
    openw,lun,outputdir+galaxy+'_output.dat',/get_lun
    printf,lun,'########################################################################################################################'
    printf,lun,'#                                                                                                                      #'
    printf,lun,'#                                      FIT KL14 PRINCIPLE TO OBSERVED GALAXY MAPS                                      #'
    printf,lun,'# Best-fitting values and derived quantities generated with the Kruijssen & Longmore (2014) uncertainty principle code #'
    printf,lun,'#           IMPORTANT: see Paper II (Kruijssen et al. 2017) for details on how these numbers were calculated           #'
    printf,lun,'#                                                                                                                      #'
    printf,lun,'#                                                  BEGIN OUTPUT FILE                                                   #'
    printf,lun,'#                                                                                                                      #'
    printf,lun,'########################################################################################################################'
    printf,lun,''
    printf,lun,''
    printf,lun,'# FUNDAMENTAL QUANTITIES (obtained directly from the fitting process)'
    for i=0,(nfit-1)*3 do printf,lun,format='(a'+f_string(varlen,0)+',3x,f12.5,3x,a'+f_string(unitlen,0)+')',fitqty(i),alog10(fit(i)),'# log10['+fitqty(i)+'/('+fitunit((i+2)/3)+')]'
    printf,lun,''
    printf,lun,''
    if map_units gt 0 then begin
        printf,lun,'# EXTERNAL QUANTITIES (obtained through secondary analysis of the maps)'
        for i=0,next*3-1 do printf,lun,format='(a'+f_string(varlen,0)+',3x,f12.5,3x,a'+f_string(unitlen,0)+')',extqty(i),alog10(ext(i)),'# log10['+extqty(i)+'/('+extunit(i/3)+')]'
        printf,lun,''
        printf,lun,''
        printf,lun,'# DERIVED QUANTITIES (from fundamental and external quantities)'
    endif else printf,lun,'# DERIVED QUANTITIES (from fundamental quantities)'
    for i=0,nder*3-1 do printf,lun,format='(a'+f_string(varlen,0)+',3x,f12.5,3x,a'+f_string(unitlen,0)+')',derqty(i),alog10(der(i)),'# log10['+derqty(i)+'/('+derunit(i/3)+')]'
    printf,lun,''
    printf,lun,''
    printf,lun,'# AUXILIARY QUANTITIES (byproduct of the fitting process)'
    for i=0,naux*2-2 do printf,lun,format='(a'+f_string(varlen,0)+',3x,f12.5,3x,a'+f_string(unitlen,0)+')',auxqty(i),alog10(aux(i)),'# log10['+auxqty(i)+'/('+auxunit(i/2)+')]'
    printf,lun,''
    printf,lun,''
    printf,lun,'########################################################################################################################'
    printf,lun,'#                                                                                                                      #'
    printf,lun,'#                                                   END OUTPUT FILE                                                    #'
    printf,lun,'#                                                                                                                      #'
    printf,lun,'########################################################################################################################'
    close,lun ;close the logical unit
    free_lun,lun ;make the logical unit available again

    print,''
    print,'         Galaxy: '+galaxy
    for i=0,0 do print,'         '+fitstrings(i)+strtrim(fit(0))
    for i=1,nfit-1 do print,'         '+fitstrings(i)+strtrim(fit(3*i-2))+strtrim(fit(3*i-1))+strtrim(fit(3*i))
    if map_units gt 0 then begin
        for i=0,next-1 do print,'         '+extstrings(i)+strtrim(ext(3*i))+strtrim(ext(3*i+1))+strtrim(ext(3*i+2))
    endif
    for i=0,nder-1 do print,'         '+derstrings(i)+strtrim(der(3*i))+strtrim(der(3*i+1))+strtrim(der(3*i+2))
    for i=0,naux-2 do print,'         '+auxstrings(i)+strtrim(aux(2*i))+strtrim(aux(2*i+1))
    print,'         '+auxstrings(i)+strtrim(aux(2*i))
    print,''
endif


;;;;;;;;;;;;;;;;;;;;;;;;
;DELETE AUXILIARY FILES;
;;;;;;;;;;;;;;;;;;;;;;;;

if cleanup then begin
    deldirs=[peakdir,griddir,arrdir,maskeddir]
    for i=0,naperture-1 do deldirs=[deldirs,rundir+'res_'+res(i)+'pc'+path_sep()]
    ndel=n_elements(deldirs)
    eligibility=dblarr(ndel)
    for i=0,ndel-1 do begin
        nfile_tot=file_search(deldirs(i)+'*',count=count_all)
        nfile_eligible=file_search(deldirs(i)+'*.{clf,dat,fits,sav}',count=count_eligible)
        if count_all eq count_eligible then eligibility(i)=1
    endfor
    eligible=where(eligibility eq 1)
    if eligible(0) ne -1 then begin
        neligible=n_elements(eligible)
        print,' THE FOLLOWING DIRECTORIES WILL BE DELETED'
        for i=0,neligible-1 do print,'    '+deldirs(eligible(i))
    endif
    ineligible=where(eligibility eq 0)
    if ineligible(0) ne -1 then begin
        nineligible=n_elements(ineligible)
        print,' THE FOLLOWING DIRECTORIES CONTAIN FILES WITH EXTENSIONS OTHER THAN .{clf,dat,fits,sav} AND WILL NOT BE DELETED'
        for i=0,nineligible-1 do print,'    '+deldirs(ineligible(i))
    endif
    ndecimals=4
    code=10.^ndecimals*randomu(systime(1))
    if cleanup eq 2 then delchoice=round(code) else read,' TYPE '+f_string(code,0)+' TO PROCEED: ',delchoice
    if delchoice eq round(code) then begin
        for i=0,neligible-1 do spawn,'rm -fr '+deldirs(eligible(i))
        print, ' All directories listed for deletion have been removed'
    endif else begin
        print, ' Invalid code, no directories were deleted'
    endelse
endif


;;;;;;;;;;;;;;
;END ANALYSIS;
;;;;;;;;;;;;;;

for i=0,2 do begin
    beep
    wait,1.
endfor
endtime=systime(1)
duration=endtime-starttime
hours=fix(duration/3600.)
minutes=fix((duration/3600.-hours)*60.)
seconds=duration-3600.*hours-60.*minutes
print,' ==> finished KL14 tuning fork analysis, date/time is ',systime(0)
print,' ==> total analysis time for galaxy '+galaxy+' was '+f_string(hours,0)+'h'+f_string(minutes,0)+'m'+f_string(seconds,1)+'s'

if autoexit then print,' ==> IDL will now exit automatically'
journal
spawn,'mv '+outputdir+logfile+' '+outputdir+galaxy+'_'+logfile
if autoexit then exit

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                        ;
;             END MAIN ROUTINE tuningfork.pro            ;
;                                                        ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
