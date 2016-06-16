;----------------------------------------------------------------------------------------
pro mask_ds9_file_convert, image_path, ds9_region_inpath, ds9_region_outpath
;----------------------------------------------------------------------------------------
; *** aim_ds9_file_convert uses ds9 to convert ds9 region files to image co-ordinates
;----------------------------------------------------------------------------------------
; *** To run, the aim_tool  requires the:
; *** 1) The IDL Astronomy User's Library http://idlastro.gsfc.nasa.gov/
;---input--------------------------------------------------------------------------------
; *** image_path        = filepath of the .fits image that the ds9 region file
; ***                     corresponds to
; *** ds9_region_inpath = filepath of the ds9 region to be read and converted
;---output-------------------------------------------------------------------------------
; *** ds9_region_outpath = filepath of the ds9 region file to be written after conversion
;----------------------------------------------------------------------------------------
  compile_opt idl2 ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults

  spawn_command = 'ds9 ' + image_path + ' -regions load ' + ds9_region_inpath +  $
                  ' -regions system image -scale log  -regions save ' + $
                  ds9_region_outpath + ' -exit'

  spawn, spawn_command

end
