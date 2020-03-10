;----------------------------------------------------------------------------------------
function get_heisenberg_output_var_struct, map_units
;----------------------------------------------------------------------------------------
; a function that returns a structure containing the expected output variables from
; running Heisenberg, along with whether the variable has associated errors, the units
; and the name of the variable for plotting given the units of the input maps
;--(input)-------------------------------------------------------------------------------
; ***  map_units  = a number from 0-3 specifying the units of the input maps
;--(output)------------------------------------------------------------------------------
; *** heisenberg_output_var_struct = a structure of the expected output variables from
;                                running Heisenberg
;-----------------------------------------------------------------------------------------
  compile_opt idl2, strictarrsubs ; enforce strict array indexing, i.e. only use of [] and not () to index and use 32bit integers rather than 16bit as defaults ; give error when subscripting one array using another array as the source of array indices that has out-of-range indices, rather than clipping into range


  ; ********************************************
  ; fit quantities (always present in output)
  heisenberg_output_var_struct =[ $
    {name:'redchi2', errors:0, plot_title_name:'redchi2',units:''}, $
    {name:'tgas', errors:1, plot_title_name:'!8t!6!Dgas!N !6',units:'Myr'}, $
    {name:'tover', errors:1, plot_title_name:'!8t!6!Dover!N !6',units:'Myr'}, $
    {name:'lambda', errors:1, plot_title_name:'!7k !6',units:'pc'}]

  ; ********************************************
  ; external quantities (differ based on map_units)
  if map_units eq 1 then begin
    heisenberg_output_var_struct =[ heisenberg_output_var_struct, $
      {name:'sfr_galaxy', errors:1, plot_title_name:'sfr_galaxy',units:''}, $
      {name:'mgas_galaxy', errors:1, plot_title_name:'mgas_galaxy',units:''}, $
      {name:'surfsfr', errors:1, plot_title_name:'surfsfr',units:''}, $
      {name:'surfgas', errors:1, plot_title_name:'surfgas',units:''}]

  endif else if map_units eq 2 then begin
    heisenberg_output_var_struct =[ heisenberg_output_var_struct, $
      {name:'mgas_galaxy1', errors:1, plot_title_name:'mgas_galaxy1',units:''}, $
      {name:'mgas_galaxy2', errors:1, plot_title_name:'mgas_galaxy2',units:''}, $
      {name:'surfgas1', errors:1, plot_title_name:'surfgas1',units:''}, $
      {name:'surfgas2', errors:1, plot_title_name:'surfgas2',units:''}]

  endif else if map_units eq 3 then begin
     heisenberg_output_var_struct =[ heisenberg_output_var_struct, $
     {name:'sfr_galaxy1', errors:1, plot_title_name:'sfr_galaxy1',units:''}, $
     {name:'sfr_galaxy2', errors:1, plot_title_name:'sfr_galaxy2',units:''}, $
     {name:'surfsfr1', errors:1, plot_title_name:'surfsfr1',units:''}, $
     {name:'surfsfr2', errors:1, plot_title_name:'surfsfr2',units:''}]

  endif



  ; ********************************************
  ;  derived quantities (part 1: variables always present in output)
  heisenberg_output_var_struct =[heisenberg_output_var_struct, $
    {name:'tstar', errors:1, plot_title_name:'tstar',units:''}, $
    {name:'ttotal', errors:1, plot_title_name:'ttotal',units:''}, $
    {name:'betastar', errors:1, plot_title_name:'betastar',units:''}, $
    {name:'betagas', errors:1, plot_title_name:'betagas',units:''}, $
    {name:'surfglobalstar', errors:1, plot_title_name:'surfglobalstar',units:''}, $
    {name:'surfglobalgas', errors:1, plot_title_name:'surfglobalgas',units:''}, $
    {name:'surfconstar', errors:1, plot_title_name:'surfconstar',units:''}, $
    {name:'surfcongas', errors:1, plot_title_name:'surfcongas',units:''}, $
    {name:'fcl', errors:1, plot_title_name:'fcl',units:''}, $
    {name:'fgmc', errors:1, plot_title_name:'fgmc',units:''}, $
    {name:'qconstar',errors:1, plot_title_name:'!6Q!Dcon,star!N',units:''}, $
    {name:'qcongas',errors:1, plot_title_name:'!6Q!Dcon,gas!N',units:''}, $
    {name:'etastar',errors:1, plot_title_name:'!7g!6!Dstar!N',units:''}, $
    {name:'etagas',errors:1, plot_title_name:'!7g!6!Dgas!N',units:''}, $
    {name:'qzetastar',errors:1, plot_title_name:'!6Q!D!7f!6,star!N',units:''}, $
    {name:'qzetagas',errors:1, plot_title_name:'!6Q!D!7f!6,gas!N',units:''}, $
    {name:'zetastar', errors:1, plot_title_name:'zetastar',units:''}, $
    {name:'zetagas', errors:1, plot_title_name:'zetagas',units:''}, $
    {name:'rpeakstar', errors:1, plot_title_name:'rpeakstar',units:''}, $
    {name:'rpeakgas', errors:1, plot_title_name:'rpeakgas',units:''}, $
    {name:'vfb', errors:1, plot_title_name:'vfb',units:''}, $
    {name:'vfbr', errors:1, plot_title_name:'vfbr',units:''}]

  ; ********************************************
  ; get derived quantities (part 2: variables only present in output if map_units gt 0)
  if map_units gt 0 then begin
    heisenberg_output_var_struct =[ heisenberg_output_var_struct, $
      {name:'tdepl', errors:1, plot_title_name:'tdepl',units:''}, $
      {name:'esf', errors:1, plot_title_name:'esf',units:''}, $
      {name:'mdotsf', errors:1, plot_title_name:'mdotsf',units:''}, $
      {name:'mdotfb', errors:1, plot_title_name:'mdotfb',units:''}, $
      {name:'etainst', errors:1, plot_title_name:'etainst',units:''}, $
      {name:'etaavg', errors:1, plot_title_name:'etaavg',units:''}, $
      {name:'pzero', errors:1, plot_title_name:'pzero',units:''}, $
      {name:'chie', errors:1, plot_title_name:'chie',units:''}, $
      {name:'chier', errors:1, plot_title_name:'chier',units:''}, $
      {name:'chip', errors:1, plot_title_name:'chip',units:''}, $
      {name:'chipr', errors:1, plot_title_name:'chip',units:''}]

  endif

  ; ********************************************
  ; get auxilliary quantities (always present in output)
  heisenberg_output_var_struct =[ heisenberg_output_var_struct, $
    {name:'npeak_star', errors:0, plot_title_name: 'npeak_star',units:''}, $
    {name:'npeak_gas', errors:0, plot_title_name: 'npeak_gas',units:''}, $
    {name:'lap_min', errors:0, plot_title_name: 'lap_min',units:''}]

  return, heisenberg_output_var_struct

end
