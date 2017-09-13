function get_kl14_output_var_struct, map_units
  ; ********************************************
  ; fit quantities (always present in output)
   kl14_output_var_struct =[ $
    {name:'redchi2', errors:0, plot_title:'redchi2'}, $
    {name:'tgas', errors:1, plot_title:'!8t!6!Dgas!N !6[Myr]'}, $
    {name:'tover', errors:1, plot_title:'!8t!6!Dover!N !6[Myr]'}, $
    {name:'lambda', errors:1, plot_title:'!7k !6[pc]'}]

  ; ********************************************
  ; external quantities (differ based on map_units)
  if map_units eq 1 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'sfr_galaxy', errors:1, plot_title:'sfr_galaxy'}, $
      {name:'mgas_galaxy', errors:1, plot_title:'mgas_galaxy'}, $
      {name:'surfsfr', errors:1, plot_title:'surfsfr'}, $
      {name:'surfgas', errors:1, plot_title:'surfgas'}]

  endif else if map_units eq 2 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'mgas_galaxy1', errors:1, plot_title:'mgas_galaxy1'}, $
      {name:'mgas_galaxy2', errors:1, plot_title:'mgas_galaxy2'}, $
      {name:'surfgas1', errors:1, plot_title:'surfgas1'}, $
      {name:'surfgas2', errors:1, plot_title:'surfgas2'}]

  endif else if map_units eq 3 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
     {name:'sfr_galaxy1', errors:1, plot_title:'sfr_galaxy1'}, $
     {name:'sfr_galaxy2', errors:1, plot_title:'sfr_galaxy2'}, $
     {name:'surfsfr1', errors:1, plot_title:'surfsfr1'}, $
     {name:'surfsfr2', errors:1, plot_title:'surfsfr2'}]

  endif



  ; ********************************************
  ;  derived quantities (part 1: variables always present in output)
   kl14_output_var_struct =[kl14_output_var_struct, $
    {name:'tstar', errors:1, plot_title:'tstar'}, $
    {name:'ttotal', errors:1, plot_title:'ttotal'}, $
    {name:'betastar', errors:1, plot_title:'betastar'}, $
    {name:'betagas', errors:1, plot_title:'betagas'}, $
    {name:'surfglobalstar', errors:1, plot_title:'surfglobalstar'}, $
    {name:'surfglobalgas', errors:1, plot_title:'surfglobalgas'}, $
    {name:'surfconstar', errors:1, plot_title:'surfconstar'}, $
    {name:'surfcongas', errors:1, plot_title:'surfcongas'}, $
    {name:'fcl', errors:1, plot_title:'fcl'}, $
    {name:'fgmc', errors:1, plot_title:'fgmc'}, $
    {name:'zetastar', errors:1, plot_title:'zetastar'}, $
    {name:'zetagas', errors:1, plot_title:'zetagas'}, $
    {name:'rpeakstar', errors:1, plot_title:'rpeakstar'}, $
    {name:'rpeakgas', errors:1, plot_title:'rpeakgas'}, $
    {name:'vfb', errors:1, plot_title:'vfb'}]

  ; ********************************************
  ; get derived quantities (part 2: variables only present in output if map_units gt 0)
  if map_units gt 0 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'tdepl', errors:1, plot_title:'tdepl'}, $
      {name:'esf', errors:1, plot_title:'esf'}, $
      {name:'mdotsf', errors:1, plot_title:'mdotsf'}, $
      {name:'mdotfb', errors:1, plot_title:'mdotfb'}, $
      {name:'etainst', errors:1, plot_title:'etainst'}, $
      {name:'etaavg', errors:1, plot_title:'etaavg'}, $
      {name:'chie', errors:1, plot_title:'chie'}, $
      {name:'chip', errors:1, plot_title:'chip'}]

  endif

  ; ********************************************
  ; get auxilliary quantities (always present in output)
   kl14_output_var_struct =[ kl14_output_var_struct, $
    {name:'npeak_star', errors:0, plot_title: 'npeak_star'}, $
    {name:'npeak_gas', errors:0, plot_title: 'npeak_gas'}, $
    {name:'lap_min', errors:0, plot_title: 'lap_min'}]

return, kl14_output_var_struct

end
