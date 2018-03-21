function get_kl14_output_var_struct, map_units
  ; ********************************************
  ; fit quantities (always present in output)
   kl14_output_var_struct =[ $
    {name:'redchi2', errors:0, plot_title_name:'redchi2',units:''}, $
    {name:'tgas', errors:1, plot_title_name:'!8t!6!Dgas!N !6',units:'Myr'}, $
    {name:'tover', errors:1, plot_title_name:'!8t!6!Dover!N !6',units:'Myr'}, $
    {name:'lambda', errors:1, plot_title_name:'!7k !6',units:'pc'}]

  ; ********************************************
  ; external quantities (differ based on map_units)
  if map_units eq 1 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'sfr_galaxy', errors:1, plot_title_name:'sfr_galaxy',units:''}, $
      {name:'mgas_galaxy', errors:1, plot_title_name:'mgas_galaxy',units:''}, $
      {name:'surfsfr', errors:1, plot_title_name:'surfsfr',units:''}, $
      {name:'surfgas', errors:1, plot_title_name:'surfgas',units:''}]

  endif else if map_units eq 2 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'mgas_galaxy1', errors:1, plot_title_name:'mgas_galaxy1',units:''}, $
      {name:'mgas_galaxy2', errors:1, plot_title_name:'mgas_galaxy2',units:''}, $
      {name:'surfgas1', errors:1, plot_title_name:'surfgas1',units:''}, $
      {name:'surfgas2', errors:1, plot_title_name:'surfgas2',units:''}]

  endif else if map_units eq 3 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
     {name:'sfr_galaxy1', errors:1, plot_title_name:'sfr_galaxy1',units:''}, $
     {name:'sfr_galaxy2', errors:1, plot_title_name:'sfr_galaxy2',units:''}, $
     {name:'surfsfr1', errors:1, plot_title_name:'surfsfr1',units:''}, $
     {name:'surfsfr2', errors:1, plot_title_name:'surfsfr2',units:''}]

  endif



  ; ********************************************
  ;  derived quantities (part 1: variables always present in output)
   kl14_output_var_struct =[kl14_output_var_struct, $
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
    {name:'vfb', errors:1, plot_title_name:'vfb',units:''}]

  ; ********************************************
  ; get derived quantities (part 2: variables only present in output if map_units gt 0)
  if map_units gt 0 then begin
     kl14_output_var_struct =[ kl14_output_var_struct, $
      {name:'tdepl', errors:1, plot_title_name:'tdepl',units:''}, $
      {name:'esf', errors:1, plot_title_name:'esf',units:''}, $
      {name:'mdotsf', errors:1, plot_title_name:'mdotsf',units:''}, $
      {name:'mdotfb', errors:1, plot_title_name:'mdotfb',units:''}, $
      {name:'etainst', errors:1, plot_title_name:'etainst',units:''}, $
      {name:'etaavg', errors:1, plot_title_name:'etaavg',units:''}, $
      {name:'chie', errors:1, plot_title_name:'chie',units:''}, $
      {name:'chip', errors:1, plot_title_name:'chip',units:''}]

  endif

  ; ********************************************
  ; get auxilliary quantities (always present in output)
   kl14_output_var_struct =[ kl14_output_var_struct, $
    {name:'npeak_star', errors:0, plot_title_name: 'npeak_star',units:''}, $
    {name:'npeak_gas', errors:0, plot_title_name: 'npeak_gas',units:''}, $
    {name:'lap_min', errors:0, plot_title_name: 'lap_min',units:''}]

return, kl14_output_var_struct

end
