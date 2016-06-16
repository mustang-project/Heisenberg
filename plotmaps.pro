;PLOT MAPS FOR DIFFERENT RESOLUTIONS


galaxies=['regrid','res_50pc','res_100pc','res_200pc','res_400pc','res_800pc']
hafiles='galaxies/m33/'+galaxies+'/m33_ha'
cofiles='galaxies/m33/'+galaxies+'/m33_nro.mom0.velmask'

plotfits_files,hafiles,'m33','ha',[6.5,5.95,6.55,6.16,5.5,5.2],0.*[1,1,1,1,1,1],1,log=1
plotfits_files,cofiles,'m33','co',[4.98,6.05,5.14,5.59,4.95,5.18],0.*[1,1,1,1,1,1],3,log=1




end



