import yt
yt.enable_parallelism()

datasets = yt.load("DD????/DD????")
for ds in datasets.piter():
#  if int(ds.basename[-4:])%5 == 0: #and int(ds.basename[-4:])>85:

     center = ds.quan(0.5, 'code_length')
     rs = ds.quan(3.5,'kpc')
     width = ds.quan(400,'kpc')

     sph = ds.sphere([center, center, center],(500,'kpc')) 
     rect = ds.region([center, center, center],
                      [center-rs/2, center-width/2, center-width/2],
                      [center+rs/2, center+width/2, center+width/2])
     

     ph = yt.PhasePlot(sph,'density','temperature','cell_mass').set_unit('cell_mass','Msun')
     ph.set_xlim(5e-31,5e-21)
     ph.set_ylim(10,1e8)
     ph.set_zlim('cell_mass',1e-2, 1e7);
     ph.save()

     p = yt.ProjectionPlot(ds,'x','density',width=width,
                           data_source=rect, weight_field='ones')
     #if ("all","particle_mass") in ds.derived_field_list: 
     #    p.annotate_particles(width=width) 
     p.set_zlim('density',1e-31,1e-24)
     p.set_colorbar_label('density','Average Density (g cm$^{-3}$)')
     p.annotate_timestamp()
     p.save() 

     p = yt.ProjectionPlot(ds,'x','temperature', width=width,
                           data_source=rect, weight_field='ones')
     p.set_zlim('temperature',1e3,1e8)
     p.set_cmap('temperature','plasma')
     p.set_colorbar_label('temperature','Average Temperature (K)')
     p.annotate_timestamp()
     p.save()


     # p = yt.ProjectionPlot(ds,'x','temperature', width=width, data_source=sph, method='mip')
     # #p.set_zlim('temperature',5e4,1e7)
     # p.set_cmap('temperature','plasma')
     # #p.set_colorbar_minorticks('temperature','on')
     # p.annotate_timestamp()
     # p.save("MIP_{}".format(ds.basename))
