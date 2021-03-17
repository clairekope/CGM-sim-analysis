import yt
yt.enable_parallelism()

datasets = yt.load("DD????/DD????")
for ds in datasets.piter(dynamic=True, ):
  #if int(ds.basename[-4:])%5 == 0: #and int(ds.basename[-4:])>85:

     center = ds.quan(0.5, 'code_length')
     rs = ds.quan(3.5,'kpc')
     width = ds.quan(400,'kpc')

     sph = ds.sphere([center, center, center],(500,'kpc')) 
     rect = ds.region([center, center, center],
                      [center-rs/2, center-width/2, center-width/2],
                      [center+rs/2, center+width/2, center+width/2])
     
     p = yt.ProjectionPlot(ds,'x',
                           ['radial_velocity','tangential_velocity'],
                           width=width, data_source=rect, weight_field='ones')

     p.set_unit('radial_velocity','km/s')
     p.set_zlim('radial_velocity', -300, 300)
     p.set_log('radial_velocity', 'log', linthresh=5)
     p.set_cmap('radial_velocity', 'coolwarm')

     p.set_cmap('tangential_velocity', 'BuPu')
     p.set_unit('tangential_velocity','km/s')
     p.set_zlim('tangential_velocity', 1e-2, 1e3)

     p.annotate_timestamp()
     p.save() 


     p = yt.ProjectionPlot(ds,'x','density',
                           width=width/2, data_source=rect, weight_field='ones')

     # if ("all","particle_mass") in ds.derived_field_list: 
     #     p.annotate_particles(width=width) 

     p.set_zlim('density',1e-31,1e-24)
     p.set_colorbar_label('density','Average Density (g cm$^{-3}$)')

     p.annotate_timestamp()
     p.save() 
