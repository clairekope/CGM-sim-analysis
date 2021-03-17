import yt
yt.enable_parallelism()

# datasets = yt.load("DD????/DD????")
# ds = datasets[-1]
# if True:
for ds in datasets.piter(dynamic=True, ):

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

     ph = yt.PhasePlot(sph,'radius','cooling_time','cell_mass')
     ph.set_unit('cell_mass','Msun')
     ph.set_unit('radius','kpc')
     ph.set_unit('cooling_time','Gyr')
     ph.save()


     ph = yt.PhasePlot(sph,'radius','temperature','cell_mass')
     ph.set_unit('cell_mass','Msun')
     ph.set_unit('radius','kpc')
     ph.save()


     ph = yt.PhasePlot(sph,'radius','entropy','cell_mass')
     ph.set_unit('cell_mass','Msun')
     ph.set_unit('radius','kpc')
     ph.save()
