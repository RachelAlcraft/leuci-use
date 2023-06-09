"""
Created 21/5/2023 - Rachel Alcraft
Updates:

--------------------------------------
This script takes 
- a pdb file
- 3 atoms
- a width
- an interpolation method
And produces and html file with an inmage of the density using plotly

"""
#### USER INPUTS ###
pdb_code = "7uly"
central_atom = "A:3@CE1.A"
linear_atom = "A:3@CZ.A"
planar_atom = "A:3@CE2.A"
interpolation = "bspline"
width = 9 #Angstrom
samples = 50
depth_samples = 10
file_output2d = "slice-2d.html"
file_output3d = "slice-3d.html"
#### ----------- ###

### CODE ###
from leuci_xyz import vectorthree as v3
from leuci_map import mapsmanager as mman
from leuci_map import mapfunctions as mfun
from leuci_map import mapplotter as mpl

# find relative path for data
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent.parent.parent )+ "/pdbdata/"
RESDIR = str(Path(__file__).resolve().parent.parent.parent )+ "/results/" 
mman.MapsManager().set_dir(DATADIR)

# downloand/upload into memory the pdb and ccp4 data (0 means skip if not there, 1 means in this thread, 2 means in another thread and don't wait)
ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)

mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation)

cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(central_atom))
ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(linear_atom))
pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(planar_atom))

# 2d plot (s)
filename = RESDIR + file_output2d
mplot = mpl.MapPlotter(mf, filename,interpolation,samples,width,cc,ll,pp,(1,1))
mplot.add_plot_slice("density",2,-1,0.9, 0.9,True,"BW","density",(1,1)) 
mplot.make_plot_slices(log_level=1)   

# 3d plot
filename = RESDIR + file_output3d
mplot = mpl.MapPlotter(mf, filename,interpolation,samples,width,cc,ll,pp,(1,1))
vals,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples)
mplot.make_plot_slice_3d("density",vals,coords,min_percent=0.9, max_percent=0.9,hue="GBR",centre=True,title="Leucippus Plot 3d")

