

######### TEST THE DATAFRAME ##############################

"""
Search criteria
method/refinement resolution < 1
Disulfide bond count >=2
Exp method x-ray
no RNA
"""
from leuci_xyz import vectorthree as v3
from leuci_map import mapsmanager as mman
from leuci_map import mapfunctions as mfun
from leuci_map import mapplothelp as mph
from leuci_geo import pdbloader as pl
from leuci_geo import pdbgeometry as pg
# find relative path for data
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent.parent.parent )+ "/pdbdata/"
RESDIR = str(Path(__file__).resolve().parent.parent.parent )+ "/results/" 
mman.MapsManager().set_dir(DATADIR)

possible_stages = "check,run,each2d,each3d,all2d,all3d,naybs"
#these_stages = "check"
#these_stages = "run,csv"
these_stages = "run,csv,each2d,all2d,all3d,naybs"
pdb_set = "ec" #xr,em or ec


pdbs = []
# these are all pdbs, followed by thre redcued to ccp4 available - done manually
###########################################################################
##pdbs.extend(["1AB1", "1AHO", "1CBN", "1DY5", "1EJG", "1ETL", "1ETM", "1ETN", "1F94", "1FN8", "1FY4", "1FY5", "1G4I", "1G66", "1G6X"])
##pdbs.extend(["1GDN", "1GDQ", "1GQV", "1GVK", "1HJ9", "1HJE", "1IC6", "1IEE", "1JXT", "1JXU", "1JXW", "1JXX", "1JXY", "1K5C", "1KTH"])
##pdbs.extend(["1L9L", "1MC2", "1OK0", "1P9G", "1PQ5", "1PQ7", "1SSX", "1V0L", "1V6P", "1VB0", "1VL9", "1X8P", "1X8Q", "1XVO", "1YWA", "1YWB", "1ZLB", "2AYW", "2B97"])
##pdbs.extend(["2BZZ", "2FMA", "2G58", "2H5C", "2H5D", "2IXT", "2NLS", "2PNE", "2PWA", "2V8B", "2VB1", "2VHK", "2VHR", "2VI3", "2VU6", "2XJH", "2XJP", "2XTT", "2XU3", "3AGN", "3AGO", "3C78"])
##pdbs.extend(["3D43", "3DW3", "3DWE", "3HGP", "3I2Y", "3I30", "3I37", "3LZT", "3M5Q", "3MFJ", "3MI4", "3NIR", "3ODV", "3PSM", "3Q8J", "3QPA", "3QPC", "3U7T", "3VLA", "3W7Y", "3WGE", "3WGX"])
##pdbs.extend(["3WL2", "3WOU", "3X2H", "3X2L", "3X2M", "3X2P", "4A7U", "4BCT", "4E3Y", "4F18", "4F19", "4F1U", "4F1V", "4HGU", "4I8G", "4I8H", "4I8J", "4I8K", "4I8L", "4LFS", "4LZT", "4M7G"])
##pdbs.extend(["4NDS", "4NSV", "4R5R", "4UNU", "4UYR", "4WKA", "4XDX", "4XOJ", "4Y9W", "4YEO", "4ZM7", "5A71", "5AVD", "5AVG", "5DJ7", "5DK1", "5DKM", "5E7W", "5E9N", "5HMV", "5HQI", "5I5B"])
##pdbs.extend(["5II6", "5KWM", "5KXV", "5MB5", "5MN1", "5MNB", "5MNC", "5MNF", "5MNG", "5MNH", "5MNK", "5MNM", "5MNN","5MNO", "5MON","5MOP", "5MOQ"])
##pdbs.extend(["5MOR", "5MOS", "5O0U", "5O2X", "5RBW","5RC2", "5RCB", "5SBQ", "5U3A", "5X9L", "5X9M", "6CNW", "6E6O","6EQE", "6ETK", "6ETL"])
##pdbs.extend(["6ETM","6ETN", "6F1O", "6RGP", "6RHH", "6RHU", "6RHX", "6RI6", "6RI8", "6RII", "6RYG","6SRY", "6SY3","6SYE", "6TN1", "6YIV"]) 
##pdbs.extend(["6YIW", "6ZSY", "7AEY", "7AF2", "7AVE", "7BCU", "7LTD","7LTI", "7LTV", "7MBO", "7OL5", "7P4R", "7P6M", "7Q5G", "7YRK"])
###########################################################################
if pdb_set == "xr":
    pdbs.extend(['1ab1', '1aho', '1ejg', '1etl', '1etm', '1etn', '1f94', '1fn8', '1fy4', '1fy5', '1g4i', '1g6x'])
    pdbs.extend(['1gdn', '1gdq', '1gqv', '1gvk', '1hj9', '1iee', '1jxt', '1jxu', '1jxw', '1jxx', '1jxy', '1k5c', '1kth'])
    pdbs.extend(['1l9l', '1mc2', '1p9g', '1pq5', '1pq7', '1v0l', '1vl9', '1x8p', '1x8q', '1xvo', '1ywa', '1ywb', '1zlb', '2ayw', '2b97'])
    pdbs.extend(['2bzz', '2fma', '2ixt', '2nls', '2pne', '2pwa', '2vb1', '2vhk', '2vhr', '2vi3', '2vu6', '2xjh', '2xjp', '2xtt', '2xu3', '3agn', '3ago', '3c78'])
    pdbs.extend(['3d43', '3dw3', '3dwe', '3hgp', '3i2y', '3i30', '3i37', '3lzt', '3m5q', '3mfj', '3mi4', '3nir', '3odv', '3psm', '3q8j', '3qpa', '3qpc', '3u7t', '3vla', '3w7y', '3wge', '3wgx'])
    pdbs.extend(['3wl2', '3wou', '3x2h', '3x2l', '3x2m', '3x2p', '4a7u', '4bct', '4e3y', '4f18', '4f19', '4f1u', '4f1v', '4hgu', '4i8g', '4i8h', '4i8j', '4i8k', '4i8l', '4lfs', '4lzt', '4m7g'])
    pdbs.extend(['4nds', '4nsv', '4r5r', '4unu', '4uyr', '4wka', '4xdx', '4xoj', '4y9w', '4yeo', '4zm7', '5a71', '5avd', '5avg', '5dj7', '5dk1', '5dkm', '5e7w', '5e9n', '5hqi', '5i5b'])
    pdbs.extend(['5ii6', '5kwm', '5kxv', '5mb5', '5mn1', '5mnb', '5mnc', '5mnf', '5mng', '5mnh', '5mnk', '5mnm', '5mnn', '5mno', '5mon', '5mop', '5moq'])
    pdbs.extend(['5mor', '5mos', '5o0u', '5o2x', '5rbw', '5rc2', '5rcb', '5sbq', '5u3a', '5x9l', '5x9m', '6cnw', '6e6o', '6eqe', '6etk', '6etl'])
    pdbs.extend(['6etm', '6etn', '6f1o', '6rgp', '6rhh', '6rhu', '6rhx', '6ri6', '6ri8', '6rii', '6ryg', '6sy3', '6sye', '6yiv'])
    pdbs.extend(['6yiw', '6zsy', '7bcu', '7ltd', '7lti', '7ltv', '7mbo', '7ol5', '7p6m', '7q5g', '7yrk'])
###########################################################################
##### Electron microscopy structures ######
## res < 2, no rna restriction
##pdbs.extend(["7A5V", "7FIX", "7JJP", "7JKC", "7MBX", "7MIJ", "7QEQ", "7UUR"])
elif pdb_set == "em":
    pdbs.extend(['7a5v', '7jjp', '7jkc', '7mbx', '7mij', '7qeq', '7uur'])
###########################################################################
##### Electron crystallography structures ######
## res < 2, no rna restriction
##pdbs.extend(["5I9S", "5K7O", "5K7R", "5K7S", "6CL7", "6H3B", "6LAV", "6LAW", "6PKP", "6PKQ", "6PKR", "6PKT", "6S2N", "6V8R", "7JSY", "7MRP", "7SKW", "7SKX", "7SW1", "7SW2", "7SW5", "7SW6", "7SW8", "7ULY", "8E53"])
elif pdb_set == "ec":
    pdbs.extend(['5i9s', '5k7o', '5k7r', '5k7s', '6cl7', '6lav', '6law', '6pkp', '6pkq', '6pkr', '6pkt', '6s2n', '6v8r', '7jsy', '7mrp', '7skw', '7skx', '7sw1', '7sw2', '7sw5', '7sw6', '7sw8', '7uly'])

###########################################################################


if "check" in these_stages:
    ccp4_pdbs = []
    print(len(pdbs),pdbs)
    for i in range(len(pdbs)):
        pdbs[i] = pdbs[i].lower()
        print("Loading ccp4",pdbs[i])
        try:
            ml = mman.MapsManager().get_or_create(pdbs[i],file=1,header=1,values=1)
            if not ml.success():
                print(pdbs[i], "has not loaded ccp4 succesfully")
            else:
                ccp4_pdbs.append(pdbs[i])
            mman.MapsManager().clear()
        except Exception as e:
            print("!!!",pdbs[i],e)

    print("########################################################################")
    print("These pdbs have electron density:", len(ccp4_pdbs))
    print(ccp4_pdbs)
    print("########################################################################")
    pdbs = []
    for pdb in ccp4_pdbs:
        pdbs.append(pdb)

if "run" in these_stages:
    # Make the geometry
    pobjs = []
    print("---")
    for pdb in pdbs:
        print("Loading",pdb)
        pla = pl.PdbLoader(pdb,DATADIR,cif=True)    
        po = pla.load_pdb()    
        pobjs.append(po)
        
    gm = pg.GeometryMaker(pobjs)

    print("----- geos -----")
    geos = ['SG:{SG@1}']
    df = gm.calculateGeometry(geos,log=0)
    print(df)

    print("----- geos -----")
    geos = ['SG:{SG@1}','SG:{(N),(O)}']
    df = gm.calculateGeometry(geos,log=0)
    print(df)

    print("----- geos -----")
    geos = ['SG:{SG@1}[dis|<3.5]','SG:{(N),(O)}[dis|<3.5]','SG:{(N),(O)}:{SG@1}']
    geos_r = ['SG:{SG@1}[dis|<3.5]','SG:{(N),(O)&1}[dis|<3.5]','SG:{(N),(O)&1}:{SG@1}'] #nearest o,n not the same residue
    geos_n = ['SG:{SG@1}[dis|<3.5]','SG:N[dis|<3.5]','SG:N:{SG@1}'] #its own N is in hb distance
    
    df = gm.calculateGeometry(geos,log=0)
    df_r = gm.calculateGeometry(geos_r,log=0)
    df_n = gm.calculateGeometry(geos_n,log=0)        
    dfs = []
    print("----- filter -----")    
    df = df.loc[(df['occ_SG:{(N),(O)}:{SG@1}'] == 1) ]        
    print("Occ=1",df)
    df2 = df.loc[(df['SG:{(N),(O)}:{SG@1}'] < 9) &(df['SG:{(N),(O)}:{SG@1}'] >= 0)]# (df['bf_SG:{(N),(O)}:{SG@1}'] < 10)  & (df['occ_SG:{(N),(O)}:{SG@1}'] == 1) ]        
    df3 = df.loc[(df['SG:{(N),(O)}:{SG@1}'] < 28) &(df['SG:{(N),(O)}:{SG@1}'] > 25)]# (df['bf_SG:{(N),(O)}:{SG@1}'] < 10)  & (df['occ_SG:{(N),(O)}:{SG@1}'] == 1) ]        
    df4 = df.loc[(df['SG:{(N),(O)}:{SG@1}'] < 66) &(df['SG:{(N),(O)}:{SG@1}'] > 64)]# (df['bf_SG:{(N),(O)}:{SG@1}'] < 10)  & (df['occ_SG:{(N),(O)}:{SG@1}'] == 1) ]        

    df2_r = df_r.loc[(df_r['SG:{(N),(O)&1}:{SG@1}'] < 9) &(df_r['SG:{(N),(O)&1}:{SG@1}'] >= 0)]
    df3_r = df_r.loc[(df_r['SG:{(N),(O)&1}:{SG@1}'] < 28) &(df_r['SG:{(N),(O)&1}:{SG@1}'] > 25)]
    df4_r = df_r.loc[(df_r['SG:{(N),(O)&1}:{SG@1}'] < 66) &(df_r['SG:{(N),(O)&1}:{SG@1}'] > 64)]

    df2_n = df_n.loc[(df_n['SG:N:{SG@1}'] < 9) &(df_n['SG:N:{SG@1}'] >= 0)]
    df3_n = df_n.loc[(df_n['SG:N:{SG@1}'] < 28) &(df_n['SG:N:{SG@1}'] > 25)]
    df4_n = df_n.loc[(df_n['SG:N:{SG@1}'] < 66) &(df_n['SG:N:{SG@1}'] > 64)]
            
    dfs.append((df2_r,"line_r",'info_SG:{(N),(O)&1}:{SG@1}'))
    dfs.append((df4_r,"65_r",'info_SG:{(N),(O)&1}:{SG@1}'))
    dfs.append((df3_r,"26_r",'info_SG:{(N),(O)&1}:{SG@1}'))

    dfs.append((df2_n,"line_n",'info_SG:N:{SG@1}'))
    dfs.append((df4_n,"65_n",'info_SG:N:{SG@1}'))
    dfs.append((df3_n,"26_n",'info_SG:N:{SG@1}'))
    
    dfs.append((df2,"line",'info_SG:{(N),(O)}:{SG@1}'))        
    dfs.append((df4,"65",'info_SG:{(N),(O)}:{SG@1}'))        
    dfs.append((df3,"26",'info_SG:{(N),(O)}:{SG@1}'))    
    
            
    for dfx, tag,atom_key in dfs:
        print("---------------",tag,"--------------------------------")
        print(dfx)
        if "csv" in these_stages:
            file_outputcsv = "sg_csv_"+tag+"_"+pdb_set
            filecsvname = RESDIR + "csv/" + file_outputcsv + ".csv"
            dfx.to_csv(filecsvname)
            print("Saved to", filecsvname)

        coords_list = []
        i_count = 0
        df_len = len(dfx.index)
        for i, row in dfx.iterrows():
            pdb_code = row['pdb_code']
            atoms = row[atom_key]    
            print(i_count, "/",df_len,i,atoms)
            i_count += 1
            clp = atoms.split("(")    
            cen = clp[1][:-1]
            lin = clp[2][:-1]
            pla = clp[3][:-1]
            cens = cen.split("|")
            lins = lin.split("|")
            plas = pla.split("|")
            #central_atom = "A:707@C.A"
            cen_str = cens[0]+":"+cens[2]+"@"+cens[3]+".A"
            lin_str = lins[0]+":"+lins[2]+"@"+lins[3]+".A"
            pla_str = plas[0]+":"+plas[2]+"@"+plas[3]+".A"
            
            coords_list.append((pdb_code,cen_str,lin_str,pla_str,atoms))
            print("Added",pdb_code,cen_str,lin_str,pla_str,atoms)

        print(coords_list)

        interpolation = "linear"
        file_output2d = tag + "/sg_2d_"+tag+"_"+pdb_set+"_"
        file_output3d = tag + "/sg_3d_"+tag+"_"+pdb_set+"_"
        samples,width,depth_samples = 50,10,10
        all_vals = None
        all_vals2d = None
        if "2d" in these_stages or "3d" in these_stages:
            count = 1
            last_pdb = ""
            for pdb_code,cen,lin,pla,atoms in coords_list:
                print(count, "/", len(coords_list),pdb_code,cen,lin,pla)
                #try:
                if pdb_code != last_pdb:
                    mman.MapsManager().clear()
                    last_pdb = pdb_code
                ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
                if not ml.success():
                    print(pdb_code, "has not loaded ccp4 succesfully")
                else:
                    mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation,as_sd=2)
                    cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(cen))
                    ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(lin))
                    pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(pla))
                    
                    # 2d plot (s)
                    filename = RESDIR + file_output2d + str(count) + "_" + pdb_code + ".html"                
                    mplot = mph.MapPlotHelp(filename)
                    vals2d = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0)                                                                        
                    if "each2d" in these_stages:
                        if "naybs" in these_stages:                            
                            mfunc = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj, "linear") #the default method is linear
                            naybs = mfunc.get_slice_neighbours(cc,ll,pp,width,samples,[0,1],log_level=1)                        
                            mplot.make_plot_slice_2d(vals2d,min_percent=0.9,max_percent=0.9,samples=samples,width=width,points=[cc,ll,pp],naybs=naybs,title=pdb_code+":"+atoms)
                        else:
                            mplot.make_plot_slice_2d(vals2d,min_percent=0.9,max_percent=0.9,samples=samples,width=width,points=[cc,ll,pp],title=pdb_code+":"+atoms)
                    
                    # 3d plot
                    filename = RESDIR + file_output3d + str(count) + "_" + pdb_code + ".html"                
                    mplot = mph.MapPlotHelp(filename)                
                    vals,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples)                
                    if "each3d" in these_stages:                    
                        mplot.make_plot_slice_3d(vals,min_percent=0.9, max_percent=0.9,title=pdb_code+":"+atoms)
                    
                    count += 1

                    if all_vals == None:
                        all_vals = vals
                        all_vals2d = vals2d
                    else:
                        i,j,k, = all_vals.shape()
                        for a in range(i):
                            for b in range(j):
                                for c in range(k):
                                    new_val = vals.get(a,b,k=c) + all_vals.get(a,b,k=c)
                                    all_vals.add(a,b,k=c,data=new_val)
                                            
                        for a in range(len(all_vals2d)):
                            for b in range(len(all_vals2d)):
                                new_val = all_vals2d[a][b]+vals2d[a][b]
                                all_vals2d[a][b] = new_val
                #except Exception as e:
                #    print(pdb_code,str(e))

            if "all3d" in these_stages and all_vals != None:
                filename = RESDIR + file_output3d + "0_all.html"                
                mplot = mph.MapPlotHelp(filename)                            
                mplot.make_plot_slice_3d(all_vals,min_percent=0.9,max_percent=0.9)

            if "all2d" in these_stages and all_vals2d != None:
                filename = RESDIR + file_output2d + "0_all.html"                
                mplot = mph.MapPlotHelp(filename)                            
                mplot.make_plot_slice_2d(all_vals2d,min_percent=0.9,max_percent=0.9)
            print("-----------------------------------------------")

                

            
            



