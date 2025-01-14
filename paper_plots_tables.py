#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 11:57:52 2021

@author: smullally
"""

import Susan 
import matplotlib.pyplot as plt
import numpy as np

outdir = "/Users/smullally/Science/tess_monitor_standards/paper/plots/"
#Directory of the data files
ddir = "/Users/smullally/Science/tess_monitor_standards/detrended_standards/good/"
#The list of names and filenames
infile = "/Users/smullally/Science/tess_monitor_standards/paper/plots/inputfilenames.csv"

filenames = np.loadtxt(infile, dtype=str, delimiter=',')

#This needs to be run first and then gets filled in below, one at a time.
stats = np.zeros((len(filenames[:,0]),5))

#%%

#for f in filenames[:,0]:
#    Susan.create_lcft_plot(ddir+f, periods = [.004,12], times=None)
#%%

i = 0
pers = [0.5,12]
times = [2174, 2230]
label = "%s\nTIC %u\nSector 32-33" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%%
i = 1
pers = [0.1,12]
times = [2204, 2214.5]
label = "%s\nTIC %u\nSector 33" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%
i = 2
pers = [0.08,12]
times = None
label = "%s\nTIC %u\nSector 5" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%%
i = 3
pers = [0.2,12]
times =[2102, 2113.5]
label = "%s\nTIC %u\nSector 29" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%
i = 4
pers = [0.015,12]
times = [1815.8, 1828]
label = "%s\nTIC %u\nSector 19" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%
i = 5
pers = [0.014,5]
times = [2406,2409]
label = "%s\nTIC %u\nSector 40" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%
i = 6
pers = [0.05,10]
times = None
label = "%s\nTIC %u\nSector 40" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%
i = 7
pers = [0.01,12]
times = [2389.5,2404.9]
label = "%s\nTIC %u\nSector 40" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 8
pers = [0.4,12]
times = [2390, 2405]
label = "%s\nTIC %u\nSector 40" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 9
pers = [0.4,12]
times = [1751,1763.5]
label = "%s\nTIC %u\nSector 16" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 10
pers = [0.2,12]
times = [1855.8,1869]
label = "%s\nTIC %u\nSector 20" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 11
pers = [0.4,12]
times = None
label = "%s\nTIC %u\nSector 21" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#
i = 12
pers = [0.04,8]
times = None
label = "%s\nTIC %u\nSector 33" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 13
pers = [0.6,14]
times = None
label = "%s\nTIC %u\nSector 1" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret
#%
i = 14
pers = [0.2,12]
times = None
label = "%s\nTIC %u\nSector 32" % (filenames[i,1], int(filenames[i,0][7:17]))
Susan.create_lcft_plot(ddir+filenames[i,0], periods = pers, times=times, label = label)
plt.savefig(outdir + filenames[i,0] + ".png")
ret = Susan.calc_variable_stats(ddir+filenames[i,0], ticid = int(filenames[i,0][7:17]), periods=pers)
stats[i,:] = ret

#%%

#Need run above cells first
ofn = "/Users/smullally/Science/tess_monitor_standards/paper/variable_stats.csv"
form = ("%u", "%5.5f", "%.4f", "%.4f", "%.4f")
np.savetxt(fname= ofn, X = stats, delimiter=",", 
           header = "TIC,period_at_max,max_amplitude,2sigma_pkpk,3pk2pk",
           fmt = form)

#%%
#Get Crowdsap for all stars.
import lightkurve as lk
import pandas as p
target_name_file = "/Users/smullally/Science/tess_monitor_standards/paper/target_names.csv"
targets = p.read_csv(target_name_file,delimiter=',', header="infer")

for i,t in targets.iterrows():
    search = lk.search_lightcurve("TIC %u" % t['ticid'], mission='TESS', 
                              author = ['SPOC','TESS-SPOC'])
    try:
        lc = search[search.author=='SPOC'][0].download()
    except:
        lc = search[0].download()    
    print(t['ticid'],lc.CROWDSAP,lc.TIMEDEL)
#%%
#Create the Variable table with Pandas
import pandas as p
target_name_file = "/Users/smullally/Science/tess_monitor_standards/paper/target_names.csv"
targets = p.read_csv(target_name_file,delimiter=',', header="infer")

var_stats = p.read_csv(ofn)
varstars = targets['var'] == 'y'
vstats = p.DataFrame({'ticid':[],'name':[],'crowdsap':[],'period_max':[],'amp_max':[],'pkpk':[]})

for i,t in targets[varstars].iterrows():
    want = var_stats['# TIC'] == t['ticid']
    pmax = float(var_stats[want]['period_at_max'])
    amax = float(var_stats[want]['max_amplitude'])
    pkpk = float(var_stats[want]['2sigma_pkpk'])
    newstat = p.DataFrame({'ticid':str(int(t['ticid'])),'name':t['name'],
                           'crowdsap':t['crowdsap'], 'pmax':[pmax],
                           'amax':[amax], 'pkpk':[pkpk]})
    vstats = vstats.append(newstat)

tabfile = "/Users/smullally/Science/tess_monitor_standards/paper/var_stats.tab"
vstats.to_latex(buf=tabfile, column_format='llrrrr', index = False,
               columns = ['name','ticid','crowdsap','pmax','amax','pkpk'],
               float_format = "%.4f")  


#%%

infile = "/Users/smullally/Science/tess_monitor_standards/detrended_standards/good/all_data.txt"
outfile = "/Users/smullally/Science/tess_monitor_standards/detrended_standards/good/all_data_stats.txt"
ddir = "/Users/smullally/Science/tess_monitor_standards/detrended_standards/good/"
filenames = np.loadtxt(infile,dtype=str)

stats = np.zeros((len(filenames), 5))

for i,f in enumerate(filenames):
    ticid = int(f[5:17])
    sector = int(f[19:22])
    
    #4 sigma long and short period limits from periodogram.
    #Bin size is used to bin for the largest variation stat.
    short, long, largest = Susan.nonvar_stats(ddir + f, period = 1, 
                                     long_bin_size = .02083)
    
    stats[i,:] = np.array([ticid, sector, 100*short, 100*long, 100*largest])

form = ("%u", "%.4f", "%.4f", "%.4f", "%.4f")
np.savetxt(outfile, stats, delimiter = ",", fmt = form, header = "TIC,sector,short,long,3siglarge")

    
#%%
import pandas as p
#The outfile here is the above csv file calculating stats for all data files.
outfile = "/Users/smullally/Science/tess_monitor_standards/detrended_standards/good/all_data_stats.txt"
novar_stats = p.read_csv(outfile)
#The following files comes from the Target_Stars_names excel sheet 1
target_name_file = "/Users/smullally/Science/tess_monitor_standards/paper/target_names.csv"
targets = p.read_csv(target_name_file,delimiter=',', header="infer")

novar_list = "/Users/smullally/Science/tess_monitor_standards/paper/" #contain TIC and actual name

stats = p.DataFrame({'ticid':[],'name':[],'crowdsap':[],'short':[],'long':[],'large':[]})

notvariable_ones = targets['var'] == 'n'

for i,t in targets[notvariable_ones].iterrows():
    want = novar_stats['# TIC'] == t['ticid']
    min_short = np.min(novar_stats[want]['short'])
    min_long = np.min(novar_stats[want]['long'])
    min_large = np.min(novar_stats[want]['3siglarge'])
    newstat = p.DataFrame({'ticid':str(int(t['ticid'])),'name':t['name'], 
                           'crowdsap':t['crowdsap'], 
                           'short':[min_short],'long':[min_long], 
                           'large':[min_large]})
    #print(newstat)
    stats = stats.append(newstat) 

tabfile = "/Users/smullally/Science/tess_monitor_standards/paper/novar_stats.tab"
stats.to_latex(buf=tabfile, column_format='lllrrr', index = False,
               columns = ['name','ticid','crowdsap','short','long','large'],
               float_format = "%.4f")


#%%
#Target Table
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
target_name_file = "/Users/smullally/Science/tess_monitor_standards/paper/target_names.csv"
targets = p.read_csv(target_name_file,delimiter=',', header="infer")

ticids = np.array(targets['ticid'])
tic_data = Catalogs.query_criteria(catalog="Tic",ID = ticids)
tic_data.add_column('skycoord')

#for i,tic in enumerate(tic_data['ID']):
sc = SkyCoord(tic_data['ra'], tic_data['dec'],unit='deg',frame='icrs')
skycoordstr = sc.to_string('hmsdms', precision=0)
tic_data['skycoord'] = skycoordstr
tic_p = tic_data.to_pandas()

#tic_data[['ID','Tmag','Teff','logg','skycoord']]

targ_tab = p.DataFrame({'ticid':[],'name':[],'coord':[],'Tmag':[],'sptype':[],'cadence':[]})
for i,t in targets.iterrows():
    want = tic_p['ID'].astype(int) == int(t['ticid'])
    tmag = "%.1f" % float(tic_p[want]['Tmag'])
    scstr = tic_p[want]['skycoord'].item()
    row = p.DataFrame({'ticid':[str(int(t['ticid']))],'name':[t['name']],
                       'coord':[scstr],'Tmag':[tmag],
                       'sptype':[t['sptype']],'cadence':[t['dtype']]})
    targ_tab = targ_tab.append(row)

targ_file = "/Users/smullally/Science/tess_monitor_standards/paper/target_table.tab"              
targ_tab.to_latex(buf = targ_file, column_format = "llllll", index=False, 
                  columns = ['name','ticid','sptype','Tmag','cadence','coord'])    