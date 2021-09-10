#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 08:20:52 2021

@author: smullally
"""

#Code to look at spectrophotometric standards.

import numpy as np
import lightkurve as lk
import matplotlib.pyplot as plt


def examine_tic(ticid):
    
    lcs = lk.search_lightcurve('TIC {}'.format(ticid), author='SPOC', 
                               exptime=120)
    
    
    if len(lcs) == 0:
        print("No two minute data data.")
        return

    stats = np.zeros(len(lcs))
    
    for i,lcf in enumerate(lcs):
        print()
        print(lcf)
        
        myfig = plt.figure(figsize=(10,13))
        ax1 = myfig.add_subplot(3,1,1)
        
        lc = lcf.download()
        lc.scatter(ax=ax1, label="raw")
        plt.legend()
        plt.title('TIC {} Sector {}'.format(ticid, lc.sector))
        
        ax2 = myfig.add_subplot(3,1,2)
        wl = int(0.3*len(lc.time.value))
        if (wl % 2) == 0:
           wl = wl+1 
        goodlc = lc.remove_nans().remove_outliers().flatten(window_length = 
                                                            wl).normalize()
        goodlc.flux = goodlc.flux*1000000
        goodlc.plot(ax=ax2, label="normalize")
        plt.legend()
        
        per = goodlc.to_periodogram(method='lombscargle', 
                                    normalization='amplitude')
        
        N=len(goodlc.flux)
    
        sigma_rms = np.nanstd(goodlc.flux)
        # mean noise level in amplitude spectrum (See A2 of Kledson and Bedding)
        sigma_amp = np.sqrt(np.pi/N)*sigma_rms
        
        stats[i] = 5 * sigma_amp
        
        ax3 = myfig.add_subplot(3,1,3)
        per.plot(ax=ax3)
        plt.ylabel('Amplitude')
        plt.hlines(5*sigma_amp, np.min(per.frequency.value), np.max(per.frequency.value),
                   color='red', linestyles='dashed', label=str(stats[i]))
        plt.legend()
        
        print('5sigma_amp:{}'.format(stats[i]))
        print(per.period_at_max_power.value)
    
    return stats


def create_detrended_lcs(ticid, 
                         ddir="/Users/smullally/Science/tess_monitor_standards/detrended_standards/"):
    """
    Create several versions of your light curve using 
    different parameters.  Write them to disk.
    """
    
    print("hello")
    lcs = lk.search_lightcurve('TIC {}'.format(ticid), author='TESS-SPOC') 
    #                           exptime=[1800,600])
    
    print(lcs)
    
    if len(lcs) == 0:
        print("No two minute data.")
        return

    
    for i,lcf in enumerate(lcs):
    
        print(i)    
        mylc = lcf.download().remove_outliers(sigma = 6)
        rootname = 'tic%014u_s%03u_' % (ticid, mylc.sector)
        
        tL = len(mylc.time)
        
        name1="flat1.fits"
        wl = int(0.5*tL)
        if (wl % 2) == 0:
           wl = wl+1 
        lc1 = mylc.flatten(window_length=wl, niters=2).normalize()
        lc1.to_fits(ddir+rootname+name1, overwrite=True)
        lc1.plot()
        plt.title("Long Window Flatten")
        
        name2 = "flat2.fits"
        wl = int(0.2*tL)
        if (wl % 2) == 0 :
            wl = wl+1
        lc2 = mylc.flatten(window_length=wl, niters=2).normalize()
        lc2.to_fits(ddir+rootname+name2, overwrite=True)
        lc2.plot()
        plt.title("Short Window Flatten")
        
        name3 = "norm1.fits"
        lc3 = mylc.normalize()
        lc3.to_fits(ddir+rootname+name3, overwrite=True)
        lc3.plot()
        plt.title("Normalize only")
    
    return lcs


def create_lcft_plot(filename, periods=[0.00278, 12], times=None,
                     label=None):
    "periods is a min max periods for the periodogram in days."
    sz=17
    
    lc = lk.read(filename)
    pgram = lc.to_periodogram(method="lombscargle", normalization="amplitude",
                           minimum_frequency=1/periods[1], maximum_frequency=1/periods[0])

    flatlc = lc.flatten(window_length=1441, polyorder=4)
    #flatlc = lc    
    #hertz = freq*0.000011574074074
    #per_hour = freq*0.041666666667

    
    plt.figure(figsize=(19,5.54))
    plt.subplot(121)
    plt.plot(lc.time.value, 100*(lc.flux.value-1.0), 'k.', label=label)
    if times != None:
        plt.xlim(times)
    plt.xticks(fontsize=sz)
    plt.yticks(fontsize=sz)
    plt.ylabel('Amplitude [percent]', fontsize=sz)
    plt.xlabel("Time - 2457000 [BJD]", fontsize=sz)
    plt.legend()
    
    plt.subplot(122)
    freq=pgram.frequency  #In per days
    microhertz=11.574074074*freq
    amp=pgram.power * 100
    
    if label == None:
        ticid = int(filename[-28:-16])
        sector = int(filename[-14:-11])
        label = 'TIC {}'.format(ticid)+'\nSector {}'.format(sector)
    else:
        label = label
        
    plt.plot(microhertz, amp, color='k', label=label)
    
    plt.legend()    
    plt.xticks(fontsize=sz)
    plt.yticks(fontsize=sz)

    plt.xlabel(r'Frequency [$\mu$Hz]', fontsize=sz)
    plt.ylabel('Amplitude [percent]', fontsize=sz)
    
    axes1 = plt.gca()
    axes2 = axes1.twiny()
    
    xloc = np.linspace(1e-4, 1e6/(24*3600*periods[0]), 7)
    period_values = 277.777778/np.array(xloc)
    xlabels= np.array(list( map( lambda x: "%.1f" % x, period_values )))
    #xlabels = np.round((1/np.array(xloc))*1e6/(3600))
    #print(xloc,xlabels)
    axes2.set_xticks(xloc[1:])
    axes2.set_xticklabels(xlabels[1:], fontsize=sz-2)
    plt.xlabel('Period [hr]', fontsize=sz-2)
    
    plt.xlim(0, 1e6/(24*3600*periods[0]))
    

    N=len(lc.flux)
    sigma_rms = np.nanstd(flatlc.flux.value)
    
    # mean noise level in amplitude spectrum (See A2 of Kledson and Bedding)
    sigma_amp = np.sqrt(np.pi/N)*sigma_rms
    #plt.hlines(5*100*sigma_amp, 0, 11.5741*np.max(pgram.frequency.value),
    #               color='red', linestyles='dashed', label=r"5$\sigma$")
    #plt.annotate("5"+r"$\sigma$", (xloc[-2], 5.2*100*sigma_amp), color='red', fontsize=sz-2)

def calc_variable_stats(filename, ticid = 111, periods=[.013,13]):
    
    lc = lk.read(filename)
    
    pgm = lc.to_periodogram(method="lombscargle", normalization="amplitude",
                           minimum_frequency=1/periods[1], 
                           maximum_frequency=1/periods[0])
    N = len(lc.time.value)
    ave_power = np.nanmean(pgm.power.value ** 2)
    ft_noise_level = (6 * ave_power )**(0.5) #Used for peak 2 peak estimate
    
    #pgm.plot()
    max_period = pgm.period_at_max_power.value
    max_amplitude = pgm.max_power.value * 100
    flux = lc.flux.value
    binflux = lc.bin(time_bin_size=(.0042), aggregate_func=np.nanmedian) #6 minutes
    flux = binflux.flux.value
    #binflux.plot()
    #peakpeak = 100*(np.nanmax(flux) - np.nanmin(flux))
    peakpeak = 100*(np.nanpercentile(flux,97.725) - np.nanpercentile(flux,2.275))
    peakpeak2 = 100 * peak2peak(lc, ft_noise_level, periods=periods)
    
    return ticid, max_period, max_amplitude, peakpeak, peakpeak2

def peak2peak(lc, min_amp, periods=[.013,13], N=4):
    """

    Parameters
    ----------
    lc : lightkurve module
        must contain time and flux.

    Returns
    -------
    peak2peak : float

    Add up amplitudes of N highest peaks.
    """
    
    pgm = lc.to_periodogram(method="lombscargle", normalization="amplitude",
                           minimum_frequency=1/periods[1], 
                           maximum_frequency=1/periods[0])

    #pgm.plot()
    #print(min_amp)

    max_amp = np.zeros(N)
    resolution = 3/(lc.time.value[-1] - lc.time.value[0])
    #print(resolution)
    
    for i in range(0,N):
        
        max_loc = np.argmax(pgm.power.value)
        max_freq = pgm.frequency.value[max_loc]
        new_amp = pgm.power.value[max_loc]
        #print(new_amp, min_amp)
        if new_amp > min_amp:
            max_amp[i] = pgm.power.value[max_loc]
        
        #print(max_loc, max_freq, max_amp[i])
        
        freq = pgm.frequency.value
        freq_range = (freq > (max_freq - resolution)) & (freq < (max_freq + resolution))
        pgm.power.value[freq_range] = 0

        #pgm.plot()
    
    pk2pk = 2 * np.nansum(max_amp)
    
    return pk2pk

def nonvar_stats(filename, period = 1, long_bin_size = .02083):
    """
    

    Parameters
    ----------
    filename : string
        DESCRIPTION.
    period_break : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    two stats
    short_period_stat : float
        3 sigma times the average power with periods shorter than period
        
    long_period_stat : float
        3 times the rms of the light curve.

    """
    
    
    lc = lk.read(filename)
    
    pgm = lc.to_periodogram(method="lombscargle", normalization="amplitude",
                           minimum_frequency=1/period)
    
    ave_power = np.nanmean(pgm.power.value ** 2)
    short_period_stat = np.nanmax(pgm.power.value)
    #short_period_stat = np.sqrt(4 * ave_power)
    
    
    pgm = lc.to_periodogram(method="lombscargle", normalization="amplitude",
                           maximum_frequency=1/period)
    ave_power = np.nanmean(pgm.power.value **2)
    #long_period_stat = np.sqrt(4 * ave_power)
    long_period_stat = np.nanmax(pgm.power.value)
    
    binflux = lc.bin(time_bin_size=long_bin_size, aggregate_func=np.nanmedian)
    #Contains 3 sigma of the points
    sig3_upper_limit = np.nanpercentile(binflux.flux.value,99.865) - np.nanpercentile(binflux.flux.value, .13495)
    #upper_limit = 3 * np.nanstd(binflux.flux.value)
    
    return short_period_stat, long_period_stat, sig3_upper_limit
    
    