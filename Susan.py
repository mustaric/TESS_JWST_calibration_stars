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
        print("No two minute data.")
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
        goodlc = lc.remove_nans().remove_outliers().flatten(window_length = 
                                                            int(0.3*len(lc.time))).normalize()
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