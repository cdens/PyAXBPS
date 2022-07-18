# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    This file is part of the Airborne eXpendable Buoy Processing System (AXBPS)
#
#    AXBPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AXBPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AXBPS.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================


###################################################################################
#                                    IMPORTS                                      #
###################################################################################


import logging

import numpy as np
from scipy import signal, interpolate, optimize, stats
from scipy.io import wavfile #for wav file reading

    
        

###################################################################################
#                         lagging box smooth filter                               #
###################################################################################

def boxsmooth_lag(data,window,startind): #filter is lagging for realtime processing
    outdata = data.copy()
    
    for i in range(startind,len(data)):
        if i < window:
            outdata[i] = np.nanmean(data[0:i+1])
        else:
            outdata[i] = np.nanmean(data[i-window:i+1])
            
    return outdata
    
    
    


###################################################################################
#                         AXCTD PCM DATA FSK DEMODULATION                         #
###################################################################################
    

def demodulate_axctd(pcm, fs, edge_buffer, sos, bitrate, f1, f2, trig1, trig2, Npcm, bit_inset, phase_error, high_bit_scale):
    
    #basic configuration- default options
    # f1 = 400 # bit 1 (mark) = 400 Hz
    # f2 = 800 # bit 0 (space) = 800 Hz
    # bitrate = 800 #symbol rate = 800 Hz
    # bit_inset = 1 #number of points after/before zero crossing where bit identification starts/ends
    # phase_error = 25
    # sos = signal.butter(6, 1200, btype='lowpass', fs=fs, output='sos') #low pass
    # N = int(np.round(fs/bitrate*(1 - phase_error/100))) #first crossing following previous that could be a full bit
    # Npcm = N - 2*bit_inset
    # trig1 = 2*np.pi*np.arange(0,Nn)/fs*f1 #trig term for power calculation
    # trig2 = 2*np.pi*np.arange(0,Nn)/fs*f2
    
    # apply filter to extract FSK data in desired frequency range (<1200 Hz) only
    pcmlow = signal.sosfilt(sos, pcm)
    
    #finding all zero crossings
    pcmsign = np.sign(pcmlow)
    pcmsign[pcmsign == 0] = 1
    zerocrossings = np.where(pcmsign[:-1] != pcmsign[1:])[0]
    
    #ignoring all zero crossings before our starting point
    zerocrossings = np.delete(zerocrossings, np.where(zerocrossings < edge_buffer)) 
    
    #identify zero crossings to separate bits within phase error percentage 
    bit_edges = [zerocrossings[0]]
    prev = zerocrossings[0]
    
    #identifying bit edges
    c = 0 #next available index that could be a bit edge
    while c < len(zerocrossings)-5:
        next_options = zerocrossings[c+1:c+5]
        c += 1 + np.argmin(abs(next_options - (zerocrossings[c] + fs/bitrate)))
        bit_edges.append(zerocrossings[c])
        
    #calculate power at FSK frequencies for each "bit"
    s1 = [] #stores power levels 
    s2 = []
    
    for e in bit_edges[:-1]:
        cdata = pcmlow[e+bit_inset:e+bit_inset+Npcm]
        s1.append(np.abs(np.sum(cdata*np.cos(trig1) + 1j*cdata*np.sin(trig1))))
        s2.append(np.abs(np.sum(cdata*np.cos(trig2) + 1j*cdata*np.sin(trig2)))*high_bit_scale)
        
    next_ind = bit_edges[-1] - 1 #first index to start demodulation on next round
        
    #determine each bit and associated confidence (power of identified freq / power of alternate freq)
    bits = []
    conf = []
    for (p1,p2) in zip(s1,s2):
        conf.append(p2/p1)
        if p1 >= p2:
            bits.append(1)
        else:
            bits.append(0)
            
    return bits, conf, bit_edges, next_ind
            
            
    
###################################################################################
#                           OPTIMIZE DEMOD SCALE FACTOR                           #
###################################################################################

def adjust_scale_factor(confs,scale_factor):
    
    Npts = len(confs)
    confs = np.asarray(confs)
    
    #calculating histogram of confidences
    bin_edges = np.arange(0.0,3,0.01)
    distribution,bin_edges = np.histogram(confs, bins=bin_edges)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges)/2
    
    #calculate cumulative percentage of confidence occurences
    cumpct = 100*np.cumsum(distribution)/Npts
    
    #slope of cumulative percentage function
    cumpctslope = np.array((cumpct[1]-cumpct[0])/(bin_centers[1] - bin_centers[0]))
    cumpctslope = np.append(cumpctslope, (cumpct[2:]-cumpct[:-2])/(bin_centers[2:] - bin_centers[:-2]))
    cumpctslope = np.append(cumpctslope, (cumpct[-1]-cumpct[-2])/(bin_centers[-1] - bin_centers[-2]))
    
    #limiting threshold to cumulative percentages between 30% and 65%
    in_range = [True if cp >= 30 and cp <= 65 else False for cp in cumpct]
    bin_centers = bin_centers[in_range]
    cumpct = cumpct[in_range]
    cumpctslope = cumpctslope[in_range]
    
    #identifying minimum slope and first/last occurences to determine where to cut off bits
    min_slope = np.min(cumpctslope)
    isminslope = np.where(cumpctslope == min_slope)[0]
    new_threshold = np.nanmean([bin_centers[isminslope[0]], bin_centers[isminslope[-1]]])
    
    
    #calcuating new scale factor
    scale_factor /= new_threshold
    
    return scale_factor










    
    