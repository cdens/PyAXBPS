#! /usr/bin/env/python3

# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    AXBPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AXBPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AXBPS.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

import numpy as np
from scipy.io import wavfile #for wav file reading
from scipy.signal import tukey #taper generation

import wave #WAV file writing

import time as timemodule
import datetime as dt

from traceback import print_exc as trace_error

from shutil import copy as shcopy
from sys import getsizeof





# #initializes dict with default settings for AXCTD processor
# #any specified settings will be overwritten after this function is called in __init__
# def init_AXCP_settings(self, settings):

#     self.settings = {}
#     self.settings["refreshrate"] = 1 #size of raw audio chunks to process, in seconds
#     self.settings["revcoil"]     = False #coil on AXCP reversed- rotates currents by 180 degrees
#     self.settings["quality"]     = 1    #profile processing quality- 1=high/slow, 2=moderate speed/quality, 3=low/fast
#     self.settings["spindown_detect_rt"] = True #realtime detection of probe spindown to avoid processing unnecessary data
    
#     #1-use zero crossings and temperature baseline freq calculated with velocity
#     #2-use FFT (window size set by settings['tempfftwindow']) once per refresh
#     self.settings["temp_mode"] = 2 
#     self.settings["tempfftwindow"] = 1.0 #FFT window for temperature in seconds (used if temp_mode > 1) 
            
#     for csetting in settings:
#         self.settings[csetting] = settings[csetting] #overwrite defaults for user specified settings 
      
#     self.revcoil = bool(self.settings["revcoil"])
#     self.spindown_detect_rt = bool(self.settings["spindown_detect_rt"])
    
#     self.quality = self.settings["quality"]
#     if self.quality <= 0:
#         self.quality = 1
#     elif self.quality >= 4:
#         self.quality = 3
    
#     self.temp_mode = self.settings["temp_mode"]
#     if self.temp_mode < 1:
#         self.temp_mode = 1
#     elif self.temp_mode > 2:
#         self.temp_mode = 2    

    
#     self.refreshrate = self.settings["refreshrate"]
#     self.tempfftwindowsec = self.settings["tempfftwindow"]
    
#     if self.tempfftwindowsec > self.refreshrate: #temperature FFT window length must be less than refresh rate
#     self.tempfftwindowsec = self.refreshrate


    
# =============================================================================
#  Audio Processor class
# =============================================================================


class AXBT_Processor:

    #initializing current thread (saving variables, reading audio data or contacting/configuring receiver)
    def __init__(self, audiofile, audiochannel=-1, timerange=[0,-1], settings={}):
        
        self.audiochannel = audiochannel
        
        self.keepgoing = True  # signal connections
        self.waittoterminate = False #whether to pause on termination of run loop for kill process to complete
        
        self.timerange = timerange
        
        self.temperature = []
        self.depth = []
        self.time = []
        self.Rp = []
        self.Sp = []
        self.fp = []
        
        self.audiofile = audiofile
        
        self.init_settings(settings)
        
        
        
    def init_settings(self,settings):
        
        self.settings = {'fftwindow':0.3, 'minfftratio':0.5, 'minsiglev':65, 'triggerfftratio':0.88, 'triggersiglev':75, 'tcoeff':[-40.0,0.02778,0.0,0.0], 'zcoeff':[0.0,1.524,0.0,0.0], 'flims':[1300,2800]}
        
        #overwriting user-defined settings
        for csetting in settings:
            self.settings[csetting] = settings[csetting]
        
        #FFT thresholds
        self.fftwindow = self.settings['fftwindow']
        self.minfftratio = self.settings['minfftratio']
        self.minsiglev = self.settings['minsiglev']
        self.triggerfftratio = self.settings['triggerfftratio']
        self.triggersiglev = self.settings['triggersiglev']
        
        #conversion coefficients + parameters
        self.tcoeff = self.settings['tcoeff']
        self.zcoeff = self.settings['zcoeff']
        self.flims = self.settings['flims']
        
        
        
        
        
        
    def init_window(self,N):
        
        # apply taper- alpha=0.25
        self.taper = tukey(N, alpha=0.25)
        
        self.N = N
        
        T = N/self.f_s
        df = 1 / T
        self.f = np.array([df * n if n < N / 2 else df * (n - N) for n in range(N)])#constraining peak frequency options to frequencies in specified band
        self.good_f_ind = np.all((np.greater_equal(self.f, self.flims[0]), np.less_equal(self.f, self.flims[1])), axis=0)
        
        self.good_f = self.f[self.good_f_ind]
                
        
        
    
    #convert time to depth, freq to temp given coefficient lists
    def btconvert(self,input,coefficients):
        output = 0
        for (i,c) in enumerate(coefficients):
            output += c*input**i
        return output
        
        
        
    
    #function to run fft here
    def dofft(self,pcmdata):
        pcmdata = self.taper * pcmdata
    
        # conducting fft, converting to real space
        fftdata = np.abs(np.fft.fft(pcmdata))
        fftdata_inrange = fftdata[self.good_f_ind]
        
        maxind = np.argmax(fftdata_inrange)
        
        #frequency of max signal within band (AXBT-transmitted frequency)
        fp = self.good_f[maxind] 
        
        #maximum signal strength in band
        Sp = 10*np.log10(fftdata_inrange[maxind])
    
        #ratio of maximum signal in band to max signal total (SNR)
        Rp = fftdata_inrange[maxind]/np.max(fftdata) 
        
        # #Signal level normalized by power at quiet frequency
        # pdead = fftdata[self.dead_freq_ind]
        # Sp = np.log10(fftdata_inrange[maxind]/pdead)
            
        return fp, Sp, Rp
        
        
        
        
    def run(self):
        
        #checking file length- wont process files with more frames than max size
        try: #exception if unable to read audio file if it doesn't exist or isn't WAV formatted
            file_info = wave.open(self.audiofile)
        except:
            print("[!] Error- unable to open audio file! Terminating...")
            return
        
        self.f_s, snd = wavfile.read(self.audiofile) #reading file
        
        #if multiple channels, sum them together
        sndshape = np.shape(snd) #array size (tuple)
        ndims = len(sndshape) #number of dimensions
        if ndims == 1: #if one channel, use that
            self.audiostream = snd
            
        elif ndims == 2: #if multiple channels, pick selected channel, otherwise sum
            if self.audiochannel >= 1:
                if sndchape[1] > sndshape [0]:
                    self.audiostream = snd[self.audiochannel-1,:]
                else:
                    self.audiostream = snd[:,self.audiochannel-1]
                    
            else:
                self.audiostream = np.sum(snd,axis=1)
                print("[+] Multiple channels in file and channel not specified- using sum of all channels")
                
        else: #if more than 2D- not a valid file
            print("[!] Error- audio file has more than two dimensions! Terminating...")
            return
                
        #trimming audio as required
        if self.timerange[1] > 0: #trim end of file first to keep this simple
            e = int(round(self.f_s * self.timerange[1]))
            self.audiostream = self.audiostream[:e]
        if self.timerange[0] > 0:
            s = int(round(self.f_s * self.timerange[0]))
            self.audiostream = self.audiostream[s:]
            
        #configuring sample times for the audio file
        self.lensignal = len(self.audiostream)
        self.maxtime = self.lensignal/self.f_s
        self.sampletimes = np.arange(0.1,self.maxtime-0.1,0.1)
        self.ntimes = len(self.sampletimes)
        
        self.istriggered = False
        self.starttime = 0
            
        # setting up thread while loop- terminates when user clicks "STOP" or audio file finishes processing
        i = -1
            
        
        N_cur = -1
        
        #MAIN PROCESSOR LOOP
        while self.keepgoing:
            i += 1 

            #kill test/audio threads once time exceeds the max time of the audio file
            #NOTE: need to do this on the cycle before hitting the max time when processing from audio because
            #       the WAV file processes faster than the thread can kill itself
            if  i >= len(self.sampletimes)-1:
                self.keepgoing = False
                print(f"[+] Processing status: Completed!")
                return
                
            #getting current time to sample from audio file
            ctime = self.sampletimes[i]

            #getting current data to sample from audio file- using indices like this is much more efficient than calculating times and using logical arrays
            ctrind = int(np.round(ctime*self.f_s))
            pmind = int(np.min([np.round(self.f_s*self.fftwindow/2),ctrind,self.lensignal-ctrind-1])) #uses minimum value so no overflow
            currentdata = self.audiostream[ctrind-pmind:ctrind+pmind]
                
            #setting up frequency array as reqd
            N = len(currentdata)
            if N_cur != N:
                self.init_window(N)
                N_cur = N

            #conducting FFT or skipping, depending on signal strength
            fp,Sp,Rp = self.dofft(currentdata)        
    
            #rounding before comparisons happen
            ctime = np.round(ctime, 1)
            fp = np.round(fp, 2)
            Sp = np.round(Sp, 2)
            Rp = np.round(Rp, 3)        
            

                
            #logic to determine whether or not profile is triggered
            if not self.istriggered and Sp >= self.triggersiglev and Rp >= self.triggerfftratio:
                self.istriggered = True
                self.firstpointtime = ctime
                    
            #logic to determine whether or not point is valid
            if self.istriggered and Sp >= self.minsiglev and Rp >= self.minfftratio:
                cdepth = self.btconvert(ctime - self.firstpointtime, self.zcoeff)
                ctemp = self.btconvert(fp, self.tcoeff)
                
                #rounding
                ctemp = np.round(ctemp, 2)
                cdepth = np.round(cdepth, 1)
                
            else:
                ctemp = np.nan
                cdepth = np.nan
                
            #saving points to lists
            self.temperature.append(ctemp)
            self.depth.append(cdepth)
            self.time.append(ctime)
            self.Rp.append(Rp)
            self.Sp.append(Sp)
            self.fp.append(fp)
                
            
            print(f"[+] Processing status: {round(100*i/self.ntimes)}%         ", end='\r')
            timemodule.sleep(0.001) #slight pause to free some resources when processing from audio
                
            
        
        

        
        
        