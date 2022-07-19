#!/usr/bin/env python3
# =============================================================================
#     Author: Casey R. Densmore,
#
#    This file is part of the AXBT Processor
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

import os, shutil
import numpy as np
from PyAXBPS import AXBTprocessor as ap
from datetime import datetime
import argparse



def writeedffile(edffile,dropdatetime,lat,lon,data,comments,QC=False):
    
    
    with open(edffile,'w') as f_out:
    
        #writing header, date and time, drop # (bad value)
        f_out.write("// This is an air-launched probe EXPORT DATA FILE  (EDF)\n")
        f_out.write("// File generated with AXBT Processor \n")
        if dropdatetime is not None:
            f_out.write(f"Date of Launch:  {datetime.strftime(dropdatetime,'%m/%d/%y')}\n")
            f_out.write(f"Time of Launch:  {datetime.strftime(dropdatetime,'%H:%M:%S')}\n")
        
        if lat is not None and lon is not None:
            #latitude and longitude in degrees + decimal minutes
            if lat >= 0:
                nsh = 'N'
            else:
                nsh = 'S'
            if lon >= 0:
                ewh = 'E'
            else:
                ewh = 'W'
            lon = abs(lon)
            lat = abs(lat)
            latdeg = int(np.floor(lat))
            londeg = int(np.floor(lon))
            latmin = (lat - latdeg)*60
            lonmin = (lon - londeg)*60
            f_out.write(f"Latitude      :  {latdeg:02d} {latmin:06.3f}{nsh}\n")
            f_out.write(f"Longitude     :  {londeg:03d} {lonmin:06.3f}{ewh}\n")        
        
        if QC:
            qcstr = " "
        else:
            qcstr = """
// This profile has not been quality-controlled.
//
"""
        
        #building drop metadata comments section
        drop_settings_info = f"""//
// Drop Settings Information:
""" + comments + """
""" + "\n".join([f"Field{i+1}  :  {ckey}" for i,ckey in enumerate(data.keys())]) + """
//""" + qcstr + """
""" + "\t".join([ckey for ckey in data.keys()]) + "\n"
        
        f_out.write(drop_settings_info)
        
        fields = list(data.keys())
        
        #determining what string format to use for each field
        field_formats = []
        for cfield in fields:
            cfieldlow = cfield.lower()
            if 'temperature' in cfieldlow or 'conductivity' in cfieldlow or 'salinity' in cfieldlow:
                field_formats.append('7.2f') #__XX.XX
            elif 'current' in cfieldlow:
                field_formats.append('6.2f') #__X.Xx
            elif 'depth' in cfieldlow or 'time' in cfieldlow or 'frequency' in cfieldlow:
                field_formats.append('9.2f') #__XXXX.XX
            else:
                field_formats.append('') #unspecified format
                
        npts = len(data[fields[0]])
        
        #writing data
        for i in range(npts):
            cline = "\t".join([f"{data[cfield][i]:{cformat}}" for cfield,cformat in zip(fields,field_formats)]) #tab-delimited
            f_out.write(cline + "\n")

            




#constant settings for audio reprocessing
fftwindow = 0.3
minfftratio = 0.5
minsiglev = 65.0
triggerfftratio = 0.88
triggersiglev = 75.0
tcoeff = [-40.0,0.02778,0.0,0.0]
zcoeff = [0.0,1.524,0.0,0.0]
flims = [1300,2800]




    
    
    
    
    
###################################################################################
#                           ARGUMENT PARSING + MAIN                               #
###################################################################################

#function to handle input arguments

def main():
    
    parser = argparse.ArgumentParser(description='Demodulate an audio file to text')
    parser.add_argument('-i', '--input', default='ERROR_NO_FILE_SPECIFIED', help='Input WAV filename')
    parser.add_argument('-c', '--audio-channel', default='-1', help='Input channel for audio file (1 for left, 2 for right, -1 to sum across all channels)')
    parser.add_argument('-o', '--output', default='output.edf', help='Output filename')
    
    parser.add_argument('-s', '--starttime', default='0', help='AXBT start time in WAV file') #13:43
    parser.add_argument('-e', '--endtime',  default='-1', help='AXBT end time in WAV file') #20:00
    
    parser.add_argument('-w','--fft-window',  default='0.3', help='Window length for each FFT (sec)')
    parser.add_argument('-m','--min-signal',  default='65', help='Minimum signal level for good data (dB)')
    parser.add_argument('-r','--min-ratio',  default='0.5', help='Minimum signal ratio for good data (0-1)')
    parser.add_argument('-a','--trigger-signal',  default='75', help='Minimum signal level to trigger profile collection (dB)')
    parser.add_argument('-b','--trigger-ratio',  default='0.88', help='Minimum signal ratio to trigger profile collection')
    
    parser.add_argument('-t','--tcoeff',  default='[-40.0,0.02778,0.0,0.0]', help='Frequency (Hz) to temperature (degC) conversion coefficients: [C0,C1,C2,C3] where T = C0 + C1*f + C2*f^2 + C3*f^3')
    parser.add_argument('-z','--zcoeff',  default='[0.0,1.524,0.0,0.0]', help='Time elapsed (sec) to depth (m) conversion coefficients: [C0,C1,C2,C3] where z = C0 + C1*t + C2*t^2 + C3*t^3')
    parser.add_argument('-f','--freq-lims',  default='[1300,2800]', help='Minimum and maximum frequencies for good AXBT data: [min_f,max_f]')
    
    
    args = parser.parse_args()    
    
    #checking for input WAV file
    if args.input == 'ERROR_NO_FILE_SPECIFIED':
        print("[!] Error- no input WAV file specified! Terminating")
        exit()
    elif not os.path.exists(args.input):
        print("[!] Specified input file does not exist! Terminating")
        exit()
    
    
    #WAV time bounds for processing
    timerange = [parse_times(args.starttime), parse_times(args.endtime)]
    if timerange[0] < 0:
        timerange[0] == 0
    if timerange[1] <= 0:
        timerange[1] = -1
        
    #pulling coefficients as required
    tcoeff = [float(c) for c in args.tcoeff.replace('[','').replace(']','').strip().split(',')]
    zcoeff = [float(c) for c in args.zcoeff.replace('[','').replace(']','').strip().split(',')]
    flims = [float(c) for c in args.freq_lims.replace('[','').replace(']','').strip().split(',')]
    
    #settings for processor
    settings = {'fftwindow': float(args.fft_window),
                'minsiglev': float(args.min_signal),
                'minfftratio': float(args.min_ratio),
                'triggersiglev': float(args.trigger_signal),
                'triggerfftratio' : float(args.trigger_ratio),
                'tcoeff' : tcoeff,
                'zcoeff' : zcoeff,
                'flims'  : flims}
    
    processAXBT(args.input, args.output, timerange, int(args.audio_channel), settings)
    

    
    
def parse_times(time_string):
    try:
        if ":" in time_string: #format is HH:MM:SS 
            t = 0
            for i,val in enumerate(reversed(time_string.split(":"))):
                if i <= 2: #only works up to hours place
                    t += int(val)*60**i
                else:
                    logging.info("[!] Warning- ignoring all end time information past the hours place (HH:MM:SS)")
        else:
            t = int(time_string)
        return t
        
    except ValueError:
        logging.info("[!] Unable to interpret specified start time- defaulting to 00:00")
        return -2

        
    
    
#process audio file, return output
def processAXBT(input_file, output_file, timerange, audiochannel, settings):
    
    #run audio processor
    axbt = ap.AXBT_Processor(input_file, timerange=timerange, audiochannel=audiochannel, fftwindow=settings['fftwindow'], minfftratio=settings['minfftratio'], minsiglev=settings['minsiglev'], triggerfftratio=settings['triggerfftratio'], triggersiglev=settings['triggersiglev'], tcoeff=settings['tcoeff'], zcoeff=settings['zcoeff'], flims=settings['flims'])
    axbt.run()
    
    temperature = axbt.temperature
    depth = axbt.depth
    time = axbt.time
    frequency = axbt.fp
    signal = axbt.Sp
    ratio = axbt.Rp
    
    
    print("Profile processing complete- writing output files")
    
    data = {'Time (s)':time, 'Frequency (Hz)':frequency, 'Sp (dB)':signal, 'Rp (%)':ratio, 'Depth (m)':depth, 'Temperature (degC)':temperature}
    
    comments = f"""Probe Type : AXBT
    Depth Coefficientd : {settings['zcoeff']}
    Temperature Coefficients : {settings['tcoeff']}
    Frequency Bounds (Hz) : {settings['flims']}
    FFT Window (sec): {settings['fftwindow']}
    Trigger signal level (dB)/ratio (%) : {settings['triggersiglev']} / {settings['triggerfftratio']}
    Minimum signal level (dB)/ratio (%) : {settings['minsiglev']} / {settings['minfftratio']}
    """
    
    writeedffile(output_file,None,None,None,data,comments,QC=False)
    

#MAIN
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    