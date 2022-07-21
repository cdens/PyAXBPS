#! /usr/bin/env python3

#   Purpose: process data from AXCTD audio files
#   Usage: via command line: -i flag specifies input audio (WAV) file and -o
#           specifies output text file containing AXCTD profile
#       e.g. "python3 processAXCTD -i inputaudiofile.wav -o outputASCIIfile.txt"
#   See README for additional usage instructions
#
#   This file is a part of AXCTDprocessor
#
#    AXCTDprocessor in this file is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AXCTDprocessor is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AXCTDprocessor.  If not, see <https://www.gnu.org/licenses/>.
#


###################################################################################
#                                    IMPORTS                                      #
###################################################################################

import argparse
import logging
import sys
import os

import numpy as np
from PyAXBPS import AXCTDprocessor

    
    
###################################################################################
#                           ARGUMENT PARSING + MAIN                               #
###################################################################################

#function to handle input arguments

def main():
    
    parser = argparse.ArgumentParser(description='Demodulate an audio file to text')
    parser.add_argument('-i', '--input', default='ERROR_NO_FILE_SPECIFIED', help='Input WAV filename')
    parser.add_argument('-o', '--output', default='output.txt', help='Output filename')
    
    parser.add_argument('-s', '--starttime', default='0', help='AXCTD start time in WAV file') #13:43
    parser.add_argument('-e', '--endtime',  default='-1', help='AXCTD end time in WAV file') #20:00
    
    parser.add_argument('-a','--autodetect-start',  default='30', help='Point at which autodetect algorithm starts scanning for profile transmission start')
    parser.add_argument('-b','--autodetect-end',  default='-1', help='Point at which autodetect algorithm stops scanning for profile transmission start')
    
    parser.add_argument('-p','--sig-threshold-400',  default='2', help='Threshold for normalized 400 Hz signal level to detect profile transmission')
    parser.add_argument('-t','--sig-threshold-7500',  default='1.5', help='Threshold for normalized 7500 Hz signal level to detect profile transmission')
    parser.add_argument('-d','--dead-freq',  default='3000', help='"Dead" (quiet) frequency used to calculate normalized signal levels (Hz)')
    parser.add_argument('-l','--pointsperloop',  default='100000', help='Number of PCM audio data points processed per iteration')
    
    parser.add_argument('-m','--mark-freq',  default='400', help='Mark (bit 1) frequency (Hz)')
    parser.add_argument('-n','--space-freq',  default='800', help='Space (bit 0) frequency (Hz)')
    parser.add_argument('-u','--use-bandpass',  action='store_true', help='Apply this flag to use a bandpass filter (100 Hz to 1200 Hz) rather than a 1200 Hz lowpass filter before demodulation')
    
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
    
    #reading settings for profile transmission detection
    triggerrange = [parse_times(args.autodetect_start), parse_times(args.autodetect_end)]
    if triggerrange[0] < 0:
        triggerrange[0] == 0
    if triggerrange[1] <= 0:
        triggerrange[1] = -1
    
    settings = {'triggerrange': triggerrange,
                'minR400': float(args.sig_threshold_400),
                'mindR7500': float(args.sig_threshold_7500),
                'deadfreq': float(args.dead_freq),
                'pointsperloop': int(args.pointsperloop),
                'mark_space_freqs' : [float(args.mark_freq), float(args.space_freq)],
                'use_bandpass' : args.use_bandpass}
    
    return processAXCTD(args.input, args.output, timerange, settings)
    

    
    
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

        
        
        
def processAXCTD(wavfile,outfile,timerange,settings):
    
    minR400 = settings['minR400']
    mindR7500 = settings['mindR7500']
    deadfreq = settings['deadfreq']
    pointsperloop = settings['pointsperloop'] #about 2 sec of data
    triggerrange = settings['triggerrange']
    mark_space_freqs = settings['mark_space_freqs']
    use_bandpass = settings['use_bandpass']
        
    #running AXCTD Processor instance, timing
    print("Processing profile")
    ap = AXCTDprocessor.AXCTD_Processor(wavfile, timerange=timerange, user_settings=settings)
    ap.run()
    print("Profile processing complete- writing output files")
    
    
    #saving output file with profile information
    with open(outfile, 'w') as f:
        
        f.write(f"AXCTD profile for {wavfile}\n")
        
        #basic file info (fs, file length, trigger times)
        fs = ap.f_s
        f.write(f'Sampling frequency (fs): {fs} Hz\n')
        f.write(f'Audio file length: {ap.numpoints/fs} sec\n')
        f.write(f'400 Hz pulse start: {ap.firstpulse400/fs} sec\n')
        f.write(f'7500 Hz tone start: {ap.profstartind/fs} sec\n')
        
        #header data
        f.write("\nAXCTD header information:\n")
        for desc,ckey in zip(['Probe Code', 'Maximum Depth (m)','Probe Serial'],['probe_code','max_depth','serial_no']):
            f.write(f"{desc}: {ap.metadata[ckey]}\n")
        f.write("Conversion equations:\n")
        for coeff,desc,symb in zip(['z','t','c'], ['Depth','Temperature','Conductivity'],['t','T','C']):
            if sum(ap.metadata[coeff + 'coeff_valid']) == 4:
                cfield = coeff + 'coeff'
                defaultstatus = ''
            else:
                cfield = coeff + 'coeff_default'
                defaultstatus = '(default)'
            cureqn = ' + '.join([f'{val}*{symb}^{i}' for i,val in enumerate(ap.metadata[cfield])])
            f.write(f'{desc}: {cureqn} {defaultstatus}\n')
        
        #processor settings
        f.write('\nProcessor Settings:\n')
        f.write(f'Time Range: {timerange[0]} sec to {timerange[1] if timerange[1] >= 0 else "N/A"} sec\n')
        f.write(f'Min. 400 Hz power ratio: {minR400}\n')
        f.write(f'Min. 7500 Hz power ratio: {mindR7500}\n')
        f.write(f'Dead frequency: {deadfreq}\n')
        f.write(f'Points per loop: {pointsperloop}\n')
        f.write(f'Trigger range: {triggerrange[0]} sec to {triggerrange[1] if triggerrange[1] >= 0 else "N/A"} sec\n')
        
        #profile data
        f.write('\nAXCTD Profile:\n')
        f.write('Time (s), Hex Frame, Depth (m), Temperature (C), Conductivity (mS/cm), Salinity (PSU)\n')
        for (t,hf,z,T,C,S) in zip(ap.time, ap.hexframes, ap.depth, ap.temperature, ap.conductivity, ap.salinity):
            f.write(f"{t:8.2f},  {hf},{z:10.2f},{T:16.2f},{C:21.2f},{S:15.2f}\n")
        
        

#MAIN
if __name__ == "__main__":
    main()



