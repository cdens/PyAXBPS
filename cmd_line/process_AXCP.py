#!/usr/bin/env python3

# =============================================================================
#     Author: Casey R. Densmore,
#
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
from PyAXBPS import AXCPprocessor as ap
from datetime import datetime, date
import argparse



def writeedffile(edffile,dropdatetime,lat,lon,data,comments,QC=False):
    
    
    with open(edffile,'w') as f_out:
    
        #writing header, date and time, drop # (bad value)
        f_out.write("// This is an air-launched probe EXPORT DATA FILE  (EDF)\n")
        f_out.write("// File generated with AXCP Processor \n")
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
    parser.add_argument('-o', '--output', default='output.edf', help='Output filename')
    
    parser.add_argument('-s', '--starttime', default='0', help='AXCP start time in WAV file') 
    parser.add_argument('-e', '--endtime',  default='-1', help='AXCP end time in WAV file') 
    
    parser.add_argument('-x', '--longitude', default='20', help='Drop longitude (E > 0)') 
    parser.add_argument('-y', '--latitude',  default='-80', help='Drop latitude (N > 0)') 
    parser.add_argument('-d', '--date',  default='-1', help='Drop date (YYYYMMDD)') 
    
    parser.add_argument('-q','--quality',  default='1', help='1 = high quality/slow speed, 2 = moderate quality/speed, 3 = low quality/fast')
    parser.add_argument('-r','--revcoil',  default='0', help='Whether or not coil direction is reversed (1=yes,0=no)')
    parser.add_argument('-p','--refreshrate',  default='0.5', help='Time (seconds) to process per iteration')
    
    parser.add_argument('-u','--spindown_detect',  default='1', help='Whether or not to try to detect probe spindown realtime and automatically stop processing (1=yes,0=no)')
    parser.add_argument('-m','--temperature_mode',  default='2', help='Whether to use zero crossings (1) or FFT (2) to calculate temperature at each depth')
    parser.add_argument('-w','--temp_fft_window',  default='1', help='Temperature FFT window (sec)- only applicable for temperature_mode=2, must be less than or equal to refreshrate.')
    
    #1-use zero crossings and temperature baseline freq calculated with velocity
    #2-use FFT (window size set by settings['tempfftwindow']) once per refresh
    self.settings["temp_mode"] = 2 
    self.settings["tempfftwindow"] = 1.0
    
    
    args = parser.parse_args()    
    
    #checking for input WAV file
    if args.input == 'ERROR_NO_FILE_SPECIFIED' and not os.path.exists(args.input):
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
        
    #pulling position and time
    droplat = float(args.latitude)
    droplon = float(args.longitude)
    
    if float(args.date) == -1:
        dropdate = date.today()
    else:
        dropdate = datetime.strptime(args.date,'%Y%m%d').date()
    
    #settings for processor
    settings = {'quality': int(args.quality),
                'revcoil': int(args.revcoil),
                'refreshrate': float(args.refreshrate),
                'spindown_detect_rt': int(args.spindown_detect),
                'temp_mode': int(args.temperature_mode),
                'tempfftwindow': float(args.temp_fft_window)}
    
    processAXCP(args.input, args.output, timerange, droplat, droplon, dropdate, settings)
    

    
    
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
def processAXCP(input_file, output_file, timerange, droplat, droplon, dropdate, settings):
    
    #run audio processor
    AXCP = ap.AXCP_Processor(input_file, timerange=timerange, lat=droplat, lon=droplon, dropdate=dropdate, settings=settings)
    AXCP.run()
    
    temperature = AXCP.TEMP
    depth = AXCP.DEPTH
    time = AXCP.TIME
    u = AXCP.U_MAG
    v = AXCP.V_MAG
    frot = AXCP.ROTF
    ut = AXCP.U_TRUE
    vt = AXCP.V_TRUE
    
    
    print("Profile processing complete- writing output files")
    
    data = {'Time (s)':time, 'Rotation Rate (Hz)':frot, 'Depth (m)':depth, 'Temperature (degC)':temperature, 'U Current (mag; m/s)':u, 'V Current (mag;m/s)':v, 'U Current (true; m/s)':ut, 'V Current (true; m/s)':vt}
    
    nshem = 'N' if droplat >= 0 else 'S'
    ewhem = 'E' if droplon >= 0 else 'W'
    revcoilstr = 'Yes' if settings['revcoil'] else 'No'
    qualitystr = 'High' if settings['quality'] <= 1 else 'Low' if settings['quality'] >= 3 else 'Medium'
    
    comments = f"""Probe Type : AXCP
    MagVar date      : {dropdate:%Y/%m/%d}
    MagVar latitude  : {abs(droplat):7.3f} {nshem}
    MagVar longitude : {abs(droplon):8.3f} {ewhem}
    Quality        : {settings['quality']} ({qualitystr})
    Reverse coil   : {revcoilstr}
    Refresh rate   : {settings['refreshrate']} sec
    """
    
    writeedffile(output_file,None,None,None,data,comments,QC=False)
    

#MAIN
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    