#!/usr/bin/env python3
# example.py

# This is an example script, demonstrating how to use PyAXBPS to process AXBT, AXCTD, 
#   and AXCP profiles from raw audio files
#
# The other imports in this file are for example purposes only and not required in your 
#   code to run PyAXBPS (datetime is imported within AXCPprocessor)
#
#
# NOTE- this package does not come with the sample probe data, you'll need to adjust the inputfile
#   variables to point to your own AXBT/AXCTD/AXCP WAV files

from PyAXBPS import AXBTprocessor, AXCTDprocessor, AXCPprocessor
from PyAXBPS import fileinteraction as io


from datetime import date, datetime 
import matplotlib.pyplot as plt #to plot output profiles
import numpy as np



#which of the probe types to process
do_axbt = False
do_axctd = False
do_axcp = True


#example date/time/position of drop (for file writing)
example_lat = 26.45 #26.45 degN
example_lon = -75.556 #75.556 degW
example_dropdatetime = datetime(2018,5,27,16,30,0) #27MAY2018 13:60 UTC

#additional example metadata required for BUFR files
example_id = 'NR23E' #aircraft tail number or station ID used for BUFR file, up to 9 characters
originatingcenter = 62 #BUFR WMO center (see BUFR tables online or in AXBPS GUI settings)






if __name__ == "__main__":

    ##################################################################################################################
    #                                             AXBT Processing Example                                            #
    ##################################################################################################################
    
    if do_axbt:
        print("Processing an AXBT file")
        inputfile = 'sample_files/AXBT_sample.WAV'
        
        #all of these values are defaults- you only need to specify them if you want something different
        timerange = [0,-1] #specify start and end times (in seconds) to process within the file, default 0 starts at the beginning and -1 means process to the end of the file
        audiochannel = -1 #specify a specific channel to process from a stereo wav file -> -1 sums across all channels
        
        #these are the default settings- you only need to include the settings you want to change in the settings dict. If you don't want to change any settings, you can pass an empty dict or not pass the argument at all
        settings = {'fftwindow':0.3, #window length in seconds to run the FFT for each temperature point
                    'minfftratio':0.5, #minimum signal-to-noise ratio (%) in AXBT temperature band for valid data
                    'minsiglev':65, #minimum signal level (dB) for valid data
                    'triggerfftratio':0.88, #minimum signal-to-noise ratio (%) to start profile collection
                    'triggersiglev':75, #minimum signal level (dB) to trigger profile collection
                    'tcoeff':[-40.0,0.02778,0.0,0.0], #frequency-to-temperature conversion coefficients
                    'zcoeff':[0.0,1.524,0.0,0.0],  #probe sink rate conversion coefficients
                    'flims':[1300,2800]} #frequency bounds for valid temperature observations
        
        #running the AXBT Processor
        axbt = AXBTprocessor.AXBT_Processor(inputfile, timerange=timerange, audiochannel=audiochannel, settings=settings)
        axbt.run()
        
        #pulling profile data after processing is complete
        temperature = np.array(axbt.temperature)
        depth = np.array(axbt.depth)
        time = np.array(axbt.time)
        frequency = np.array(axbt.fp)
        signal = np.array(axbt.Sp)
        snr = np.array(axbt.Rp)
        
        print("AXBT profile processing complete- writing output files!")
        
        
        #LOG DTA file- temperature only
        io.writelogfile('sample_AXBT_output.DTA', example_dropdatetime, depth, temperature, frequency, time, 'AXBT')
        
        
        #removing all NaNs
        isgood = [False if np.isnan(ct) else True for ct in temperature]
        temperature = temperature[isgood]
        time = time[isgood]
        frequency = frequency[isgood]
        signal = signal[isgood]
        snr = snr[isgood]
        depth = depth[isgood]
        
        edf_data = {'Time (s)':time, 'Frequency (Hz)':frequency, 'Sp (dB)':signal, 'Rp (%)':snr, 'Depth (m)':depth, 'Temperature (degC)':temperature} #format to give data to an EDF file- the dict key is the column header and the dict value should be a list with the column's contents
        
        #comments to include in the EDF file
        edf_comments = f"""Probe Type : AXBT
        Depth Coefficientd : {settings['zcoeff']}
        Temperature Coefficients : {settings['tcoeff']}
        Frequency Bounds (Hz) : {settings['flims']}
        FFT Window (sec): {settings['fftwindow']}
        Trigger signal level (dB)/ratio (%) : {settings['triggersiglev']} / {settings['triggerfftratio']}
        Minimum signal level (dB)/ratio (%) : {settings['minsiglev']} / {settings['minfftratio']}
        """
        
        
        #EDF file- variable format, custom fields and header with comments/drop metadata
        io.writeedffile('sample_AXBT_output.edf', example_dropdatetime, example_lat, example_lon, edf_data, edf_comments)
        
        #NVO file
        io.writefinfile('sample_AXBT_output.nvo', example_dropdatetime, example_lat, example_lon, 99, depth, temperature)
        
        
        #plotting temperature vs. depth
        fig = plt.figure() 
        fig.clear()
        ax = fig.add_axes([0.1,0.1,0.85,0.85])
        ax.plot(temperature, depth, 'r', linewidth=2)
        ax.set_xlabel('Temperature (^oC)')
        ax.set_ylabel('Depth (m)')
        ax.set_title('AXBT Example', fontweight="bold")
        ax.grid()
        ax.invert_yaxis()
        fig.savefig('AXBT_plot.png', format='png')
        plt.close('fig')
    
    
    
    
    
    
    
    
    
    
    ##################################################################################################################
    #                                             AXCTD Processing Example                                           #
    ##################################################################################################################
    
    
    if do_axctd:
    
        inputfile = 'sample_files/AXCTD_sample.WAV'
        
        #all of these values are defaults- you only need to specify them if you want something different
        timerange = [0,-1] #specify start and end times (in seconds) to process within the file, default 0 starts at the beginning and -1 means process to the end of the file
        settings = {}
        settings['minr400']          = 2.0         #header pulse power threshold
        settings['mindr7500']        = 1.5         #profile tone power threshold
        settings['deadfreq']         = 3000        #"quiet" frequency to use as baseline for power calcs
        settings['triggerrange']     = [30,-1]    #bounds for how soon to trigger profile after first pulse
        settings['mark_space_freqs'] = [400,800]   #mark and space bit frequencies (shouldn't be adjusted)
        settings['bitrate']          = 800         #symbol rate (shouldn't be adjusted)
        settings['bit_inset']        = 1           #how many PCM points to inset for demod bit identification
        settings['phase_error']      = 25          #phase error to consider for zero crossing trackind in demod
        settings['usebandpass']      = False       #use bandpass filter instead of lowpass?
        settings['refreshrate']      = 2.0         #how many seconds of data to process per loop
        
        #default conversion coefficients to use if headers can't be decoded
        settings["zcoeff_axctd"]     = [0.72, 2.76124, -0.000238007, 0]
        settings["tcoeff_axctd"]     = [-0.053328, 0.994372, 0.0, 0.0]
        settings["ccoeff_axctd"]     = [-0.0622192, 1.04584, 0.0, 0.0]
        
        settings["tlims_axctd"]      = [-10,50] #min/max temperature for a good datapoint
        settings["slims_axctd"]      = [-1,100] #min/max salinity for a good datapoint
        
            
        #running AXCTD Processor instance, timing
        print("Processing an AXCTD profile!")
        axctd = AXCTDprocessor.AXCTD_Processor(inputfile, timerange=timerange, user_settings=settings)
        axctd.run()
        print("AXCTD profile processing complete- writing output files!")
        
        
        
        edf_data = {'Time (s)': axctd.time, 'Frame (hex)':axctd.hexframes, 'Depth (m)':axctd.depth, 'Temperature (degC)':axctd.temperature, 'Conductivity (mS/cm)':axctd.conductivity, 'Salinity (PSU)':axctd.salinity}
        
        fs = axctd.f_s #sampling frequency of audio file
        
        metadata = axctd.metadata
        coeffs = {}
        coeffops = ['z','t','c']
        for c in coeffops:
            if sum(metadata[c + 'coeff_valid']) == 4:
                coeffs[c] = metadata[c+'coeff']
            else:
                coeffs[c] = metadata[c+'coeff_default']
        edf_comments = f"""Probe Type       :  AXCTD
            Serial Number      :  {metadata['serial_no'] if metadata['serial_no'] is not None else 'Not Provided'}
            Terminal Depth (m) :  {metadata['max_depth'] if metadata['max_depth'] is not None else 'Not Provided'}
            Depth Coeff. 1     :  {coeffs['z'][0]}
            Depth Coeff. 2     :  {coeffs['z'][1]}
            Depth Coeff. 3     :  {coeffs['z'][2]}
            Depth Coeff. 4     :  {coeffs['z'][3]}
            Pressure Pt Correction:  N/A
            Temp. Coeff. 1     :  {coeffs['t'][0]}
            Temp. Coeff. 2     :  {coeffs['t'][1]}
            Temp. Coeff. 3     :  {coeffs['t'][2]}
            Temp. Coeff. 4     :  {coeffs['t'][3]}
            Cond. Coeff. 1     :  {coeffs['c'][0]}
            Cond. Coeff. 2     :  {coeffs['c'][1]}
            Cond. Coeff. 3     :  {coeffs['c'][2]}
            Cond. Coeff. 4     :  {coeffs['c'][3]}
            Sampling frequency (fs): {fs} Hz
            Audio file length: {axctd.numpoints/fs} sec
            400 Hz pulse start: {axctd.firstpulse400/fs} sec
            7500 Hz tone start: {axctd.profstartind/fs} sec
            """
            
        #EDF file- variable format, custom fields and header with comments/drop metadata
        io.writeedffile('sample_AXCTD_output.edf', example_dropdatetime, example_lat, example_lon, edf_data, edf_comments)
        
        
        #pulling profile data after processing is complete
        depth = np.array(axctd.depth)
        temperature = np.array(axctd.temperature)
        salinity = np.array(axctd.salinity)
        
        #removing all NaNs
        isgood = [False if np.isnan(ct*cz*cs) else True for ct,cz,cs in zip(temperature,depth,salinity)]
        depth = depth[isgood]
        temperature = temperature[isgood]
        salinity = salinity[isgood]
        
        
        #plotting temperature and salinity vs. depth
        fig = plt.figure() 
        fig.clear()
        ax = fig.add_axes([0.1,0.1,0.85,0.8])
        axS = ax.twiny()
        
        ax.set_xlabel('Temperature (^oC)')
        ax.set_ylabel('Depth (m)')
        axS.set_xlabel('Salinity (PSU)', fontsize=12)
        ax.xaxis.label.set_color("red") #temperature axis red
        ax.tick_params(axis='x', colors='red')
        axS.xaxis.label.set_color("blue") #salinity axis blue
        axS.tick_params(axis='x', colors='blue')
        ax.invert_yaxis() #only need to invert one zxis since twiny causes the other to copy it
        
        
        ax.plot(temperature, depth, 'r', linewidth=2)
        axS.plot(salinity, depth, 'b', linewidth=2)
        
        ax.set_title('AXCTD Example', fontweight="bold")
        ax.grid()
        fig.savefig('AXCTD_plot.png', format='png')
        plt.close('fig')
    
    
    
    
    
    
    
    
    ##################################################################################################################
    #                                             AXCP Processing Example                                            #
    ##################################################################################################################
    
    if do_axcp:
        
        inputfile = 'sample_files/AXCP_sample.WAV'
        
        # date/lat/lon are optional but necessary to calculate the horizontal and vertical components
        #   of Earth's magnetic field at the drop location, which go into the current speed calculations
        #   Having a somewhat accurate estimate of position (within a few degrees) is necessary for an accurate profile
        # This info is also used to calculate declination and convert the U/V components of the currents from 
        #   degrees magnetic to true
        
        dropdate = example_dropdatetime.date()
        droplat = example_lat
        droplon = example_lon
        
        #all of these values are defaults- you only need to specify them if you want something different
        timerange = [0,-1] #specify start and end times (in seconds) to process within the file, default 0 starts at the beginning and -1 means process to the end of the file
        settings = {'quality': 1, #1=high/slow, 2=moderate quality/speed, 3=low quality/fast
                    'revcoil': 0, #whether the probe's coil is reversed (false unless you know about it!)
                    'refreshrate': 1.0, #how many seconds of data to process at once
                    'spindown_detect_rt': 1, #whether or not to auto-detect probe spindown and stop processing at file end
                    'temp_mode': 2, #1-> use zero-crossings to ID temperature frequency, 2->use FFT (slower, more accurate)
                    'tempfftwindow': 1.0} #FFT window for temperature if temp_mode = 2 (tempfftwindow must be <= refreshrate)
        
        print("Processing an AXCP profile!")
        axcp = AXCPprocessor.AXCP_Processor(inputfile, timerange=timerange, lat=droplat, lon=droplon, dropdate=dropdate, settings=settings)
        axcp.run()
        
        
        
        print("AXCP profile processing complete- writing output files!")
        
        edf_comments = f"""Probe Type       :  AXCP
            Spinup Time (s)         : {axcp.tspinup}
            Spindown Time (s)       : {axcp.tspindown}
        """
        
        edf_data = {'Time (s)': axcp.TIME, 'Rotation Rate (Hz)':axcp.ROTF, 'Depth (m)':axcp.DEPTH, 'Temperature (C)':axcp.TEMP, 'Zonal Current (m/s)':axcp.U_TRUE, 'Meridional Current (m/s)':axcp.V_TRUE}
        
        
        
        #EDF file- variable format, custom fields and header with comments/drop metadata
        io.writeedffile('sample_AXCP_output.edf', example_dropdatetime, example_lat, example_lon, edf_data, edf_comments)
        
        
        
        #plotting temperature and currents vs. depth
        fig = plt.figure() 
        fig.clear()
        ax = fig.add_axes([0.1,0.1,0.85,0.8])
        axC = ax.twiny()
        
        ax.set_xlabel('Temperature (^oC)')
        ax.set_ylabel('Depth (m)')
        axC.set_xlabel('Current Speed (m/s)', fontsize=12)
        ax.xaxis.label.set_color("red") #temperature axis red
        ax.tick_params(axis='x', colors='red')
        axC.xaxis.label.set_color("black") #current axis black
        axC.tick_params(axis='x', colors='black')
        ax.invert_yaxis() #only need to invert one zxis since twiny causes the other to copy it
        
        
        axC.plot(axcp.U_TRUE, axcp.DEPTH, 'b', linewidth=1, label='U')
        axC.plot(axcp.V_TRUE, axcp.DEPTH, 'g', linewidth=1, label='V')
        ax.plot(axcp.TEMP, axcp.DEPTH, 'r', linewidth=2, label='Temperature')
        
        ax.set_title('AXCP Example', fontweight="bold")
        ax.grid()
        axC.legend()
        fig.savefig('AXCP_plot.png', format='png')
        plt.close('fig')
    
    
    
    
    
    
    
    
    
    
    