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

import os
import sys
from collections import namedtuple


import numpy as np

import gsw
    


###################################################################################
#                               BITSTREAM PARSING                                 #
###################################################################################

# hexframes,times,depths,temps,conds,psals,next_buffer_ind = parse.parse_bitstream_to_profile(self.binary_buffer, binbufftimes, self.r7500_buffer, self.masks)
def parse_bitstream_to_profile(bitstream, times, r400_in, r7500_in, tempLUT, tcoeff, ccoeff, zcoeff):
    
    hexframes = [] # hexadecimal representation of frame
    proftime = [] #time (post-profile start) corresponding to each ob
    z = [] #depth
    T = [] #temperature
    C = [] #conductivity
    S = [] #salinity
    r400 = []
    r7500 = []
    
    #initializing fields for loop
    s = 0 #starting bit
    numbits = len(bitstream)
    
    #looping through bits until finished parsing
    while s < numbits - 32:
        
        foundMatch = False
        
        if s >= numbits - 32: #won't overrun bitstream
            break
            
        #pulling current segment
        frame = bitstream[s:s+32]
        
        #verifying that frame meets requirements, otherwise increasing start bit and rechecking
        if frame[0:2] != [1, 0] or not check_crc(frame) or not r7500_in[s] > 0:
            s += 1
            
        else: #good profile point
            ctime = times[s] #current profile time
            
            #converting frame to T/C/S/z
            Tint, Cint = convertFrameToInt(frame)
            cT, cC, cS, cz = convertIntsToFloats(Tint, Cint, ctime, tempLUT, tcoeff, ccoeff, zcoeff)
            
            #storing frame/time/profile data
            hexframes.append(binListToHex(frame))
            proftime.append(ctime)
            z.append(cz)
            T.append(cT)
            C.append(cC)
            S.append(cS)
            
            r400.append(r400_in[s])
            r7500.append(r7500_in[s])
            
            s += 32 #increase start bit by 32 to search for next frame

    # End parse bitstream
    return hexframes, proftime, z, T, C, S, r400, r7500, s
    

    
    
    

###################################################################################
#          FRAME CONVERSION TO TEMPERATURE/CONDUCTIVITY/SALINITY/DEPTH            #
###################################################################################

def convertFrameToInt(frame):
    """ Convert a frame to integer fields
    frame is a list of bits """
    Tint = binListToInt(frame[14:26])
    Cint = binListToInt(frame[2:14])
    
    return Tint, Cint
    

    
def convertIntsToFloats(Tint, Cint, time, tempLUT, tcoeff, ccoeff, zcoeff):
    """ Convert a list of integer data fields to observations (floats) """
    
    #depth from time
    z = dataconvert(time, zcoeff)
    
    #uncalibrated temperature and conductivity from integers
    if Tint >= 0 and Tint <= len(tempLUT)-1:
        Tuncal = tempLUT[Tint]
    else:
        Tuncal = np.NaN
        
    Cuncal = Cint * 60 / 4096
    
    #calibrated temperature and conductivity from uncalibrated
    T = dataconvert(Tuncal, tcoeff)
    C = dataconvert(Cuncal, ccoeff)
    
    #salinity from temperature/conductivity/depth
    S = gsw.SP_from_C(C,T,z) #assumes pressure (mbar) approx. equals depth (m)
        
    return T, C, S, z
        
        
        
        
def read_temp_LUT(filename):
    tempLUT = []
    with open(filename) as f:
        for line in f.readlines():
            cline = line.strip().split(',')
            if len(cline) >= 2:
                tempLUT.append(float(cline[1]))
                
    return tempLUT
    
    
    
###################################################################################
#                                 HEADER PARSING                                  #
###################################################################################   
    
    
    
def trim_header(bits_in):
    
    bits = bits_in.copy()
    bits[:25] = [True for _ in range(25)]
    
    last_index_pulse = 0 #index of last 1 bit preceded by 7 other 1 bits
    n_ones_in_last_25 = 0 #number of bit 1 in the last 25 bits
    
    for i,b in enumerate(bits):
        
        if b: #counting number of 1 bits in last 25 bits
            n_ones_in_last_25 += 1
            if i > 10:
                if np.sum(bits[i-7:i+1]) == 8: #if last 8 bits are all 1
                    last_index_pulse = i #save index
                    
        if i > 24: #first 25 bits (indices 0 to 24) are all ones
            if bits[i-25]: #subtract a bit from the total if the bit on the backside of the sliding window is a 1
                n_ones_in_last_25 -= 1
                
            if i >= 400 and n_ones_in_last_25 <= 20: # >0.5 sec into header and >1/5 of last 25 bits =0, PULSE END
                break
    
                
    header_bits = bits[last_index_pulse : last_index_pulse + 32*75] #pad with an extra three frames worth of bits to parse
    
    return header_bits
    
    
    
def initialize_axctd_metadata():
    
    #dict containing metadata (conversions, other header info) for the AXCTD
    axctd_metadata = {'tcoeff': [0,1,0,0], 'ccoeff': [0,1,0,0], 'zcoeff': [1,1,1,1], 'serial_no':None, 'probe_code':None, 'max_depth':None, 'misc':None, 'tcoeff_hex':['','','',''], 'ccoeff_hex':['','','',''], 'zcoeff_hex':['','','',''], 'tcoeff_valid':[False] * 4, 'ccoeff_valid':[False] * 4, 'zcoeff_valid':[False] * 4}
    
    return axctd_metadata
    
    
    
    
def parse_header(bits):
    
    #data format:
    #bits 0-1 = '10'
    #bits 2-9 = counter (increments 0-63 then '11111' followed by counter 0-7)
    #bits 10-25 = data 
    #bits 26-31 = CRC
    
    #logical list for whether each frame was detected
    counter_found = [False] * 72
    
    #whether or not conversion coefficients were successfully pulled for T/C/z
    tcoeff_valid = [False] * 4
    ccoeff_valid = [False] * 4
    zcoeff_valid = [False] * 4
    
    #dict containing metadata (conversions, other header info) for the AXCTD
    axctd_metadata = initialize_axctd_metadata()

    #parsing header into (up to) 72 frames
    frame_data = [None] * 72
    lastframe = -1
    Nbits = len(bits)
    s = 0
    
    past_seven = False
    while lastframe < 71 and s < Nbits-32:
        if bits[s:s+2] != [1, 0] or not check_crc(bits[s:s+32]):
            s += 1
            
        else: #frame
        
            #get counter, acknowledge it's found
            counter_bits = bits[s+2:s+10]
            if counter_bits[:5] == [1,1,1,1,1]:
                curcounter = binListToInt(counter_bits[5:]) + 64 #last 8 frames have counter 0-7 preceded by 11111
            else:
                curcounter = binListToInt(counter_bits) #first 64 frames just have counter
                
            if curcounter <= 71:
                counter_found[curcounter] = True
                lastframe = curcounter
                
                #convert frame data to hexidecimal format
                frame_data[curcounter] = binListToHex(bits[s+10:s+26])
            
            s += 32
    
            
    
    #pulling frame data into metadata dict
    if sum(counter_found[4:6]) == 2:
        axctd_metadata['serial_no'] = frame_data[4] + frame_data[5]
    
    if counter_found[6]: #frame 6: probe max depth
        axctd_metadata['max_depth'] = frame_data[6]
        
    if counter_found[7]: #frame 7: probe type (AXCTD should be A000)
        axctd_metadata['probe_code'] = frame_data[7]
    
    #frames 12-23: depth coefficients D3,D2,D1,D0 in 3-frame sequences
    for i,cf in enumerate(range(21,11,-3)): 
        if sum(counter_found[cf:cf+3]) == 3:
            axctd_metadata['zcoeff_hex'][i] = ''.join(frame_data[cf:cf+3])
    
    #frames 24-35: temp coefficients D3,D2,D1,D0 in 3-frame sequences
    for i,cf in enumerate(range(33,23,-3)): 
        if sum(counter_found[cf:cf+3]) == 3:
            axctd_metadata['tcoeff_hex'][i] = ''.join(frame_data[cf:cf+3])
            
    #frames 36-47: cond. coefficients D3,D2,D1,D0 in 3-frame sequences
    for i,cf in enumerate(range(45,35,-3)): 
        if sum(counter_found[cf:cf+3]) == 3:
            axctd_metadata['ccoeff_hex'][i] = ''.join(frame_data[cf:cf+3])
    
    #calculating conversion coefficients
    coeff_types = ['t','c','z']
    for coeff in coeff_types:
        for i in range(4):
            if axctd_metadata[coeff+'coeff_hex'][i] != '':
                chex = axctd_metadata[coeff+'coeff_hex'][i].upper().replace('B','+').replace('D','-')
                axctd_metadata[coeff+'coeff'][i] = int(chex[:9])/1E7 * 10**int(chex[9:])
                axctd_metadata[coeff+'coeff_valid'][i] = True
                
    #saving frame data to metadata
    axctd_metadata['frame_data'] = frame_data
    axctd_metadata['counter_found'] = counter_found
    
    return axctd_metadata
    
    
    
    


###################################################################################
#                   GENERAL LIST FORMAT DATA CONVERSION                           #
###################################################################################
    
#conversion: coefficients=C,  D_out = C[0] + C[1]*D_in + C[2]*D_in^2 + C[3]*D_in^3 + ...
def dataconvert(data_in,coefficients):
    output = 0
    for (i,c) in enumerate(coefficients):
        output += c*data_in**i
    return output
    
    


###################################################################################
#                   CYCLIC REDUNCANCY CHECK FOR FRAME                             #
###################################################################################

def check_crc(bits):
    
    divisor = [1,1,0,0,1,0,1]
    result = bits.copy()
    
    for k in range(26):
        if result[k]:
            for i in range(7):
                result[i+k] = result[i+k] != divisor[i]
            
    isGood = not sum(result) #returns 1 if sum is 0, otherwise 0
    
    return isGood

    

###################################################################################
#                         BINARY LIST / INTEGER CONVERSION                       #
###################################################################################


def binListToInt(binary):
    """ Convert a list of bits into a binary number """
    x = 0
    mask = 0x1 << (len(binary) - 1)
    for b in binary:
        if b:
            x |= mask
        mask >>= 1

    return x


def intToBinList(cInt, masklen):
    """ Convert a number into a list of bits with length
    at least masklen  """
    x = cInt
    i = 0
    bin_list = [0] * masklen
    while x:
        try:
            bin_list[i] = x & 1
        except IndexError: # if there are more bits than masklen
            bin_list.append(x & 1)
        x >>= 1
        i += 1

    bin_list.reverse()
    return bin_list

    
    
    
def binListToHex(binary):
    prof_starts = [i for i in range(0,len(binary),4)]
    prof_ends = [i for i in range(4,len(binary)+1,4)]
    
    letters = 'abcdef'
    
    hex_str = ''
    for s,e in zip(prof_starts, prof_ends):
        cbinnum = binListToInt(binary[s:e])
        if cbinnum <= 9: #convert decimal integer to hex digit
            cstr = str(cbinnum)
        else:
            cstr = letters[cbinnum-10]
            
        hex_str += cstr #append current hex digit to hex string
        
    return hex_str
    
    
    
            
