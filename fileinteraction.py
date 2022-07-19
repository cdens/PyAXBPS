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

import numpy as np
from datetime import date, datetime, timedelta
import chardet

    
    
    
    
def writelogfile(logfile, dropdatetime, depth, tempC, frequency, time, probetype):
    
    header = f""" Probe Type = {probetype}     
       Date = {dropdatetime:%Y/%m/%d}
       Time = {dropdatetime:%H:%M:%S}
 
    Time     Depth    Frequency    (C)       (F) 
    """
    
    with open(logfile,'w') as f_out:
        f_out.write(header + "\n")
        
        for (ct,cf,cd,cT) in zip(time,frequency,depth,tempC):
            
            if np.isnan(cf*cT*cd) or cf <= 0: #invalid point
                cline = f'{ct:7.1f}{-10:10.1f}{0:11.1f}    ******    ******'
                
            else: #valid data point
                cTf = cT*9/5 + 32 #degF
                cline = f'{ct:7.1f}{cd:10.1f}{cf:11.1f}{cT:10.2f}{cTf:10.2f}'
                
            f_out.write(cline + "\n")
        
            



def writeedffile(edffile,dropdatetime,lat,lon,data,comments,QC=False):
    
    
    with open(edffile,'w') as f_out:
    
        #writing header, date and time, drop # (bad value)
        f_out.write("// This is an air-launched probe EXPORT DATA FILE  (EDF)\n")
        f_out.write("// File generated with the Airborne eXpendable Buoy Processing System (AXBPS)\n")
        f_out.write(f"Date of Launch:  {datetime.strftime(dropdatetime,'%m/%d/%y')}\n")
        f_out.write(f"Time of Launch:  {datetime.strftime(dropdatetime,'%H:%M:%S')}\n")
        
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

    
            
            




#write data to FIN file
def writefinfile(finfile,cdtg,lat,lon,num,depth,temperature,salinity=None,U=None,V=None):
    
    hasSal = False
    hasCurrents = False
    if salinity is not None: 
        hasSal = True
    elif U is not None and V is not None:
        hasCurrent = True
    
        
    with open(finfile,'w') as f_out:
    
        dayofyear = date.toordinal(date(cdtg.year,cdtg.month,cdtg.day)) - date.toordinal(date(cdtg.year-1,12,31))

        #formatting latitude string
        if lat >= 0:
            latsign = ' '
        else:
            latsign = '-'
            lat = abs(lat)

        #formatting longitude string
        if lon >= 0:
            lonsign = ' '
        else:
            lonsign = '-'
            lon = abs(lon)

        #writing header data
        f_out.write(f"{cdtg.year}   {dayofyear:03d}   {datetime.strftime(cdtg,'%H%M')}   {latsign}{lat:06.3f}   {lonsign}{lon:07.3f}   {num:02d}   6   {len(depth)}   0   0   \n")
        
        #creating list of data points
        Npts = len(depth)
        outdata = []
        if hasSal:
            obsperline = 3
            nobtypes = 3
            for (t,d,s) in zip(temperature,depth,salinity):
                outdata.extend([f"{t: 8.3f}",f"{d: 8.1f}",f"{s: 8.3f}"])
        elif hasCurrents:
            obsperline = 3
            nobtypes = 4
            for (t,d,u,v) in zip(temperature,depth,U,V):
                outdata.extend([f"{t: 8.3f}",f"{d: 8.1f}",f"{u: 8.3f}",f"{v: 8.3f}"])
        else:
            nobtypes = 2
            obsperline = 5
            for (t,d) in zip(temperature,depth):
                outdata.extend([f"{t: 8.3f}",f"{d: 8.1f}"])
                
        pointsperline = nobtypes*obsperline

        #writing profile data
        i = 0
        while i < Npts:
            if i+pointsperline < Npts:
                pointstopull = pointsperline
            else:
                pointstopull = Npts - i
                
            line = "".join(outdata[i:i+pointstopull]) + "\n"
            f_out.write(line)
            i += pointstopull
            

            
            
            
            


#cdtg is formatted as datetime.datetime object
def writebufrfile(bufrfile, cdtg, lon, lat, identifier, originatingcenter, depth, temperature, salinity=None, U=None, V=None, optionalinfo=None):

    binarytype = 'big'  # big-endian
    reserved = 0
    version = 4  # BUFR version number (3 or 4 supported)

    # section 1 info
    if version == 3:
        sxn1len = 18
    elif version == 4:
        sxn1len = 22

    mastertable = 0  # For standard WMO FM 94-IX BUFR table
    originatingsubcenter = 0
    updatesequencenum = 0  # first iteration
    
    datacategory = 31  # oceanographic data (Table A)
    datasubcategory = 3 #bathy (JJVV) message
    versionofmaster = 32
    versionoflocal = 0
    yearofcentury = int(cdtg.year - 100 * np.floor(cdtg.year / 100))  # year of the current century

    
    # section 2 info
    if optionalinfo:
        hasoptionalinfo = True
        hasoptionalsectionnum = int('10000000', 2)
        sxn2len = len(optionalinfo) + 4
    else:
        sxn2len = 0
        hasoptionalinfo = False
        hasoptionalsectionnum = int('00000000', 2)
        
        
    hasSal = False
    hasCurrents = False
    if salinity is not None: 
        hasSal = True
    elif U is not None and V is not None:
        hasCurrent = True
        
        
    # Section 3 info
    sxn3len = 25  # 3 length + 1 reserved + 2 numsubsets + 1 flag + 9*2 FXY = 25 octets
    if hasSal:
        sxn3len = 27 #add 2 octets for salinity FXY code
    numdatasubsets = 1
    
    # whether data is observed, compressed (bits 1/2), bits 3-8 reserved (=0)
    s3oct7 = int('10000000', 2)
    
    #replication factor explanation:
    #F=1 for replication, X=number of descriptors (2 for d/T, 3 for d/T/S), Y=#times replicated (1)
    
    fxy_all = [int('0000000100001011',2), #0/01/011: identifier (72b, 9 bytes ascii)
    int('1100000100001011',2), #3/01/011: year/month/day (12b/4b/6b = 22b total)
    int('1100000100001100',2), #3/01/012: hour/minute (5b/6b = 11b total)
    int('1100000100010111',2), #3/01/023: lat/lon coarse accuracy (15b lat/16b lon = 31b total)
    int('0000001000100000',2), #0/02/032: digitization indicator (2b, 00)
    int('0100001000000000',2), #1/02/000: replication of 2 descriptors
    int('0001111100000010',2), #0/31/002: extended delayed decriptor replication factor (16b, # obs)
    int('0000011100111110',2), #0/07/062: depth (17b, m*10)
    int('0001011000101011',2),] #0/22/043: temperature (15b, degK*100)    
    
    if hasSal: #change replicator and append salinity 
        fxy_all[5] = int('0100001100000000',2) #change replicator (index 5) to delayed rep of 3 descriptors, code 1/03/000
        fxy_all.append(int('0001011000111110',2)) #0/22/062: salinity (14b, 100*PPT (=PSU))
    elif hasCurrents:
        fxy_all[5] = int('0100010000000000',2) #change replicator (index 5) to delayed rep of 4 descriptors, code 1/04/000
        fxy_all.append(int('0001011000000100',2)) #0/22/004: current direction (9b, degrees True/towards)
        fxy_all.append(int('0001011000011111',2)) #0/22/031: current speed (13b, 100*vel in m/s)
        

    # Section 4 info (data)
    identifier = identifier[:9] #concatenates identifier if necessary
    idtobuffer = 9 - len(identifier) # padding and encoding station identifier
    id_utf = identifier.encode('utf-8')
    for i in range(0, idtobuffer):
        id_utf = id_utf + b'\0' #pads with null character (\0 in python)

    # setting up data array
    bufrarray = ''

    # year/month/day (3,01,011)
    bufrarray = bufrarray + format(cdtg.year, '012b')  # year (0,04,001)
    bufrarray = bufrarray + format(cdtg.month, '04b')  # month (0,04,002)
    bufrarray = bufrarray + format(cdtg.day, '06b')  # day (0,04,003)

    # hour/minute (3,01,012)
    bufrarray = bufrarray + format(cdtg.hour, '05b')  # hour (0,04,004)
    bufrarray = bufrarray + format(cdtg.minute, '06b')  # min (0,04,005)
    
    # lat/lon (3,01,023)
    bufrarray = bufrarray + format(int(np.round((lat * 100)) + 9000), '015b')  # lat (0,05,002)
    bufrarray = bufrarray + format(int(np.round((lon * 100)) + 18000), '016b')  # lon (0,06,002)
    

    # temperature-depth profile (3,06,001)
    bufrarray = bufrarray + '00'  # indicator for digitization (0,02,032): 'values at selected depths' = 0
    bufrarray = bufrarray + format(len(temperature),'016b')  # delayed descripter replication factor(0,31,001) = length
    
        
    #converting temperature, salinity (as req'd.), and depth and writing
    if hasSal:
        for t,d,s in zip(temperature,depth,salinity):
            d_in = int(np.round(d*10)) # depth (0,07,062)
            bufrarray = bufrarray + format(d_in,'017b')
            t_in = int(np.round(100 * (t + 273.15)))  # temperature (0,22,043)
            bufrarray = bufrarray + format(t_in,'015b')
            s_in = int(np.round(100 * s))  # salinity (0,22,062)
            bufrarray = bufrarray + format(s_in,'014b')
        
    elif hasCurrents:
        for t,d,u,v in zip(temperature,depth,U,V):
            d_in = int(np.round(d*10)) # depth (0,07,062)
            bufrarray = bufrarray + format(d_in,'017b')
            t_in = int(np.round(100 * (t + 273.15)))  # temperature (0,22,043)
            bufrarray = bufrarray + format(t_in,'015b')
            
            #calculating direction and velocity
            curvel = np.sqrt(u**2 + v**2)
            curdir = 90 - (180/np.pi) * np.arctan2(v,u) #degrees true, north 0/clockwise
            
            curdir_in = int(np.round(curdir))  #0/22/004: current direction (9b, degrees True/towards)
            bufrarray = bufrarray + format(curdir_in,'09b')
            curvel_in = int(np.round(100 * curvel))  #0/22/031: current speed (13b, 100*vel in m/s)
            bufrarray = bufrarray + format(curvel_in,'013b')
        
        
    else:
        for t,d in zip(temperature,depth):
            d_in = int(np.round(d*10)) # depth (0,07,062)
            bufrarray = bufrarray + format(d_in,'017b')
            t_in = int(np.round(100 * (t + 273.15)))  # temperature (0,22,042)
            bufrarray = bufrarray + format(t_in,'015b')

    #padding zeroes to get even octet number, determining total length
    bufrrem = len(bufrarray)%8
    for i in range(8-bufrrem):
        bufrarray = bufrarray + '0'
    bufrarraylen = int(len(bufrarray)/8)
    sxn4len = 4 + 9 + bufrarraylen  # length/reserved + identifier + bufrarraylen (lat/lon, dtg, t/d profile)
    
    # total length of file in octets
    num_octets = 8 + sxn1len + sxn2len + sxn3len + sxn4len + 4  # sxn's 0 and 5 always have 8 and 4 octets, respectively
        
    # writing the file
    with open(bufrfile, 'wb') as bufr:

        # Section 0 (indicator)
        bufr.write(b'BUFR')  # BUFR
        bufr.write(num_octets.to_bytes(3, byteorder=binarytype, signed=False))  # length (in octets)
        bufr.write(version.to_bytes(1, byteorder=binarytype, signed=False))

        # Section 1 (identifier) ***BUFR version 3 or 4 ****
        bufr.write(sxn1len.to_bytes(3, byteorder=binarytype, signed=False))
        bufr.write(mastertable.to_bytes(1, byteorder=binarytype, signed=False))
        if version == 3:
            bufr.write(originatingsubcenter.to_bytes(1, byteorder=binarytype, signed=False))
            bufr.write(originatingcenter.to_bytes(1, byteorder=binarytype, signed=False))
        elif version == 4:
            bufr.write(originatingcenter.to_bytes(2, byteorder=binarytype, signed=False))
            bufr.write(originatingsubcenter.to_bytes(2, byteorder=binarytype, signed=False))
        bufr.write(updatesequencenum.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(hasoptionalsectionnum.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(datacategory.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(datasubcategory.to_bytes(1, byteorder=binarytype, signed=False))
        if version == 4:
            bufr.write(datasubcategory.to_bytes(1, byteorder=binarytype, signed=False)) #write again for local data subcategory in version 4
        bufr.write(versionofmaster.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(versionoflocal.to_bytes(1, byteorder=binarytype, signed=False))
        if version == 3:
            bufr.write(yearofcentury.to_bytes(1, byteorder=binarytype, signed=False))
        elif version == 4:
            bufr.write(cdtg.year.to_bytes(2, byteorder=binarytype, signed=False))
        bufr.write(cdtg.month.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(cdtg.day.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(cdtg.hour.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(cdtg.minute.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(int(0).to_bytes(1, byteorder=binarytype, signed=False)) #seconds for v4, oct18 = 0 (reserved) for v3

        # Section 2 (optional)
        if hasoptionalinfo:
            bufr.write(sxn2len.to_bytes(3, byteorder=binarytype, signed=False))
            bufr.write(reserved.to_bytes(1, byteorder=binarytype, signed=False))
            bufr.write(optionalinfo)

        # Section 3 (Data description)
        bufr.write(sxn3len.to_bytes(3, byteorder=binarytype, signed=False))
        bufr.write(reserved.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(numdatasubsets.to_bytes(2, byteorder=binarytype, signed=False))
        bufr.write(s3oct7.to_bytes(1, byteorder=binarytype, signed=False))
        for fxy in fxy_all:
            bufr.write(fxy.to_bytes(2, byteorder=binarytype, signed=False))

        # Section 4
        bufr.write(sxn4len.to_bytes(3, byteorder=binarytype, signed=False))
        bufr.write(reserved.to_bytes(1, byteorder=binarytype, signed=False))
        bufr.write(id_utf)

        #writing profile data- better option- converts and writes 1 byte at a time
        for cbit in np.arange(0,len(bufrarray),8):
            bufr.write(int(bufrarray[cbit:cbit+8],2).to_bytes(1,byteorder=binarytype,signed=False)) #writing profile data

        # Section 5 (End)
        bufr.write(b'7777')

        