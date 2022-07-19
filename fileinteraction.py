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



#read raw temperature/depth profile from LOGXX.DTA file
def readlogfile(logfile):

    with open(logfile,'r') as f_in:
        depth = []
        temperature = []

        for line in f_in:

            line = line.strip().split() #get rid of \n character, split by spaces

            try:
                curfreq = np.double(line[2])
                cdepth = np.double(line[1])
                ctemp = np.double(line[3])
            except:
                curfreq = np.NaN

            if ~np.isnan(curfreq) and curfreq != 0: #if it is a valid frequency, save the datapoint
                depth.append(cdepth)
                temperature.append(ctemp)
    
    #convert to numpy arrays
    depth = np.asarray(depth)
    temperature = np.asarray(temperature)
    
    return [temperature,depth]
    
    
    
    
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
        
        
        
    
            
    
    
#code to parse a variety of date formats to python datetime variable
def parse_date(line):
    day = month = year = False #initializing
    
    if "/" in line and len(line) <= 8: #mm/dd/yy format
        line = line.split('/')
        month = int(line[0])
        day = int(line[1])
        year = int(line[2]) + 2000
    elif "/" in line and len(line) <= 10: #mm/dd/yyyy, or yyyy/mm/dd (assuming not dd/mm/yyyy)
        line = line.split('/')
        if len(line[0]) == 4:
            year = int(line[0])
            month = int(line[1])
            day = int(line[2])
        elif len(line[2]) == 4:
            month = int(line[0])
            day = int(line[1])
            year = int(line[2])
            
    elif "-" in line and len(line) <= 8: #mm-dd-yy format
        line = line.split('-')
        month = int(line[0])
        day = int(line[1])
        year = int(line[2]) + 2000
    elif "-" in line and len(line) <= 10: #mm-dd-yyyy, or yyyy-mm-dd (assuming not dd-mm-yyyy)
        line = line.split('-')
        if len(line[0]) == 4:
            year = int(line[0])
            month = int(line[1])
            day = int(line[2])
        elif len(line[2]) == 4:
            year = int(line[2])
            month = int(line[1])
            day = int(line[0])
    
    else: #trying YYYYMMDD format instead
        year = int(line[:4])
        month = int(line[4:6])
        day = int(line[6:8])
        
    return year,month,day
    
    
    
    
#parsing a variety of lats/lons into valid floats with N and E positive
def parse_lat_lon(line):
    line = line.strip().lower() #sanitizing
    
    #identifying hemisphere
    hemsign = 1 #assuming N/E hemisphere unless proven otherwise
    if 'n' in line or 's' in line or 'e' in line or 'w' in line: #hemisphere given
        if 's' in line or 'w' in line:
            hemsign = -1
        line = line[:-1] #drop hemisphere from line
            
    else: #lat/lon is positive or negative for E/N or W/S hemisphere respectively
        if line[0] == '-':
            hemsign = -1
            line = line[1:]
    
    #sanitize and split line
    line_split = line.replace(' ',':')
    line_split = line_split.split(':')
        
    #calculating lon/lat in decimal degrees
    value = 0
    for i,valstr in enumerate(line_split):
        value += hemsign*float(valstr)/60**i
    
    return value
    
    
    
    
    
    
    
def readdatfile(datfile):
    
    data = {'temperature':[], 'salinity':[], 'depth':[]} #initializing data lists
    
    with open(datfile) as f_in:
        
        for i,cline in enumerate(f_in): 
            cline = [i for i in cline.strip().split(' ') if len(i) > 0] #space delimited
            
            if i < 3: #first 3 lines are trash
                pass
            elif i == 3: #header
                dtgstr = cline[0]+cline[1]
                if len(dtgstr) == 12:
                    dtgformat = '%Y%m%d%H%M%S'
                else:
                    dtgformat = '%Y%m%d%H%M' #correct to best guess of format
                    dtgstr = dtgstr[:10]
                    
                dropdatetime = datetime.strptime(dtgstr,dtgformat)
                lat = float(cline[2])
                lon = float(cline[3])
            
            else:
                if len(cline) > 1:
                    data['depth'].append(float(cline[0]))
                    data['temperature'].append(float(cline[1]))
                    
                    try: #try/except so dat format is compatible with z/T columns only
                        cpsal = float(cline[2])
                    except IndexError:
                        cpsal = np.NaN
                        
                    if cpsal >= 0:
                        data['salinity'].append(cpsal)
                    else:
                        data['salinity'].append(np.NaN)
        
    
    for cfield in ['depth','temperature','salinity']:  
        data[cfield] = np.asarray(data[cfield])   
        
    return dropdatetime,lat,lon,data
    
    
    
    
#TODO: understand format of 1st 4 lines -> what is second time, what is 12, what are first two lines for?
# ****0000005762****
# SOFX01 KWBC 172308

# 20210817 225634  17.230  -79.962 N43RF AL072021 12 2021-229-23:08:29
    
def writedatfile(datfile,dropdatetime,lat,lon,headerstart,tailnum,missionnum,depth,temperature,salinity=None):
    
    with open(datfile,'w') as f_out:
        
        header = headerstart + f"\n{dropdatetime:%Y%m%d} {dropdatetime:%H%M%S} {lat:7.3f} {lon:8.3f} {tailnum} {missionnum} 12 {datetime.utcnow():%Y-%j-%H:%M:%S}\n"
        
        f_out.write(header)
        
        if salinity is None:
            salinity = [np.NaN]*len(depth)
        if len(salinity) != len(depth): #now salinity is definitely a list, can call len on it
            salinity = [np.NaN]*len(depth)
            
        salinity = [-9.99 if np.isnan(i) else i for i in salinity] #replace any NaNs with negative value
        
        for (cd,cT,cS) in zip(depth,temperature,salinity):
            if not np.isnan(cd):
                f_out.write(f'\n{cd:6.1f}{cT:6.2f}{cS:6.2f}\n')
    
        f_out.write('\n\n') #ends with two more empty lines
    
    
    
    
    

    
def readedffile(edffile):
    
    encoding = 'utf-8'
    
    lon = lat = day = month = year = hour = minute = second = False #variables will be returned as 0 if unsuccessfully parsed
    
    data = {'depth':[],'temperature':[],'salinity':[]} #initializing each data field as [] so if it isn't empty upon function completion the user knows the file had a match for that field
    
    fields = ['depth','temperature','salinity']
    fieldcolumns = [0,1,-1]
    
    with open(edffile,'rb') as f_in:
        
        for cline in f_in:
            try:
                cline = cline.decode(encoding).strip()
            except:
                fileinfo = chardet.detect(cline)
                encoding = fileinfo['encoding']
                cline = cline.decode(encoding).strip()
                
            try:
                if ":" in cline: #input parameter- parse appropriately
                    line = cline.strip().split(':')
                    
                    if "time" in line[0].lower(): #assumes time is in "HH", "HH:MM", or "HH:MM:SS" format
                        hour = int(line[1].strip())
                        minute = int(line[2].strip())
                        second = int(line[3].strip())
                        
                    elif "date" in line[0].lower(): #TODO: fix date and time reading
                        line = line[1].strip() #should be eight digits long
                        year,month,day = parse_date(line)
                            
                    
                    elif "latitude" in line[0].lower(): 
                        startind = cline.find(':')
                        lat = parse_lat_lon(cline[startind+1:])
                            
                    elif "longitude" in line[0].lower():
                        startind = cline.find(':')
                        lon = parse_lat_lon(cline[startind+1:])
                        
                    elif "field" in line[0].lower(): #specifying which column contains temperature, depth, salinity (optional)
                        matched = False
                        curlinefield = line[1].lower()
                        curfieldcolumn = int(line[0].strip()[-1]) - 1
                        for i,field in enumerate(fields):
                            if field in curlinefield:
                                fieldcolumns[i] = curfieldcolumn
                                matched = True
                        if not matched:
                            fields.append(curlinefield)
                            fieldcolumns.append(curfieldcolumn)
                            data[curlinefield] = []
                        
                                
                #space or tab delimited data (if number of datapoints matches fields and line doesnt start with //)
                elif cline[:2] != '//': 
                    curdata = cline.strip().split()
                    #if salinity is included, #of columns should match number of fields, if not, # columns should be # fields minus 1 (no salinity data)
                    if (fieldcolumns[2] >= 0 and len(curdata) == len(fields)) or (fieldcolumns[2] == -1 and len(curdata) == len(fields)-1):
                        for i,cfield in enumerate(fields):
                            data[cfield].append(curdata[fieldcolumns[i]])
                    
            except (ValueError, IndexError, AttributeError):
                pass
    
    #determining datetime
    try:
        dropdatetime = datetime(year,month,day,hour,minute,second)
    except:
        dropdatetime = False #return false for datetime if inputs are invalid
        
    #checking if each field can be converted from string to float and doing so as possible
    for field in data.keys():
        try:
            curdata_float = [float(val) for val in data[field]] #if this doesn't raise an error, the float conversion worked
            data[field] = curdata_float
        except (ValueError, TypeError): #raised if data can't be converted to float
            pass
    
    for cfield in ['depth','temperature','salinity']:  
        data[cfield] = np.asarray(data[cfield])   
        
    return data,dropdatetime,lat,lon

    
    
    

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

    
    
    
#read data from JJVV file (AXBT and/or temperature data only)
def readjjvvfile(jjvvfile):
    with open(jjvvfile,'r') as f_in:

        depth = []
        temperature = []

        line = f_in.readline()
        line = line.strip().split()
        
        #date and time info
        datestr = line[1]
        day = int(datestr[:2])
        month = int(datestr[2:4])
        yeardigit = int(datestr[4])
        # time = int(line[2][:4])
        hour = int(line[2][:2])
        minute = int(line[2][2:4])
        
        #determining year (drop must have been within last 10 years)
        curyear = datetime.utcnow().year
        decade = int(np.floor(curyear/10)*10)
        if yeardigit + decade > curyear:
            decade -= 1
        year = decade + yeardigit
        
        try:
            dropdatetime = datetime(year,month,day,hour,minute,0)
        except:
            dropdatetime = False #return false for datetime if inputs are invalid

        #latitude and longitude
        latstr = line[3]
        lonstr = line[4]
        quad = int(latstr[0])
        lat = np.double(latstr[1:3]) + np.double(latstr[3:])/10**(len(latstr)-3)
        lon = np.double(lonstr[:3]) + np.double(lonstr[3:])/10**(len(lonstr)-3)
        if quad == 3:#hemisphere (if quad == 1, no need to change anything)
            lat = -1*lat
        elif quad == 5:
            lon = -1*lon
            lat = -1*lat
        elif quad == 7:
            lon = -1*lon
            
        lastdepth = -1
        hundreds = 0
        l = 0

        identifier = 'UNKNOWN'

        for line in f_in:
            l = l + 1

            #removing non-data entry from first column, 2nd line of JJVV
            line = line.strip().split()
            if l == 1: line = line[1:]

            for curentry in line:

                try:
                    int(curentry) #won't execute if curentry has non-numbers in it (e.g. the current entry is the identifier)

                    if int(curentry[:3]) == 999 and int(curentry[3:])*100 == hundreds + 100:
                        hundreds = hundreds + 100
                    else:
                        if int(curentry[:2]) + hundreds != lastdepth:
                            cdepth = int(curentry[:2]) + hundreds
                            lastdepth = cdepth
                            depth.append(cdepth)
                            temperature.append(np.double(curentry[2:])/10)

                except: identifier = curentry
    
    #converting to numpy arrays
    depth = np.asarray(depth)
    temperature = np.asarray(temperature)
    
    #correcting negative temperatures (encoded as abs(T) + 50)
    ind = temperature >= 50
    temperature[ind] = -1*(temperature[ind]-50)
    
    return [temperature,depth,dropdatetime,lat,lon,identifier]




#write data to JJVV file (AXBT/temperature only)
def writejjvvfile(jjvvfile,temperature,depth,cdtg,lat,lon,identifier,isbtmstrike):
    
    #open file for writing
    with open(jjvvfile,'w') as f_out:
    
        #first line- header information
        if lon >= 0 and lat >= 0:
            quad = '1'
        elif lon >= 0 and lat < 0:
            quad = '3'
        elif lon < 0 and lat >= 0:
            quad = '7'
        else:
            quad = '5'
            
        #getting ones digit from year
        yeardigit = cdtg.year - int(np.floor(cdtg.year/10)*10)

        f_out.write(f"JJVV {datetime.strftime(cdtg,'%d%m')}{str(yeardigit)} {datetime.strftime(cdtg,'%H%M')}/ {quad}{int(abs(lat)*1000):05d} {int(abs(lon)*1000):06d} 88888\n")

        #create a list with all of the entries for the file
        filestrings = []
        filestrings.append('51099')
        hundreds = 0
        i = 0
        lastdepth = -1

        # appending data to list, adding hundreds increment counters where necessary (while loop necessary in case a gap > 100m exists)
        while i < len(depth):
            cdepth = round(depth[i])
            ctemp = temperature[i]
            
            if ctemp < 0: #correcting if value is negative
                ctemp = np.abs(ctemp) + 50
            
            if cdepth - hundreds > 99:  # need to increment hundreds counter in file
                hundreds = hundreds + 100
                filestrings.append(f'999{int(hundreds/100):02d}')
            else:
                if cdepth - lastdepth >= 1 and cdepth - hundreds <= 99:  # depth must be increasing, not outside of current hundreds range
                    lastdepth = cdepth
                    filestrings.append(f"{int(round(cdepth-hundreds)):02d}{int(round(ctemp,1)*10):03d}")
                    
                i += 1

        if isbtmstrike: #note if the profile struck the bottom
            filestrings.append('00000')

        identifier = identifier[:5] #concatenates if larger than 5 digits
        filestrings.append(identifier) #tack identifier onto end of file entries

        #writing all data to file
        i = 0
        while i < len(filestrings):
            if i == 0 and len(filestrings) >= 6: #first line has six columns (only if there is enough data)
                line = f"{filestrings[i]} {filestrings[i+1]} {filestrings[i+2]} {filestrings[i+3]} {filestrings[i+4]} {filestrings[i+5]}\n"
                i += 6
            elif i+5 < len(filestrings): #remaining full lines have five columns
                line = f"{filestrings[i]} {filestrings[i+1]} {filestrings[i+2]} {filestrings[i+3]} {filestrings[i+4]}\n"
                i += 5
            else: #remaining data on last line
                line = ''
                while i < len(filestrings):
                    line += filestrings[i]
                    if i == len(filestrings) - 1:
                        line += '\n'
                    else:
                        line += ' '
                    i += 1
            f_out.write(line)
            



#read data from FIN file
def readfinfile(finfile, hasSal=False):
    
    with open(finfile,'r') as f_in:

        line = f_in.readline()
        line = line.strip().split()

        #pulling relevant header information
        year = int(line[0])
        dayofyear = int(line[1])
        curdate = date.fromordinal(date.toordinal(date(year-1,12,31)) + dayofyear)
        dropdatetime = datetime.fromisoformat(curdate.isoformat())
        dropdatetime += timedelta(hours=int(line[2][:2]),minutes=int(line[2][2:4]))
        lat = np.double(line[3])
        lon = np.double(line[4])
        num = int(line[5])

        #setting up lists
        temperature = []
        salinity = []
        depth = []

        #reading data
        all_data = []
        for line in f_in:
            all_data.extend([np.double(item) for item in line.strip().split()])
        
            
        #parsing into profiles
        if hasSal:
            for i in range(0,len(all_data),3):
                temperature.append(all_data[i])
                depth.append(all_data[i+1])
                salinity.append(all_data[i+2])
                
        else:
            for i in range(0,len(all_data),2):
                temperature.append(all_data[i])
                depth.append(all_data[i+1])
            
            
    #converting data to arrays and returning
    data = {"depth":np.asarray(depth), "temperature":np.asarray(temperature), "salinity":np.asarray(salinity)}
    return [data,dropdatetime,lat,lon,num]




#write data to FIN file
def writefinfile(finfile,cdtg,lat,lon,num,depth,temperature,salinity=None):
    
    hasSal = False
    if salinity is not None:
        hasSal = True
    
        
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
        
        
    if salinity is not None: #TODO: add ability to write currents as well
        hasSal = True
    else:
        hasSal = False
        
        
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
    
    if hasSal: #TODO: change replicator and append salinity 
        fxy_all[5] = int('0100001100000000',2) #change replicator (index 5) to delayed rep of 3 descriptors, code 1/03/000
        fxy_all.append(int('0001011000111110',2)) #0/22/062: salinity (14b, 100*PPT (=PSU))

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
    bufrarray = bufrarray + format(len(temperature),'16b')  # delayed descripter replication factor(0,31,001) = length
    
    #converting temperature, salinity (as req'd.), and depth and writing
    if hasSal:
        for t,d,s in zip(temperature,depth,salinity):
            d_in = int(np.round(d*10)) # depth (0,07,062)
            bufrarray = bufrarray + format(d_in,'017b')
            t_in = int(np.round(100 * (t + 273.15)))  # temperature (0,22,043)
            bufrarray = bufrarray + format(t_in,'015b')
            s_in = int(np.round(100 * s))  # salinity (0,22,062)
            bufrarray = bufrarray + format(s_in,'014b')
        
        
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

        