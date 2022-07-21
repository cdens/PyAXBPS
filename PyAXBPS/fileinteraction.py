# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    This file is part of the Airborne eXpendable Buoy Processing System
#
#    ARES is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ARES is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ARES.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

import numpy as np
from datetime import date, datetime, timedelta
import chardet






# =============================================================================
#                   Date/Time, Lat/Lon Parsing
# =============================================================================
            
    
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
    
    
    
    
    
    
    
    

# =============================================================================
#                   Log (DTA) file read/write
# =============================================================================


#read raw temperature/depth profile from LOGXX.DTA file
#ignoring date/time as it isn't trusted over user input (bad computer time settings on aircraft kits)
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
    
    
    
    
    

#writing data to log file
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
            
            
    
                
                
                
                
                
                
            
                

# =============================================================================
#                   EDF file reading/writing
# =============================================================================
    
def readedffile(edffile):
    
    encoding = 'utf-8'
    
    lon = lat = day = month = year = hour = minute = second = False #variables will be returned as 0 if unsuccessfully parsed
    
    data = {'depth':[],'temperature':[],'salinity':[],'U':[],'V':[]} #initializing each data field as None so if it isn't none upon function completion the user knows the file had a match for that field
    
    fields_match = [['depth'], ['temperature','temp'], ['salinity','sal'], 
        ['zonal','east/west'], #U current
                ['meridional','north/south']] #V current
                
    fields = ['depth', 'temperature', 'salinity','U','V'] #V current
    nfields = 2 #initial number of fields- temperature and depth only
    fieldcolumns = [0,1,-1,-1,-1]
    
    fields_init = False
    
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
                        curlinefield = line[1].lower().strip()
                        curfieldcolumn = int(line[0].strip()[-1]) - 1
                        for i,fieldops in enumerate(fields_match):
                            for field in fieldops:
                                # print(f"{field=}, {curlinefield=}, ")
                                if field in curlinefield:# and fieldcolumns[i] > -1:
                                    fieldcolumns[i] = curfieldcolumn
                                    matched = True
                                    break #don't need to keep checking
                        if not matched:
                            curlinefield = line[1].strip()
                            fields.append(curlinefield)
                            fieldcolumns.append(curfieldcolumn)
                            data[curlinefield] = []
                            
                        #recalculate each time after adding a new field
                        nfields = np.sum([1 if cc >= 0 else 0 for cc in fieldcolumns]) #total number of fields in file
                        corresponding_field_column = [-1 for _ in range(nfields)]
                        for i,ind in enumerate(fieldcolumns):
                            if ind > -1:
                                corresponding_field_column[ind] = i
                        
                        fields_init = True
                        
                #space or tab delimited data (if number of datapoints matches fields and line doesnt start with //)
                elif cline[:2] != '//' and fields_init: 
                    curdata = cline.strip().split()
                    if len(curdata) == nfields:
                        for i,cdata in enumerate(curdata):
                            cfield = fields[corresponding_field_column[i]]
                            data[cfield].append(cdata)
                    
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
    
    for cfield in ['depth','temperature','salinity','U','V']:  
        data[cfield] = np.asarray(data[cfield])   
        
    return data,dropdatetime,lat,lon

    
    
    

def writeedffile(edffile,dropdatetime,lat,lon,data,comments,field_formats=[],QC=False):
    
    
    with open(edffile,'w') as f_out:
        
        
        #identifying fields to write
        fields = list(data.keys())
        
        #determining what string format to use for each field if default no specification is given
        if field_formats == []:
            for cfield in fields:
                cfieldlow = cfield.lower()
                if 'temperature' in cfieldlow or 'conductivity' in cfieldlow or 'salinity' in cfieldlow:
                    field_formats.append('7.2f') #__XX.XX
                elif 'current' in cfieldlow:
                    field_formats.append('8.3f') #__X.Xx
                elif 'depth' in cfieldlow or 'time' in cfieldlow or 'frequency' in cfieldlow:
                    field_formats.append('9.2f') #__XXXX.XX
                else:
                    field_formats.append('') #unspecified format
                
        npts = len(data[fields[0]])
        
        #identifying which fields to write to the file
        field_inds_write = [i for i,ckey in enumerate(data) if len(data[ckey]) > 0 ]
        fields_write = [fields[i] for i in field_inds_write]
        
        field_formats_write = [field_formats[i] for i in field_inds_write]
    
        #writing header, date and time, drop # (bad value)
        f_out.write("// This is an air-launched probe EXPORT DATA FILE  (EDF)\n")
        f_out.write("// File generated with the AXBT Realtime Editing System (ARES)\n")
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
""" + "\n".join([f"Field{i+1}  :  {ckey}" for i,ckey in enumerate(fields_write)]) + """
//""" + qcstr + """
""" + "\t".join([ckey for ckey in fields_write]) + "\n"
        
        f_out.write(drop_settings_info)
        
        #writing data
        for i in range(npts):
            cline = "\t".join([f"{data[cfield][i]:{cformat}}" for cfield,cformat in zip(fields_write,field_formats_write)]) #tab-delimited
            f_out.write(cline + "\n")

    
            
            
            
            
            
            
            
            
            

                

# =============================================================================
#                       DAT (NOAA) file reading/writing
# =============================================================================            
            
            
def readdatfile(datfile):

    data = {'temperature':[], 'salinity':[], 'depth':[], 'U':[], 'V':[]} #initializing data lists
    header_read = False

    with open(datfile) as f_in:

        for i,cline in enumerate(f_in): 
            cline = [i for i in cline.strip().split(' ') if len(i) > 0] #space delimited
            
            if not header_read:
                try:#header
                    dtgstr = cline[0]+cline[1]
                    if len(dtgstr) == 14:
                        dtgformat = '%Y%m%d%H%M%S'
                    else:
                        dtgformat = '%Y%m%d%H%M' #correct to best guess of format
                        dtgstr = dtgstr[:12]
    
                    dropdatetime = datetime.strptime(dtgstr,dtgformat)
                    lat = float(cline[2])
                    lon = float(cline[3])
                    
                    header_read = True
                    
                except: #not the header- format is wrong
                    pass
                    
            else: #header is already read
                if len(cline) > 1:
                    data['depth'].append(float(cline[0]))
                    data['temperature'].append(float(cline[1]))
                    
                    hasSal = False
                    hasCurrent = False
                    if len(cline) >= 3: #has salinity 
                        hasSal = True
                    if len(cline) >= 5: #has currents
                        hasCurrent = True
                        
                    if hasSal:
                        cpsal = float(cline[2])
                        if cpsal >= 0:
                            data['salinity'].append(cpsal)
                        else:
                            data['salinity'].append(np.NaN)
                    if hasCurrent:
                        u = float(cline[3])
                        if u > -999:
                            data['U'].append(u)
                        else:
                            data['U'].append(np.NaN)
                        v = float(cline[4])
                        if v > -999:
                            data['V'].append(cpsal)
                        else:
                            data['V'].append(np.NaN)
                            

    for cfield in ['depth','temperature','salinity']:  
        data[cfield] = np.asarray(data[cfield])   

    return dropdatetime,lat,lon,data




#TODO: understand format of 1st 4 lines -> what is second time, what is 12, what are first two lines for?
# ****0000005762****
# SOFX01 KWBC 172308
#
# 20210817 225634  17.230  -79.962 N43RF AL072021 12 2021-229-23:08:29

def writedatfile(datfile, dropdatetime, lat, lon, headerstart, tailnum, missionnum, depth, temperature, salinity=None, U=None, V=None):

    with open(datfile,'w') as f_out:

        header = headerstart + f"\n{dropdatetime:%Y%m%d} {dropdatetime:%H%M%S} {lat:7.3f} {lon:8.3f} {tailnum} {missionnum} 12 {datetime.utcnow():%Y-%j-%H:%M:%S}\n"

        f_out.write(header)

        if salinity is None:
            salinity = [np.NaN]*len(depth)
        if len(salinity) != len(depth): #now salinity is definitely a list, can call len on it
            salinity = [np.NaN]*len(depth)
            
        if U is not None and V is not None:
            hasCurrent = True
        else:
            hasCurrent = False
            U = [np.NaN]*len(depth)
            V = [np.NaN]*len(depth)
        
        #replace any NaNs with negative value
        # temperature = [-999 if np.isnan(i) else i for i in temperature]
        salinity = [-9.99 if np.isnan(i) else i for i in salinity] 
        U = [-999 if np.isnan(i) else i for i in U]
        V = [-999 if np.isnan(i) else i for i in V]

        for (cd,cT,cS,cU,cV) in zip(depth,temperature,salinity,U,V):
            if not np.isnan(cd*cT):
                if hasCurrent:
                    f_out.write(f'\n{cd:6.1f}{cT:6.2f}{cS:6.2f}{cU:8.3f}{cV:8.3f}\n')
                else:
                    f_out.write(f'\n{cd:6.1f}{cT:6.2f}{cS:6.2f}\n')

        f_out.write('\n\n') #ends with two more empty lines
            
            
        
            
            
            
            
            

                

# =============================================================================
#                   FIN/NVO/TXT file reading/writing
# =============================================================================

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
        num_obs = int(line[7])
        
        #reading data
        all_data = []
        for line in f_in:
            all_data.extend([np.double(item) for item in line.strip().split()])
        
    #file closes here, figure out what data is included now
    #total types of obs (e.g. temp/depth) = # of datapoints in file divided by # of obs reported in header 
    Npts = len(all_data) #number of points in file
    nobtypes = Npts/num_obs
    
    #raise an error if the number of reported profile points in the header isn't evenly divisble by 
    #the total amount of data in the file. This is critical because we need that resulting quotient (# ob types)
    #to determine what kind of data is included (2=T/z, 3=T/z/S, 4=T/z/u/v, 5=T/z/S/u/v)
    try:
        assert(nobtypes == np.round(nobtypes))
        nobtypes = int(nobtypes) #as long as they match, switch it to an integer
    except AssertionError as e:
        raise Exception('Number of reported profile points should be evenly divisible by number of datapoints in file') from e
    
    #setting up lists
    temperature = []
    salinity = []
    depth = []
    U = []
    V = []
    
    #determining whether or not salinity and currents are included from # obs per depth value
    #2=T/z, 3=T/z/S, 4=T/z/u/v, 5=T/z/S/u/v
    hasSal = False
    hasCurrent = False
    if nobtypes >= 4:
        hasCurrent = True
    if nobtypes == 3 or nobtypes == 5:
        hasSal = True
        
    #parsing into profiles
    for i in range(0,len(all_data)-1,nobtypes):
        temperature.append(all_data[i]) #1=T
        depth.append(all_data[i+1]) #2=z
        if hasSal: #has salinity- 3=S
            salinity.append(all_data[i+2])
        elif hasCurrent: #no salinity so 1=T,2=z,3=u,4=v
            U.append(all_data[i+2])
            V.append(all_data[i+3])
        if hasSal and hasCurrent: # salinity and currents, so 1=T,2=z,3=S,4=u,5=v
            U.append(all_data[i+3])
            V.append(all_data[i+4])
            
            
    #converting data to numpy arrays, storing in dict and returning
    data = {"depth":np.asarray(depth), "temperature":np.asarray(temperature), "salinity":np.asarray(salinity), "U":np.asarray(U), "V":np.asarray(V)}
    
    return [data,dropdatetime,lat,lon,num]




#write data to FIN file
def writefinfile(finfile,cdtg,lat,lon,num,depth,temperature,salinity=None,U=None,V=None):
    
    #determining whether to write salinity/current data, adjusting S/U/V values if not to function with code
    hasSal = False
    hasCurrent = False
    if salinity is not None:
        hasSal = True
    else:
        salinity = np.array([None] * len(depth))
    if U is not None and V is not None:
        hasCurrent = True
    else:
        U = np.array([None] * len(depth))
        V = np.array([None] * len(depth))
    
    #getting rid of any NaNs
    isgood = np.array([False if np.isnan(cT*cz) else True for (cz,cT) in zip(depth, temperature)])
    depth = depth[isgood]
    temperature = temperature[isgood]
    salinity = salinity[isgood]
    U = U[isgood]
    V = V[isgood]
    
    #setting number of observation vars (e.g. T,z,S) and datapoints per line 
    if hasSal and hasCurrent: #temperature, salinity, U, V vs. depth
        obsperline = 2
        nobtypes = 5
    elif hasCurrent: #temperature, U, V vs. depth (AXCP)
        obsperline = 3
        nobtypes = 4
    elif hasSal: #temperature and salinity vs. depth (AXCTD)
        obsperline = 3
        nobtypes = 3
    else: #temperature vs. depth only (AXBT)
        obsperline = 5
        nobtypes = 2
    
    #setting up data_out array
    #creating list of data point
    Npts = len(depth)
    outdata = []
    for (t,d,s,u,v) in zip(temperature,depth,salinity,U,V):
        outdata.extend([f"{t: 8.3f}",f"{d: 8.1f}"]) #Temp (C) to 3 dec., depth (m) to 1 dec.
        if s is not None:
            outdata.append(f"{s: 8.3f}") #salinity (PSU) to 3 dec
        if u is not None and v is not None:
            outdata.extend([f"{u: 7.3f}",f"{v: 7.3f}"]) #U/V (m/s) to 3 dec
    
    
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
        f_out.write(f"{cdtg.year}   {dayofyear:03d}   {datetime.strftime(cdtg,'%H%M')}   {latsign}{lat:06.3f}   {lonsign}{lon:07.3f}   {num:02d}   6   {Npts}   0   0   \n")
                
        pointsperline = nobtypes*obsperline

        #writing profile data
        i = 0
        Ndata = nobtypes * Npts
        while i < Ndata:
            if i+pointsperline < Ndata:
                pointstopull = pointsperline
            else:
                pointstopull = Ndata - i
                
            line = "".join(outdata[i:i+pointstopull]) + "\n"
            f_out.write(line)
            i += pointstopull
            

            
        
        
        
            
                

# =============================================================================
#                   JJVV file reading/writing
# =============================================================================            
    
    
#read data from JJVV file (AXBT and/or temperature data only)
#relyear is a year in the decade of the drop, otherwise the code assumes the drop was within the last 
#10 years and not in the future by more than a day (to correct for timezone issues)
def readjjvvfile(jjvvfile,relyear=None):
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
        hour = int(line[2][:2])
        minute = int(line[2][2:4])
        
        
        curdatetime = datetime.utcnow()        
        if relyear is not None: #determining year (drop must have been in same decade as relyear)
            decade = int(np.floor(relyear/10)*10)
            year = decade + yeardigit
            
        else: #determining year (drop must have been within last 10 years)
            curyear = curdatetime.year #assuming drop was in current decade
            year = int(np.floor(curyear/10)*10) + yeardigit #getting year of drop
        
        #building guess of drop date/time
        dropdatetime = datetime(year,month,day,hour,minute,0)
            
        #making sure drop date/time are within current decade, but not more than 1 month after 
        #current computer date/time (1 month is leeway for date/time issues), if your computer is 
        #more than a month off the actual date/time, that's on you friend
        dt_days = (1/(24 * 3600)) * (curdatetime - dropdatetime).total_seconds() #negative = drop is after current time
        year_correction = 0
        if dt_days < -30: #drop is more than 30 days in the future
            year_correction = -10 #send the drop back a decade
        elif relyear is None and dt_days >= 365.25 * 10: #more than a decade in the past without relyear to specify
            year_correcion = 10
            
        #do this rather than add timedelta because the slight rounding issue for # of seconds in a year will change the drop time otherwise even using days=3652.5 for a decade
        old = dropdatetime
        dropdatetime = datetime(old.year + year_correction, old.month, old.day, old.hour, old.minute, old.second)

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
    
    #getting rid of any NaNs
    isgood = np.array([False if np.isnan(cT*cz) else True for (cz,cT) in zip(depth, temperature)])
    depth = depth[isgood]
    temperature = temperature[isgood]
    
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
            
            
            

            

            
# =============================================================================
#                   BUFR File writing only
# =============================================================================            

        

def writebufrfile(bufrfile, cdtg, lon, lat, identifier, originatingcenter, depth, temperature, salinity=None, U=None, V=None, optionalinfo=None):
    
    
    hasSal = False
    hasCurrent = False
    
    #removing NaNs from profile
    nancheck = depth * temperature
    if hasSal:
        nancheck *= salinity
    if hasCurrent:
        nancheck *= U * V
        
    isgood = np.array([False if np.isnan(cnancheck) else True for cnancheck in nancheck])
    depth = depth[isgood]
    temperature = temperature[isgood]
    
    
    binarytype = 'big'  # big-endian
    reserved = 0
    version = 4  # BUFR version number (3 or 4 supported)

    # section 1 info
    if version == 3:
        sxn1len = 18
    elif version == 4:
        sxn1len = 22
        
    if optionalinfo is None:
        hasoptionalsection = False
    else:
        hasoptionalsection = True

    mastertable = 0  # For standard WMO FM 94-IX BUFR table
    originatingsubcenter = 0
    updatesequencenum = 0  # first iteration
    if hasoptionalsection:
        hasoptionalsectionnum = int('10000000', 2)
    else:
        hasoptionalsectionnum = int('00000000', 2)
    datacategory = 31  # oceanographic data (Table A)
    datasubcategory = 3 #bathy (JJVV) message
    versionofmaster = 32
    versionoflocal = 0
    yearofcentury = int(cdtg.year - 100 * np.floor(cdtg.year / 100))  # year of the current century

    # section 2 info
    if hasoptionalsection:
        sxn2len = len(optionalinfo) + 4
    else:
        sxn2len = 0

    # Section 3 info
    sxn3len = 25  # 3 length + 1 reserved + 2 numsubsets + 1 flag + 2 FXY = 9 octets
    numdatasubsets = 1
    
    # whether data is observed, compressed (bits 1/2), bits 3-8 reserved (=0)
    s3oct7 = int('10000000', 2)
    
    
    
    fxy_all = [int('0000000100001011', 2), #0/01/011: identifier (72b, 9 bytes ascii)
    int('1100000100001011', 2), #3/01/011: year/month/day (12b/4b/6b = 22b total)
    int('1100000100001100', 2), #3/01/012: hour/minute (5b/6b = 11b total)
    int('1100000100010111', 2), #3/01/023: lat/lon coarse accuracy (15b lat/16b lon = 31b total)
    int('0000001000100000', 2), #0/02/032: digitization indicator (2b, 00)
    int('0100001000000000', 2), #1/02/000: replication of 2 descriptors
    int('0001111100000010', 2), #0/31/002: extended delayed decriptor replication factor (16b, # obs)
    int('0000011100111110', 2), #0/07/062: depth (17b, m*10)
    int('0001011000101010', 2),] #0/22/043: temperature (15b, degK*100)    
    
    #replication factor explanation:
    #F=1 for replication, X=number of descriptors (2 for d/T, 3 for d/T/S), Y=#times replicated (1) 
        
    
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
    
    
    #converting temperature and depth and writing
    for t,d in zip(temperature,depth):
        d_in = int(np.round(d*10)) # depth (0,07,062)
        bufrarray = bufrarray + format(d_in,'017b')
        t_in = int(np.round(10 * (t + 273.15)))  # temperature (0,22,042)
        bufrarray = bufrarray + format(t_in,'012b')
    

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
        if hasoptionalsection:
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
        
        
        
        
        
        