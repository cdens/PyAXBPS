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
            

            
            
            
            


            

        