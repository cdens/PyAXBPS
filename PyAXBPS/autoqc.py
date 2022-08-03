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

def autoqc(rawdata,rawdepth,smoothlev,profres,maxdev,checkforgaps):
    
    #remove NaNs
    isgood = [False if (np.isnan(rd*rz)) else True for rd,rz in zip(rawdata,rawdepth)]
    rawdepth = rawdepth[isgood]
    rawdata = rawdata[isgood]
    

    #Step 1: Find and remove gaps due to VHF interference kicking off early Mk21 start (should be AXBT issue only)
    if checkforgaps and len(rawdepth) > 0:
        rawdepth,rawdata = removegaps(rawdepth,rawdata)
    
    #Step 2: remove spikes using running standard deviation filter
    if len(rawdepth) > 0:
        depth_despike,data_despike = rundespiker(rawdepth,rawdata,maxdev)
    else:
        depth_despike = data_despike = []
        
    #Step 3: smooth the despiked profile
    if len(depth_despike) > 0:
        depth_smooth,data_smooth = runsmoother(depth_despike,data_despike,smoothlev)
    else:
        depth_smooth = data_smooth = []
    
    #Step 4: pull critical points from profile to save (note- returned variables are lists)
    if len(depth_smooth) > 0:
        depth,data = subsample_profile(depth_smooth,data_smooth,profres)
    else:
        depth = data = []

    #add surface value if one doesn't exist
    if len(depth) > 0:
        if depth[0] != 0:
            sst = data[0]
            depth.insert(0,0)
            data.insert(0,sst)
    
    #convert back to numpy arrays
    data = np.array(data)
    depth = np.array(depth)
    
    return [data,depth]
    
    
    
    
    
#function to identify and remove gaps due to false starts from interference
def removegaps(rawdepth,rawdata):
    donegapcheck = False
    while not donegapcheck:
        maxcheckdepth = 50 #only checks the upper 50m of the profile
        maxgapdiff = 10 #if gap is larger than this range (m), correct profile
        isgap = [0]
        for i in range(1,len(rawdepth)):
            if (rawdepth[i] >= rawdepth[i-1]+ maxgapdiff) and (rawdepth[i-1] <= maxcheckdepth): #if there is a gap of sufficient size to correct
                isgap.append(1) #NOTE: the logical 1 is placed at the first depth AFTER the gap
            else:
                isgap.append(0)
        
        #if there are gaps, find the deepest one and correct t/d profile with that depth as the surface (only works with linear fall rate equation)
        if np.sum(isgap) > 0:
            lastgap = np.max(np.argwhere(isgap))
            realstartdepth = rawdepth[lastgap]
            rawdata = rawdata[lastgap:]
            rawdepth = rawdepth[lastgap:]-realstartdepth
        else: #otherwise, exit loop
            donegapcheck = True
        
    return rawdepth,rawdata
    
    
    
    
#removes spikes from profile with depth-based standard deviation filter
def rundespiker(rawdepth,rawdata,maxdev):
    data_despike = np.array([])
    depth_despike = np.array([])
    
    depthwin = 30 #range of spiker is +/- 5 meters
    maxdepth = np.max(rawdepth)
    
    for n,cdepth in enumerate(rawdepth):
        
        #assigning region for running standard deviation filter
        if cdepth <= depthwin:
            goodindex = np.less_equal(rawdepth,depthwin)
        elif cdepth >= maxdepth - depthwin:
            goodindex = np.greater_equal(rawdepth, maxdepth-depthwin)
        else:
            ge = np.greater_equal(rawdepth, cdepth - depthwin) #all depths above bottom threshold
            le = np.less_equal(rawdepth, cdepth + depthwin) #all depths below top threshold
            goodindex = np.all([ge,le],axis=0) #all depths with both requirements satisfied
            
        #pulling subset
        dataspike = rawdata[goodindex]
        
        #mean and standard deviation of current range
        curmean = np.mean(dataspike)
        curstd = np.std(dataspike)
        
        #only retain values within +/- 1 standard deviation of running mean or top 10 m
        if abs(rawdata[n]-curmean) <= maxdev*curstd or rawdepth[n] < 10:
            depth_despike = np.append(depth_despike,rawdepth[n])
            data_despike = np.append(data_despike,rawdata[n])

    return depth_despike,data_despike    
            
            
    
            
#run depth-based smoother- ensures that first and last datapoints match original profile
def runsmoother(depth_despike,data_despike,smoothlev):
    
    data_smooth = np.array([])
    depth_smooth = depth_despike.copy()
    mindepth = np.min(depth_despike)
    maxdepth = np.max(depth_despike)

    for n,cdepth in enumerate(depth_despike):
        if cdepth == mindepth or cdepth == maxdepth: #if first or last point in profile, append current dataerature
            data_smooth = np.append(data_smooth,data_despike[n])
            
        else: #otherwise- average a range of points
            if cdepth <= smoothlev/2: #in top of profile
                cursmoothlev = 2*cdepth
            elif cdepth >= maxdepth - smoothlev/2: #in bottom of profile
                cursmoothlev = 2*(maxdepth - cdepth)
            else: #in middle of profile
                cursmoothlev = smoothlev
                
            ge = np.greater_equal(depth_despike, cdepth - cursmoothlev/2)
            le = np.less_equal(depth_despike,cdepth + cursmoothlev/2)            
            goodindex = np.all([ge,le],axis=0)
                
            #append mean of all points in range as next datapoint
            data_smooth = np.append(data_smooth,np.mean(data_despike[goodindex]))
            
    return depth_smooth,data_smooth
    

    
    
#subsample profile
def subsample_profile(depth_smooth,data_smooth,profres):
    
    dtdz = [] #calculating profile slope (append 0 at start and end so length matches that of depth_smooth)
    dtdz.append(0)
    for i in range(1,len(data_smooth)-1): #need to use range here because we are only interested in a subset of indices
        #dtdz = (t3 - t1)/(z3 - z1): centered on z2
        dtdz.append(((data_smooth[i+1] - data_smooth[i-1])/ 
              (depth_smooth[i+1]-depth_smooth[i-1])))
    dtdz.append(0)

    depth = []
    dataerature = []
    lastdepth = -100000 #large enough value so subsampler will always grab first datapoint
    for i,cdepth in enumerate(depth_smooth):
        #constraint on first derivative of data with depth (no unrealistic spikes), if the point is a critical value given the selected resolution level, append to output profile
        if dtdz[i] <= 0.5 and cdepth-lastdepth >= profres:
            depth.append(cdepth)
            dataerature.append(data_smooth[i])
            lastdepth = cdepth
    
    return depth,dataerature