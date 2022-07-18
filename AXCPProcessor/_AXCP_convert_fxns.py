# =============================================================================
#     Author: Casey R. Densmore, 12FEB2022
#
#    AXBPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AXBPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AXBPS.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================

# =============================================================================
#   This code is translated/updated from the MATLAB xcpdsp.m script written by 
#       John Dunlap, University of Washington APL, 4 November 2009
# =============================================================================

import numpy as np


#conversion: coefficients=C,  D_out = C[0] + C[1]*D_in + C[2]*D_in^2 + C[3]*D_in^3 + ...
def dataconvert(data_in,coefficients):
    
    datatype = 1 #integer or float
    if type(data_in) == list:
        datatype = 2
    elif type(data_in) == np.ndarray: #numpy array
        dataype = 3
        
    if datatype == 1:
        data_in = [data_in]
        
    output = []
    for cur_data_in in data_in:
        cur_output = 0
        for (i,c) in enumerate(coefficients):
            cur_output += c*cur_data_in**i
        output.append(cur_output)
        
    if datatype == 1: #convert back from list to int/float
        output = output[0]
    elif datatype == 3: #convert to np array
        output = np.asarray(output)
            
    return output
    
    


def calc_temp_from_freq(self, freq, depth):
    teres = 1.0/(4.4*freq*self.tcalcap)-self.tcalrs
    if teres>0:
        ln = np.log(teres) 
    else:
        ln = np.NaN 

    temp=1.0/(self.tcal[0]+ln*(self.tcal[1]+ln*(self.tcal[2]+ ln*self.tcal[3]))) - 273.15
    
    # adjust temperature to match mendo.rt processing
    temp = temp + self.tcor[0] + self.tcor[1] * depth
    
    return temp
    



def calc_vel_components(self, fcss, fess, tss):

    # tz to get rotation period
    # interpolate fcss to get positive zero crossing times, tz
    j = np.where(np.isfinite(fcss))[0]
    x = fcss - np.mean(fcss[j]) 
    r = np.ones(len(x))
    r[x < 0] = -1
    jp = np.where(np.diff(r) > 0)[0]
    dxjp =   x[jp+1] -   x[jp]
    dtjp = tss[jp+1] - tss[jp]
    tzp = tss[jp] - x[jp] * dtjp / dxjp
    per = np.diff(tzp)
    
    if np.sum(np.isfinite(per)) > 1:
        peravg = np.nanmean(per)
        perrms = np.nanstd(per)
    else:
        peravg = np.NaN
        perrms = np.NaN
    
    
    rotfavg = (1/peravg)
    rotfrms = rotfavg * perrms / peravg
    
    # make phase
    phase = np.NaN * np.ones(len(tss))
    for jper,_ in enumerate(per):
        j = np.where((tzp[jper]<tss) & (tss <= tzp[jper+1]))[0]
        phase[j] = 2*np.pi*(tss[j] - tzp[jper]) / per[jper]
    
    # sinusoidal fitting
    j=np.where(np.isfinite(phase))[0]
    nfit = len(j)
    
    #using matrix math to combine individual compass coil (direction) and EF (current speed) measurements for current
    #datapoint and calculate frequency amplitudes and phases for velocity estimates
    if nfit > len(phase)/2:
        aprxcc = np.stack([np.cos(phase[j]), np.sin(phase[j]), np.ones(len(j))])
        coefcc, _,_,_ = np.linalg.lstsq(aprxcc.T, fcss[j]) #coefcc = aprxcc \ fcss(j)
        rescc = fcss[j] - np.matmul(aprxcc.T, coefcc) # rescc = fcss[j] - aprxcc * coefcc
        fccr = np.nanstd(rescc)
        fcca = np.sqrt(coefcc[0]**2+coefcc[1]**2)
        fccp = np.arctan2(-coefcc[1],coefcc[0])
        ccbl = coefcc[2]
        
        aprxef = np.append(aprxcc, np.array([np.linspace(-1,1,len(j))]), axis=0) #aprxef = [aprxcc,  linspace(-1,1,length(j))']
        coefef, _,_,_ = np.linalg.lstsq(aprxef.T, fess[j]) #coefef = aprxef \ fess(j)'
        resef = fess[j] - np.matmul(aprxef.T, coefef) # resef = fess[j] - aprxef * coefef
        fefr = np.nanstd(resef)
        fefa = np.sqrt(coefef[0]**2+coefef[1]**2)
        fefp = np.arctan2(-coefef[1],coefef[0])
        efbl = coefef[2]
        
    else:
        ccbl = np.NaN
        efbl = np.NaN
        fccr = np.NaN
        fefr = np.NaN
        fcca = np.NaN
        fefa = np.NaN
        fccp = np.NaN
        fefp = np.NaN
        
    
    # probe gain and phase angles as function of rotation frequency
    gcca  = dataconvert(rotfavg,self.gcca_poly)
    gccp  = dataconvert(rotfavg,self.gccp_poly)
    gcora = dataconvert(rotfavg,self.gcora_poly)
    gcorp = dataconvert(rotfavg,self.gcorp_poly)
    gefa  = dataconvert(rotfavg,self.gefa_poly)
    gefp  = dataconvert(rotfavg,self.gefp_poly)
    
    # convert frequency amp and phase to velocity estimates
    vc0a = fcca / self.gcvfa / gcca * 1e6
    vc0p = fccp - self.gcvfp - gccp
    
    ve0a1 = fefa / self.gevfa / gefa * 1e6
    ve0p1 = fefp - self.gevfp - gefp
    
    ve0q1 = ve0a1 * np.cos(ve0p1)
    ve0i1 = ve0a1 * np.sin(ve0p1)
    
    ve0a2 = fcca / self.gcvfa / gcca * gcora / gefa * 1e6
    ve0p2 = fccp - self.gcvfp - gccp + gcorp - gefp
    
    ve0q2 = ve0a2 * np.cos(ve0p2)
    ve0i2 = ve0a2 * np.sin(ve0p2)
    
    ve0q = ve0q1 - ve0q2
    ve0i = ve0i1 - ve0i2
    
    ve0a = np.sqrt(ve0q**2+ve0i**2)
    ve0p = np.arctan2(ve0i,ve0q)
    
    if self.revcoil:
        ve0p += np.pi
        
    return rotfavg, rotfrms, efbl, ccbl, fefr, fccr, vc0a, vc0p, ve0a, ve0p, gcca, gefa
    
    
    
    
def calc_currents(self, rotfavg, fccr, fefr, vc0a, vc0p, ve0a, ve0p, gcca, gefa, nindep, w):

    #calculate current U/V components, velocity error, coil area/error
    area = vc0a / rotfavg * self.sfa
    aerr = np.abs(fccr/self.gcvfa/gcca*1e6/rotfavg*self.sfa)/np.sqrt(nindep)
    
    psi = -ve0p + vc0p + np.pi/2 + np.pi
    uw = w * (area/self.amean_rough) * self.sfw
    
    umag =  ve0a * np.cos(psi) * self.sfv
    vmag = -ve0a * np.sin(psi) * self.sfv + uw
    
    verr = np.abs(fefr * self.sfv / self.gevfa / gefa * 1e6) / np.sqrt(nindep)
    
    return area, aerr, umag, vmag, verr
    
    
    
    
    
    
    
    
    
    
    