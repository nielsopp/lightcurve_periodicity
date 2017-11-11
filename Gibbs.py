# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 10:07:06 2014

@author: niels
"""

import sys
from nifty import *
import os

#name of directory to store results in:
directory = sys.argv[1]

#name of file containing data:
data_file = sys.argv[2]

#number of bins for the resonstruciton:
Npix = int(sys.argv[3])

#factor by which reconstruction will be extended beyond observational baseline:
extension_factor = float(sys.argv[4])

#length of the burn-in phase:
burn_in = int(sys.argv[5])

#correlation length:
corr_len = int(sys.argv[6])

#logarithmic difference between the means of the two chains defining convergence:
conv_crit = float(sys.argv[7])


#setting random seed
np.random.seed(102)

#create directory for output:
if directory[-1] != '/':
    directory += '/'
if not os.path.exists(directory):
    os.makedirs(directory)

#read data:
d = np.genfromtxt(data_file)
time = d[:,0]
mag = d[:,1]
meanmag = mag.mean()
sigma = d[:,2]

#set up binning and nifty spaces:
maxperiod = extension_factor*(time.max() - time.min())
res = (maxperiod)/Npix
timebins = np.arange(Npix)*res + time.min() - (extension_factor - 1.)/2.*(time.max() - time.min())
assign = np.digitize(time,bins=timebins)

#position space
px = rg_space(num=Npix,zerocenter=True,purelyreal=True,dist=res,fourier=False)
#auxiliary position space
px_aux = rg_space(num=Npix,zerocenter=True,purelyreal=False,dist=res,fourier=False,hermitian=False)
#harmonic analog
lm = px.get_codomain()
#auxiliary harmonic space
lm_aux = px_aux.get_codomain()
#data space
ds = point_space(num=len(time))

#write out discretized time axis
np.save(directory + 'timebins.npy',timebins)


def WF(Dinv,R,N,d):
    """
        Calculate the Wiener filter estimate.
    """
    out = N.inverse_times(d)
    out = R.adjoint_times(out)
    out = Dinv.inverse_times(out)
    return out
    
    
def draw_signal(S,R,N,Dinv,data):
    """
        Draw lightcurve from posterior with fixed power spectrum.
    """
    stilde = S.get_random_field(domain=px)
    ntilde = N.get_random_field(domain=ds)
    dtilde = R(stilde) + ntilde
    d = data - dtilde
    spr = WF(Dinv,R,N,d)
    return stilde + spr
    

def draw_power(s,rho):
    """
        Draw power spectrum from posterior with fixed lightcurve.
    """
    totsigpow = s.power(bare=True)*rho
    q = np.random.chisquare(rho)
    #Threshold power spectrum to avoid getting too close to zero:
    power = np.maximum(totsigpow/q, 1.e-5)
    return power

    
def Gibbs(data,N,R,pindex,pundex,kindex,rho,burn_in,corr_len,conv_crit):
    """
        Run two Gibbs sampling chains.
    """
    M_aux = diagonal_operator(px_aux,diag=(R.adjoint_times(N.inverse_times(R(1)))).val,
                          bare=False)
    
    #set up explicit M-matrix in Fourier space:
    M_mat = np.zeros((lm.dim(),lm.dim()),dtype=complex)
    for jj in range(lm.dim()):
        xi = field(lm_aux,val=0.)
        xi[jj] = 1.
        M_mat[:,jj] = (M_aux(xi)).val/lm_aux.vol
    
    pow1 = np.random.normal(0.,5.,len(kindex))
    pow2 = np.random.normal(0.,5.,len(kindex))

    #save logarithmic power spectra
    np.save(directory + 'power1.npy',pow1)
    np.save(directory + 'power2.npy',pow2)
    power1 = np.exp(pow1)
    power2 = np.exp(pow2)

    Sinv1_mat = np.diag(1./power1[pindex]/lm.vol**2)
    Sinv2_mat = np.diag(1./power2[pindex]/lm.vol**2)
    S1 = power_operator(lm,spec=power1,bare=True,pindex=pindex)
    S2 = power_operator(lm,spec=power2,bare=True,pindex=pindex)
    
    Dinv_matrix1 = Sinv1_mat + M_mat
    Dinv_matrix2 = Sinv2_mat + M_mat
    Dinv1 = explicit_operator(lm,matrix=Dinv_matrix1,bare=True,sym=True)
    Dinv2 = explicit_operator(lm,matrix=Dinv_matrix2,bare=True,sym=True)
    m1 = WF(Dinv1,R,N,data)
    m2 = WF(Dinv2,R,N,data)
    count = 0
    kk = 0
    while np.max(np.abs(pow1 - pow2)) > conv_crit:
        m1 = draw_signal(S1,R,N,Dinv1,data)
        m2 = draw_signal(S2,R,N,Dinv2,data)
        power1 = draw_power(m1,rho)
        power2 = draw_power(m2,rho)
        Sinv1_mat = np.diag(1./power1[pindex]/lm.vol**2)
        Sinv2_mat = np.diag(1./power2[pindex]/lm.vol**2)
        S1.set_power(power1,bare=True,pindex=pindex)
        S2.set_power(power2,bare=True,pindex=pindex)
        Dinv_matrix1 = Sinv1_mat + M_mat
        Dinv_matrix2 = Sinv2_mat + M_mat
        Dinv1 = explicit_operator(lm,matrix=Dinv_matrix1,bare=True,sym=True)
        Dinv2 = explicit_operator(lm,matrix=Dinv_matrix2,bare=True,sym=True)
        if count > burn_in - 1:
            if count % corr_len == corr_len - 1:
                print(count, np.max(np.abs(pow1 - pow2)))
                pow1 *= kk
                pow2 *= kk
                pow1 += np.log(power1)
                pow2 += np.log(power2)
                pow1 /= (kk + 1.)
                pow2 /= (kk + 1.)
                np.save(directory + 'power1.npy',pow1)
                np.save(directory + 'power1_%05i.npy'%kk,power1)
                np.save(directory + 'power2.npy',pow2)
                np.save(directory + 'power2_%05i.npy'%kk,power2)
                np.save(directory + 'm1_%05i.npy'%kk,m1.val)
                np.save(directory + 'm2_%05i.npy'%kk,m2.val)
                kk += 1
        count += 1
    np.save(directory + 'power1.npy',pow1)
    np.save(directory + 'power2.npy',pow2)
    return
    

def run():
    """
        Running the whole shebang.
    """
    kindex, rho, pindex, pundex = lm.get_power_indices()

    #saving frequency grid:
    np.save(directory + 'freqs.npy',kindex)

    #set up noise-covariance and response operators:
    N = diagonal_operator(ds,sigma**2)
    R = response_operator(domain=px,target=ds,den=False,assign=assign)
            
    Gibbs(mag - meanmag,N,R,pindex,pundex,kindex,rho,burn_in,corr_len,conv_crit)

    return 0


if __name__ == '__main__':
    run()
