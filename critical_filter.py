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

#logarithmic difference between two consecutive iterations defining convergence:
conv_crit = float(sys.argv[5])

#strength of the spectral smoothness prior:
smooth = float(sys.argv[6])


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
    

def reconstruct(d,N,R,pindex,pundex,kindex,rho,conv_crit,smoothness):
    """
        Reconstruct the lightcurve and power spectrum using the critical filter.
    """

    M_aux = diagonal_operator(px_aux,diag=(R.adjoint_times(N.inverse_times(R(1)))).val,
                          bare=False)
    
    #set up explicit M-matrix in Fourier space
    M_mat = np.zeros((lm.dim(),lm.dim()),dtype=complex)
    for jj in range(lm.dim()):
        xi = field(lm_aux,val=0.)
        xi[jj] = 1.
        M_mat[:,jj] = (M_aux(xi)).val/lm_aux.vol
    
    pow1 = np.ones(kindex.shape)*8.

    #save logarithmic power spectrum
    np.save(directory + 'power.npy',pow1)
    powspec = np.exp(pow1)

    Sinv_mat = np.diag(1./powspec[pindex]/lm.vol**2)
    
    Dinv_matrix = Sinv_mat + M_mat
    Dinv = explicit_operator(lm,matrix=Dinv_matrix,bare=True,sym=True)
    m = WF(Dinv,R,N,d)
    
    repeat = True

    count = 0
    while repeat:
        powspecold = powspec.copy()
        D_mat = np.linalg.inv(Dinv_matrix)/lm.vol**2
        D = explicit_operator(lm,matrix=D_mat,bare=True,sym=True)
        powspec = infer_power(m,domain=lm,D=D,pindex=pindex,
                              pundex=pundex,kindex=kindex,rho=rho,q=0,
                              alpha=1,smoothness=True,var=smoothness,force=True,
                              bare=True)
        
        #replace monopole with dipole value
        powspec[0] = powspec[1]
        
        print count, np.max(np.abs(np.log10(powspec) - np.log10(powspecold)))
        if np.max(np.abs(np.log10(powspec) - np.log10(powspecold))) < conv_crit:
            repeat = False
        Sinv_mat = np.diag(1./powspec[pindex]/lm.vol**2)
        Dinv_matrix = Sinv_mat + M_mat
        Dinv = explicit_operator(lm,matrix=Dinv_matrix,bare=True,sym=True)
        m = WF(Dinv,R,N,d)
        count += 1
        if count == 500:
            repeat = False
            print('No convergence in 500 iterations')
        np.save(directory + 'power.npy',powspec)
        np.save(directory + 'm.npy',m.val)

    #calculate local uncertainty for lightcurve:    
    Dhat = np.zeros(px.dim())
    for ii in range(Npix):
        xi = field(px,val=0.)
        xi[ii] = 1.
        Dhat[ii] = D(xi)[ii]/px.vol
    np.save(directory + 'Dhat.npy',Dhat)

    #calculate local uncertainty estimate for power spectrum:    
    m = m.transform(target=lm)
    Hess = np.zeros((len(kindex),len(kindex)))
    L,I = nifty_power._calc_laplace(kindex)
    Amem = np.dot(L.T,np.dot(np.diag(I,k=0),L,out=None),out=None)
    T2 = 2/1.e5*Amem
    Tp = np.dot(T2,np.log(powspec))
    for ii in range(len(kindex)):
        for jj in range(len(kindex)):
            for q in np.arange(lm.dim())[pindex == ii]:
                for k in np.arange(lm.dim())[pindex == jj]:
                    mq = m.val[q]
                    mk = m.val[k]
                    Dqk = D_mat[q,k]
                    Hess[ii,jj] += 2.*np.real(mq*np.conjugate(mk)*Dqk) + np.abs(Dqk)**2
            Hess[ii,jj] /= 2.*powspec[ii]*powspec[jj]
            Hess[ii,jj] += T2[ii,jj]
        Hess[ii,ii] += 0.5*rho[ii] + Tp[ii]
    invHess = np.linalg.inv(Hess)
    invHessdiag = np.diagonal(invHess)
    np.save(directory + 'invHessdiag.npy',invHessdiag)
                    
    return None
    
    

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
            
    reconstruct(mag - meanmag,N,R,pindex,pundex,kindex,rho,conv_crit,smooth)

    return 0


if __name__ == '__main__':
    run()
