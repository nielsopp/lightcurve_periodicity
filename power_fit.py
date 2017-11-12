#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 15:24:45 2017

@author: niels
"""


import matplotlib.pyplot as plt
from nifty import *
import scipy.optimize as so
import sys


#name of directory that results are stored in:
directory = sys.argv[1]

#name of file containing data:
data_file = sys.argv[2]


if directory[-1] != '/':
    directory += '/'


#read lightcurve Gibbs samples and summarize as mean and standard deviation:
timevals = np.load(directory + 'timebins.npy')
try:
    true_lightcurve = np.load(directory + 'true_lightcurve.npy')
except:
    true_lightcurve = np.zeros(timevals.shape)*np.nan
data = np.genfromtxt(data_file)
s = np.zeros(timevals.shape)
ssig = np.zeros(timevals.shape)
count = 0
repeat = True
while repeat:
    try:
        s += np.load(directory + 'm1_%05i.npy'%count)
        s += np.load(directory + 'm2_%05i.npy'%count)
        ssig += (np.load(directory + 'm1_%05i.npy'%count))**2
        ssig += (np.load(directory + 'm2_%05i.npy'%count))**2
        count += 1
    except:
        repeat = False
s /= (2.*count)
ssig /= (2.*count)
ssig -= s**2
ssig = ssig**0.5
s_field = field(rg_space(naxes=1,num=len(s),dist=timevals[1] - timevals[0]),val=s)
s_power = s_field.power(bare=True)


#read critical filter results for lightcurve:
try:
    m_critical = np.load(directory + 'm.npy')
except:
    m_critical = np.zeros(true_lightcurve.shape)*np.nan
try:
    Dhat = np.load(directory + 'Dhat.npy')
except:
    Dhat = np.zeros(true_lightcurve.shape)*np.nan


#load frequencies
freqs = np.load(directory + 'freqs.npy')


#load various power spectra:
try:
    pow_realization = np.load(directory + 'power_realization.npy')
except:
    pow_realization = np.zeros(freqs.shape)*np.nan
try:
    pow_critical = np.load(directory + 'power.npy')
except:
    pow_critical = np.zeros(freqs.shape)*np.nan
try:
    invHessdiag = np.load(directory + 'invHessdiag.npy')
except:
    invHessdiag = np.zeros(freqs.shape)*np.nan


def get_interval(percent,data):
    """Calculate a confidence interval from 0-dimensional posterior samples."""
    dat = data.copy()
    dat.sort()
    numdat = len(dat)
    numdat = numdat // 2
    numdat = numdat - int(np.round(numdat/100.*percent))
    return dat[numdat - 1], dat[-numdat]

    
def get_lower_and_upper(percent,data):
    """Calculate a point-wise confidence interval for functional posterior
    samples."""
    lower = np.zeros(data.shape[1])
    upper = np.zeros(data.shape[1])
    for ii in range(data.shape[1]):
        lower[ii], upper[ii] = get_interval(percent,data[:,ii])
    return lower, upper


#load Gibbs power spectrum samples:
powers = np.zeros((2*count,len(freqs)))
for ii in range(count):
    powers[ii] = np.load(directory + 'power1_%05i.npy'%ii)
    powers[count + ii] = np.load(directory + 'power2_%05i.npy'%ii)
powers = np.array(powers,dtype=float)

#get confidence intervals:
log_lower68, log_upper68 = get_lower_and_upper(68.,np.log(powers))
log_lower95, log_upper95 = get_lower_and_upper(95.,np.log(powers))
meanpow = np.mean(powers,axis=0)


#FITTING:
#parameterized model
def logmodel(params):
    func = freqs.copy()
    func[1:] = params[0]*freqs[1:]**(-params[1]) + params[2]
    func[0] = func[1]
    return np.log(func)

#model fit
logpow = np.median(np.log(powers),axis=0)
fit_params = so.minimize(lambda x: np.sum((logmodel(x) - logpow)**2),x0=[1.e3,1.,1.e3],method='L-BFGS-B',bounds=[(1.e-5,None),(0.1,None),(1.e-5,None)]).x
repeat = True
while repeat:
    indices = ((logmodel(fit_params) < log_upper68) & (logmodel(fit_params) > log_lower68))
    fit_params_new = so.minimize(lambda x: np.sum((logmodel(x)[indices] - logpow[indices])**2),x0=fit_params,method='L-BFGS-B',bounds=[(1.e-5,None),(0.1,None),(1.e-5,None)]).x
    if np.max(np.abs(fit_params_new - fit_params)/fit_params) < 0.01:
        repeat = False
    fit_params = fit_params_new.copy()

#find fraction of posterior samples above a threshold:
def find_frac(ind,thresh):
    return (1.*np.sum(powers[:,ind] > thresh))/powers.shape[0]


#find special frequencies and write them out
f = open(directory + 'special_frequencies.tsv','w')
special_indices = np.arange(len(freqs))[logmodel(fit_params) < log_lower68]
print >> f, 'frequency [per 100 days]\tperiod [days]\tsignificance score'
for ind in special_indices:
    print('frequency %03.3f per 100 days, period %03.3f days, significance score %03.3f'%(freqs[ind]*100,1./freqs[ind],(find_frac(ind,np.exp(logmodel(fit_params)[ind])) - 0.84)/0.16))
    print >> f, '%03.3f\t%03.3f\t%03.3f'%(freqs[ind]*100,1./freqs[ind],(find_frac(ind,np.exp(logmodel(fit_params)[ind])) - 0.84)/0.16)
f.close()

    
#####################################################################################
# Plot is generated below; comment or uncomment the features that you want to include
#####################################################################################
    

#PLOT LIGHTCURVE
meanmag = np.mean(data[:,1])
fig = plt.figure(figsize=(3.3,6.6))
ax = fig.add_axes([0.15,0.57,0.845,0.425])

#plot data
ax.errorbar(data[:,0],data[:,1],data[:,2],linestyle='none',color='green')

#plot true lightcurve
ax.plot(timevals,true_lightcurve,linestyle='-',lw=2,color='cyan',alpha=0.8)

#plot 68% Gibbs sampling interval
ax.fill_between(timevals,s - ssig + meanmag,s + ssig + meanmag,color='blue',alpha=0.3)

#plot Gibbs sampling expectation values
#ax.plot(timevals,s + meanmag,color='blue',lw=3,alpha=0.8)

#plot critical filter estimate
ax.plot(timevals,m_critical + meanmag,color='magenta',lw=3,alpha=0.8)

#plot one-sigma uncertainty estimate for critical-filter estimate
#ax.fill_between(timevals,m_critical + meanmag - Dhat**0.5,m_critical + meanmag + Dhat**0.5,color='magenta',alpha=0.3)

#adjust ticks and labels:
xlocs = np.arange(np.ceil((data[:,0].min() - 500.)/500.)*500,(np.floor((data[:,0].max() + 500.)/500.) + 1)*500,500)
ys = ax.get_yticks()
if np.max(ys) - np.min(ys) < 3.:
    factor = 10.
else:
    factor = 1.
ylocs = np.arange(np.ceil(np.min(ys)*factor),np.floor(np.max(ys)*factor) + 1,2)/factor
plt.xticks(xlocs,fontsize=9)
plt.xlim(data[:,0].min() - 500,data[:,0].max() + 500)
plt.xlabel('time [days]',fontsize=9)
plt.yticks(ylocs,rotation='vertical',fontsize=9)
plt.ylabel('brightness [mag]',fontsize=9)
plt.gca().invert_yaxis()


#PLOT POWER SPECTRUM
ax = fig.add_axes([0.15,0.07,0.845,0.425])

#plot power of the true lightcurve realization
#plt.plot(np.log10(freqs[1:]),np.log10(pow_realization[1:]),color='magenta',linestyle='--')

#plot power of the posterior-mean reconstruction from the Gibbs samples
#plt.plot(np.log10(freqs[1:]),np.log10(s_power[1:]),color='green',linestyle='--')

#plot 95%-confidence interval from Gibbs samples
#plt.fill_between(np.log10(freqs[1:]),np.log10(np.exp(log_lower95[1:])),np.log10(np.exp(log_upper95[1:])),color='blue',alpha=0.15)

#plot 68%-confidence interval from Gibbs samples
ax.fill_between(np.log10(freqs[1:]),np.log10(np.exp(log_lower68[1:])),np.log10(np.exp(log_upper68[1:])),color='blue',alpha=0.3)

#plot posterior-mean power spectrum
#plt.plot(np.log10(freqs[1:]),np.log10(meanpow[1:]),color='blue',lw=3,alpha=0.8)

#plot model fit
ax.plot(np.log10(freqs[1:]),np.log10(np.exp(logmodel(fit_params))[1:]),color='green',lw=3,alpha=0.8,ls='--')

#plot critical-filter reconstruction
ax.plot(np.log10(freqs[1:]),np.log10(pow_critical[1:]),color='magenta',lw=3,alpha=0.8)

#plot critical-filter uncertainty estimate
#ax.fill_between(np.log10(freqs[1:]),np.log10(np.maximum(np.exp(np.log(pow_critical[1:]) - invHessdiag[1:]**0.5),1.e-6)),np.log10(np.exp(np.log(pow_critical[1:] + invHessdiag[1:]**0.5))),color='magenta',alpha=0.3)

#adjust ticks and labels:
xlocs = np.arange(np.ceil(np.min(np.log10(freqs[1:]))),np.max(np.log10(freqs[1:])) + 1,1)
xlabs = [r'$10^{%i}$'%x for x in xlocs]
plt.xticks(xlocs,xlabs,fontsize=9)
plt.xlim(-np.log10(data[:,0].max() - data[:,0].min()),np.max(np.log10(freqs[1:])))
plt.xlabel(r'frequency [1/day]',fontsize=9)
ylocs = np.arange(np.floor(np.min(np.log10(np.exp(logpow[1:])))),np.ceil(np.max(np.log10((meanpow[1:])))) + 1,2)
ylabs = [r'$10^{%i}$'%y for y in ylocs]
plt.yticks(ylocs,ylabs,rotation='vertical',fontsize=9)
plt.ylim(ylocs.min(),ylocs.max())
plt.ylabel(r'power density [mag$^2$ day]',fontsize=9)
for ind in special_indices:
    plt.axvline(np.log10(freqs[ind]),color='grey',lw=5,alpha=0.8,ymin=0,ymax=0.3)
#for ii in [1]:#range(2):
#    plt.axvline(np.log10(ii*4.652e-3),color='brown',lw=2,ls='--',alpha=0.7)
plt.savefig(directory + 'twopanelplot.pdf')
#plt.show()



##########################################
# Summaries of the results are saved below
##########################################


#posterior mean for lightcurve
np.save(directory + 'posterior_mean_lightcurve.npy',s + meanmag)

#one-sigma posterior uncertainty for lightcurve
np.save(directory + 'posterior_uncertainty_lightcurve.npy',ssig)

#critical-filter estimate for lightcurve
np.save(directory + 'critical_filter_lightcurve.npy',m_critical + meanmag)

#one-sigma uncertainty estimate for critical-filter estimate
np.save(directory + 'uncertainty_critical_filter_lightcurve.npy',Dhat**0.5)

#power in the posterior-mean reconstruction
np.save(directory + 'power_of_posterior_mean_lightcurve.npy',s_power)

#lower and upper bounds of 68- and 95-% confidence intervals for power
np.save(directory + 'lower_95perc_power.npy',np.exp(log_lower95))
np.save(directory + 'lower_68perc_power.npy',np.exp(log_lower68))
np.save(directory + 'upper_95perc_power.npy',np.exp(log_upper95))
np.save(directory + 'upper_68perc_power.npy',np.exp(log_upper68))

#posterior-mean power spectrum
np.save(directory + 'post_mean_power.npy',meanpow)

#model fit to the power spectrum
np.save(directory + 'power_spectrum_model_fit.npy',np.exp(logmodel(fit_params)))

#critical-filter reconstruction of the power spectrum
np.save(directory + 'critical_filter_power.npy',pow_critical)

#critical-filter uncertainty estimate for the power-spectrum estimate
np.save(directory + 'critical_filter_power_uncertainty.npy',invHessdiag[1:]**0.5)
