# lightcurve_periodicity

Code for analyzing power spectra of lightcurves via Gibbs sampling and IFT's critical-filter formalism to extract special frequencies.

The implemented algorithms are described in [arXiv:2017.????](http://www.arxiv.org/abs/????). Please reference this paper when making use of the code.


## Prerequisites

The python scripts contained in this repository rely heavily on the following packages:
* [numpy](http://www.numpy.org)
* [scipy](https://www.scipy.org)
* [nifty](http://wwwmpa.mpa-garching.mpg.de/ift/nifty/)
* [matplotlib](http://matplotlib.org) (for final plotting purposes only)

Once these are installed, they should run out-of-the box under python 2.6 and above.


## Code pieces

The code consists of three separate python scripts, each taking parameters described below.
* Gibbs.py runs the Gibbs-sampling analysis of a data file.
* critical_filter.py runs the critical-filter analysis of a data file.
* power_fit.py extracts special periods from the Gibbs reconstruction as described in the paper and produces plots in the style of the paper figures. This assumes that Gibbs sampling has been run; the critical-filter results are optional.

An example bash script that runs on the example data contained in the tsv file is given as lightcurve_analysis.sh. Note that, in its default version, this script deletes the samples of the Gibbs-sampling analysis after they have been analyzed.


## Parameters

* GIBBS_DIR: The analysis scripts will create a directory by this name and store their output there.
* GIBBS_DATA_FILE: The name of the file containing the data to be analyzed (see below for its structure).
* GIBBS_NPIX: The number of time bins to use for the reconstructed lightcurve.
* GIBBS_EXTENSION: This is a factor by which the time interval for the reconstruction is extended beyond the length of the observations. This is done to reduce the effect of perodicity enforced through the use of Fast Fourier Trasforms.
* GIBBS_BURN_IN: The number of samples to discard at the beginning of the two Gibbs-sampling chains. This parameter should be estimated whenever studying a new setup.
* GIBBS_CORR_LEN: An estimated correlation length for the Gibbs chains. The analysis will regard samples spaced this number apart as independent. This parameter should also be estimated anew when studying a new observational setup or drastically changing any other paramter.
* GIBBS_CONV_CRIT: This number defines the stopping criterion for the Gibbs sampling. Two chains are started at different random locations. The mean natural logarithm of the power spectra for each chain is calculated after each correlation length. Once these two means differ by less than GIBSS_CONV_CRIT for each frequency, the sampling is stopped.
* CRIT_CONV_CRIT: This number defines the stopping criterion for the critical-filter iteration. The power spectrum is updated in each iteration. Once the difference between the natural logarithms of two consecutive power spectra is smaller than CRIT_CONV_CRIT, the iteration is stopped.
* CRIT_SMOOTH: Strength of the spectral-smoothness prior that prevents the power spectrum from dropping to zero in the critical filter. This number is akin to a prior variance, so larger numbers mean a weaker enforcement of smoothness.


## Input

The data input is provided as a text file. Each row represents an observation and should consist of three numbers, separated by whitespace or commas:
* The time of the observations (assumed to be in days for plotting purposes; otherwise any units are fine).
* The observed brightness (assumed to be in magnitudes for plotting purposes; otherwise any units are fine).
* The error bar on the observed brightness (assumed to be in magnitudes for plotting purposes; otherwise the unit just has to match that of the observed brightness).
Note that the observational uncertainties are assumed to be Gaussian in the units provided and independent between observations. 


## Output

