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

An example bash script that runs on the example data contained in the txt file (corresponding to Fig. 5 of the paper) is given as lightcurve_analysis.sh. Note that, in its default version, this script deletes the samples of the Gibbs-sampling analysis after they have been analyzed.


## Parameters

* GIBBS_DIR: The analysis scripts will create a directory by this name and store their output there.
* GIBBS_DATA_FILE: The name of the file containing the data to be analyzed (see below for its structure).
* GIBBS_NPIX: The number of time bins to use for the reconstructed lightcurve.
* GIBBS_EXTENSION: This is a factor by which the time interval for the reconstruction is extended beyond the length of the observations. This is done to reduce the effect of perodicity enforced through the use of Fast Fourier Trasforms.
* GIBBS_BURN_IN: The number of samples to discard at the beginning of the two Gibbs-sampling chains. This parameter should be estimated whenever studying a new setup.
* GIBBS_CORR_LEN: An estimated correlation length for the Gibbs chains. The analysis will regard samples spaced this number apart as independent. This parameter should also be estimated anew when studying a new observational setup or drastically changing any other paramter.
* GIBBS_CONV_CRIT: This number defines the stopping criterion for the Gibbs sampling. Two chains are started at different random locations. The mean natural logarithm of the power spectra for each chain is calculated after each correlation length. Once these two means differ by less than GIBSS_CONV_CRIT for each frequency, the sampling is stopped. While the sampling runs, the code will print out the sample number and the current value of this quantity.
* CRIT_CONV_CRIT: This number defines the stopping criterion for the critical-filter iteration. The power spectrum is updated in each iteration. Once the difference between the natural logarithms of two consecutive power spectra is smaller than CRIT_CONV_CRIT, the iteration is stopped. While the iteration runs, the code will print out the iteration number and the current value of this quantity.
* CRIT_SMOOTH: Strength of the spectral-smoothness prior that prevents the power spectrum from dropping to zero in the critical filter. This number is akin to a prior variance, so larger numbers mean a weaker enforcement of smoothness.


## Input

The data input is provided as a text file. Each row represents an observation and should consist of three numbers, separated by whitespace or commas:
* The time of the observations (assumed to be in days for plotting purposes; otherwise any units are fine).
* The observed brightness (assumed to be in magnitudes for plotting purposes; otherwise any units are fine).
* The error bar on the observed brightness (assumed to be in magnitudes for plotting purposes; otherwise the unit just has to match that of the observed brightness).
Note that the observational uncertainties are assumed to be Gaussian in the units provided and independent between observations. 


## Output

The python scripts produce output in *.npy files, each containing a one-dimensional array.

### Output of Gibbs.py

* timebins.npy: Time values corresponding to the grid points used in the reconstruction.
* freqs.npy: Frequencies corresponding to the Fourier grid points.
* power1.npy (power2.npy): Current mean of the natural logarithm of the power spectrum from the first (second) Gibbs-sampling chain. This is updated while the code runs. The comparison of the two means defines convergence (see above).
* power1_xxxxx.npy (power2_xxxxx.npy): Power-spectrum samples from the first (second) chain. Only the samples beyond the burn-in phase and spaced a correlation length apart are saved.
* m1_xxxxx.npy (m2_xxxxx.npy): Lightcurve samples from the first (second) chain. Only the samples beyond the burn-in phase and spaced a correlation length apart are saved. Note that the mean of the data is subtracted.


### Output of critical_filter.py

* timebins.npy: Time values corresponding to the grid points used in the reconstruction.
* freqs.npy: Frequencies corresponding to the Fourier grid points.
* power.npy: Critical-filter estimate of the power spectrum. This is updated as the iteration progresses.
* m.npy: Critical-filter estimate of the lightcurve. This is updated as the iteration progresses. Note htat the mean of the data is subtracted.
* Dhat.npy: Point-wise posterior variance estimate for the lightcurve.
* invHessdiag.npy: Point-wise posterior variance estimate for the power spectrum.


### Output of power_fit.py

* posterior_mean_lightcurve.npy: Posterior mean of the lightcurve (based on Gibbs sampling).
* posterior_uncertainty_lightcurve.npy: Point-wise one-sigma uncertainty estimate for the lightcurve (based on Gibbs sampling).
* critical_filter_lightcurve.npy: Lightcurve estimate from the critical filter.
* uncertainty_critical_filter_lightcurve.npy: One-sigma uncertainty estimate for the lightcurve based on the critical filter.
* power_of_posterior_mean_lightcurve.npy: The power contained in the posterior-mean lightcurve estimate (based on Gibbs sampling). This will in general have less power than the estimated power spectrum.
* post_mean_power.npy: Posterior mean of the power spectrum (based on Gibbs sampling).
* lower_95perc_power.npy (lower_68perc_power.npy): Lower bound of the frequency-wise 95-percent (68-percent) confidence interval for the power spectrum (based on Gibbs sampling).
* upper_95perc_power.npy (upper_68perc_power.npy): Upper bound of the frequency-wide 95-percent (68-percent) confidence interval for the power spectrum (based on Gibbs sampling).
* power_spectrum_model_fit.npy: Model fit to the power spectrum results of the Gibbs sampling, as described in the paper.
* critical_filter_power.npy: Critical-filter estimate for the power spectrum.
* critical_filter_power_uncertainty.npy: Frequency-wise uncertainty estimate for the critical-filter power spectrum.


## Contact

For questions/comments/bugs, please contact niels@cita.utoronto.ca
