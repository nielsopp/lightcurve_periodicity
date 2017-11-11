#!/bin/bash

#PARAMETERS FOR BOTH GIBBS AND CRITICAL-FILTER RECONSTRUCTIONS:

#directory name for output
export GIBBS_DIR=test

#data file name
export GIBBS_DATA_FILE=test_data.txt

#number of time bins for reconstruction
export GIBBS_NPIX=256

#factor by which to extend reconstruction beyond baseline
export GIBBS_EXTENSION=4


#PARAMETERS SPECIFIC TO GIBBS SAMPLING:

#length of burn-in phase in samples
export GIBBS_BURN_IN=200

#correlation length in samples
export GIBBS_CORR_LEN=50

#convergence criterion quantity
export GIBBS_CONV_CRIT=2.0


#PARAMETERS SPECIFIC TO CRITICAL-FILTER RECONSTRUCTION:

#convergence criterion quantity
export CRIT_CONV_CRIT=0.01

#strength of the spectral smoothness prior
export CRIT_SMOOTH=100000.0


#############################################################

python Gibbs.py $GIBBS_DIR $GIBBS_DATA_FILE $GIBBS_NPIX $GIBBS_EXTENSION $GIBBS_BURN_IN $GIBBS_CORR_LEN $GIBBS_CONV_CRIT

#############################################################

python critical_filter.py $GIBBS_DIR $GIBBS_DATA_FILE $GIBBS_NPIX $GIBBS_EXTENSION $CRIT_CONV_CRIT $CRIT_SMOOTH

#############################################################

python power_fit.py $GIBBS_DIR $GIBBS_DATA_FILE

#############################################################

#CLEAN UP SAMPLES
rm $GIBBS_DIR/m1_?????.npy
rm $GIBBS_DIR/m2_?????.npy
rm $GIBBS_DIR/power1_?????.npy
rm $GIBBS_DIR/power2_?????.npy
rm $GIBBS_DIR/power1.npy
rm $GIBBS_DIR/power2.npy
