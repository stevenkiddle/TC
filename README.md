# tempCluster

Temporal Clustering Documentation

Steven Kiddle (steven.kiddle@kcl.ac.uk)



1 - Repeating paper analysis

To repeat the analysis in Kiddle et al., you first need to download the data. This cannot be shared as it belongs to the relevant consortia, whose terms and conditions you will need to accept. Therefore section (1a) provides instructions on accessing this data, while section (1b) provides instructions on repeating the analysis.


1a - access to AAC data

AAC is a combination of ADNI, AIBL and CAMD cohorts. CAMD is itself a combination of placebo arms of several trials. To get access to AAC data you therefore need you to apply to ADNI, AIBL and CAMD seperately for data. It may take several weeks to be approved to access these datasets. 

When downloaded the data can be combined using a script that I have provided. This file expects the data to be in the temporalClustering/data folder, in folders called:

adni_sep_2015
aibl_28apr_2015
camd_june_2015


1a(i) ADNI and AIBL

Apply for access to ADNI and AIBL at:

http://adni.loni.usc.edu/data-samples/access-data/

From ADNI download:

APOERES.csv
DXSUM_PDXCONV_ADNIALL.csv
MMSE.csv
PTDEMOG.csv
UPENNBIOMK2.csv

From AIBL download 

aibl_ptdemog_28-Apr-2015.csv
aibl_mmse_28-Apr-2015.csv


1a(ii) CAMD

Apply for access to CAMD at:

http://c-path.org/programs/camd/

By clicking on ‘Tools and teams’, and then on ‘Request a login’.

It may be easiest to download the whole database, but you need:

qs.csv
dm.csv


1b - Running analysis

When the relevant files are arranged in the relevant folders you can then navigate to:

temporalClustering/scripts

and run:

source(’TCanalysis.R’)

if this command doesn’t work then try replacing the quotation marks with your own and repeating.

This re-runs most of the analyses, with random seed set to 1 at the beginning, making it reproducible. Elements of the code that were run on the cluster are given, but for runtime we provide pre-computed results that are summarised by the code above. Due to the nature of running this on the cluster in parellel, set seed was not used for these (K2 bootstrap, and all simulations).

Results are generated in the temporalClustering/results folder. Total runtime may be several hours long.# TC
