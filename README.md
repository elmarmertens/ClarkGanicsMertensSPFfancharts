# CLARK-GANICS-MERTENS REPLICATION FILES

This readme file describes the set of replication files (“the replication set“) for  “Constructing Fan Charts from the Ragged Edge of SPF Forecasts.“ 

Replication code for an earlier version of this project, alongside codes for a companion paper, titled “What Is the Predictive Value of SPF Point and Density Forecasts?“ can be found in the branch "EntropicTiltingPaper2022" of this repository.

The replication set contains code as well as all of our input data in raw form as obtained from their original sources described further below. 

The project is work in progress, and all results are to be considered preliminary.  The materials provided do not necessarily reflect the views of the Federal Reserve Bank of Cleveland, the Federal Reserve System, the Banco de España, the Deutsche Bundesbank, or the Eurosystem.

Recent copies of the paper and its supplementary online appendix are posted at [www.elmarmertens.com](https://www.elmarmertens.com)

## Authors

- Todd E. Clark (Federal Reserve Bank of Cleveland)
- Gergely Ganics (Banco de España)
- Elmar Mertens (Deutsche Bundesbank) [^em] 

[^em]: Corresponding author: [em@elmarmertens.com](mailto:em@elmarmertens.com)

## Data files

Data used for this project comprises SPF survey responses as well as realized values for up to five different macroeconomic variables. Additional details are also described in Sections 3 and 2 of the StateSpace and EntropicTilting papers, respectively. All data were obtained from two, publicly available online sources: 

1. The [Real-Time Data Research Center (RTDRC)](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research)[^rtdrc] at the Federal Reserve Bank of Philadelphia.

2. The [FRED database](https://fred.stlouisfed.org)[^fred] hosted by the Federal Reserve Bank of St. Louis.

[^rtdrc]: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research
[^fred]: https://fred.stlouisfed.org

From the RTDRC we obtained SPF mean responses for the following variables (with SPF mnemonics in parentheses as listed on the [SPF website](https://www.philadelphiafed.org/surveys-and-data/data-files)[^spf]):

[^spf]: https://www.philadelphiafed.org/surveys-and-data/data-files

- Level of real GDP/GNP (RGDP)
- Level of the price index for GDP/GNP (PGDP)
- Civilian unemployment rate (UNEMP)
- CPI inflation rate (CPI)


In addition, for GDP growth and GDP inflation, we collected real-time measures of quarter $t-1$ data as these data were publicly available in quarter $t$ from the quarterly files of real-time data in the Philadelphia Fed's [Real-Time Data Set for Macroeconomists (RTDSM)](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history)[^rtdsm].

[^rtdsm]: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history

From FRED, we collected realized values for CPI and UNRATE and  (mnemonics: CPIAUCSL and UNRATE) using ”final” vintage data.

For forecast evaluation, we measure the outcomes of GDP growth and GDP inflation with the RTDSM vintage published two quarters after the outcome date (that is, we use the quarterly vintage in $t + h + 2$ to evaluate forecasts for $t+h$ made in $t$; this is the second estimate available in the RTDSM’s vintages). Because revisions to quarterly data are relatively small for the unemployment rate and CPI inflation and non-existent
for the T-bill rate, we simply use the historical time series available in the St. Louis Fed’s FRED database to measure the outcomes and corresponding forecast errors for these variables.

The replication set includes copies of the raw input files in the folder [`rawdataKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/rawdataKensington). Below we also describe code that transforms the input data before further processing by our main estimation routines.


### Code

All code used for this project was written in MATLAB. The code has been run on various, recent MATLAB versions (Versions 2022a through 2024a) as well as different operating systems (Linux, Windows and macOS) without the need for any particular adjustments across platforms. The code uses MATLAB’s Statistics and Machine Learning Toolbox toolbox as well as (optionally) the Parallel Computing Toolbox. 

The MATLAB code also creates LaTeX files collecting tables and figures produced by the MATLAB code. If a LaTeX installation is present (and if the “pdflatex” command is available on the command line via MATLAB’s “system” command), the LaTeX files will also be compiled into PDF files.

The replication code comes in the following subdirectories: 

- [`rawdataKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/rawdataKensington) contains raw data files obtained from the sources described above as well as MATLAB files for transforming the raw inputs into a set of `.mat` files (one for each of the variables). 
- [`matdataKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/matdataKensington) is a placeholder, to be populated by the `.mat` files generated by the scripts in `rawdataKensington`. 
- [`mcmcKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/mcmcKensington) contains various scripts and functions to perform MCMC estimation of our models.
- [`tablesandfiguresKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/tablesandfiguresKensington) contains various scripts to collect results (as generated by code provided in the `mcmcKensington` directory) and produce tables and figures. 
- [`matlabtoolbox`](https://github.com/elmarmertens/em-matlabbox) contains a GitHub submodule with various folders providing different auxiliary `.m`-files (MATLAB scripts and functions) also available at https://github.com/elmarmertens/em-matlabbox, with further documentation found there. The scripts used in this repository load the appropriate toolbox folders automatically onto the MATLAB path as needed.

### To prepare input data files for estimation

The directory `rawdataKensington` contains all of the raw data files obtained from  RTDRC and FRED described above. In addition, the directory contains several `.m`-file scripts to transform  the raw data input files for further processing by our main estimation routines contained in `mcmcKensington`. All `.m`-file scripts create `.mat` data files in MATLAB format.

To process raw data for all variables (RGDP, PGDP, CPI, and UNRATE) please run `collectAllData.m`. For each variable, a data file is created and stored in MATLAB’s `.mat` format. (The resulting `.mat` should be copied into the `matdataKensington` directory for further use by the estimation routines described below.) 

In case of updating the data, please adapt the end of the sample in the appropriate scripts as indicated by comments therein (see variable `SPFsamstop` in the script `collectData.m` and its variants). 

### To estimate the models

Our main estimation routines are contained in `mcmcKensington`.

- `kensingtonMDStrendHScycleSVt2blockNoiseHS.m` generates out-of-sample estimates from our MDS model with noise in measurement equations for annual forecasts.

- `kensingtonVAR0trendHScycleSVt2blockNoiseHS.m` generates out-of-sample estimates from our VAR model with noise in measurement equations for annual forecasts.

- `kensingtonMDStrendHScycleSVt2block.m` and `kensingtonVAR0trendHScycleSVt2block.m` generates out-of-sample estimates from the noise-free versions of MDS and VAR models.

- Various `go*.m` files allow to estimate MDS and VAR models over a single sample, or to compare output from estimates of two model versions.


### To generate tables and figures

The directory [`tablesandfiguresKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/tablesandfiguresKensington) provides scripts to generate LaTeX tables and figures (as shown in our papers and the appendices). These scripts assume that MCMC results generated by the code in [`mcmcKensington`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/mcmcKensington) and stored in `.mat` files are stored in a common directory (indicated by the variable `resultsdir` in each script).



