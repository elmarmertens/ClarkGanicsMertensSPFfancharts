# CLARK-GANICS-MERTENS REPLICATION FILES

This readme file describes the set of replication files (“the replication set“) for two papers:

1. “Constructing Fan Charts from the Ragged Edge of SPF Forecasts“ (henceforth the StateSpace paper)[^DOIstsp]

2. “What Is the Predictive Value of SPF Point and Density Forecasts?“ (henceforth the EntropicTilting paper)  [^DOIet]

[^DOIstsp]: [latest draft](https://drive.google.com/file/d/1CxA4ho1lBi3s0piCgbbD-vPL92tGv8nq/view?usp=share_link), with [supplementary appendix](https://drive.google.com/file/d/1c34QIQQw08Avnojou-9ypSMu5a2LBUXG/view?usp=share_link), also available as Federal Reserve Bank of Cleveland Working Paper 22-36 [https://doi.org/10.26509/frbc-wp-202236](https://doi.org/10.26509/frbc-wp-202236)

[^DOIet]: [latest draft](https://drive.google.com/file/d/1Ofrqan6wXHHAMQHVA3TtnUQLIq7AdNcA/view?usp=share_link), with [supplementary appendix](https://drive.google.com/file/d/1gdoTaOB_MW3rN_CKDGb7P9YnS9-N5AI1/view?usp=share_link), also available as Federal Reserve Bank of Cleveland Working Paper 22-37 [https://doi.org/10.26509/frbc-wp-202237](https://doi.org/10.26509/frbc-wp-202237)

The replication set contains code as well as all of our input data in raw form as obtained from their original sources described further below. The data and code used by both papers overlap. Specifically, the EntropicTilting paper builds on estimates generated as part of the StateSpace paper. However, due to the availability of SPF density forecasts (needed for the EntropicTilting paper) the StateSpace paper covers also more economic variables than the EntropicTilting paper. 

The project is work in progress, and all results are to be considered preliminary.  The materials provided do not necessarily reflect the views of the Federal Reserve Bank of Cleveland, the Federal Reserve System, the Banco de España, the Deutsche Bundesbank, or the Eurosystem.

## Authors

- Todd E. Clark (Federal Reserve Bank of Cleveland)
- Gergely Ganics (Banco de España and John von Neumann University)
- Elmar Mertens (Deutsche Bundesbank) [^em] 

[^em]: Corresponding author: [em@elmarmertens.com](mailto:em@elmarmertens.com)

## Overview (both papers)

### Data

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
- 3-month Treasury bill rate (TBILL)

In addition, for GDP growth and GDP inflation, we collected real-time measures of quarter $t-1$ data as these data were publicly available in quarter $t$ from the quarterly files of real-time data in the Philadelphia Fed's [Real-Time Data Set for Macroeconomists (RTDSM)](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history)[^rtdsm].

[^rtdsm]: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history

From FRED, we collected realized values for CPI, UNRATE and TBILL (mnemonics: CPIAUCSL, UNRATE, TB3MS) using ”final” vintage data.

For forecast evaluation, we measure the outcomes of GDP growth and GDP inflation with the RTDSM vintage published two quarters after the outcome date (that is, we use the quarterly vintage in $t + h + 2$ to evaluate forecasts for $t+h$ made in $t$; this is the second estimate available in the RTDSM’s vintages). Because revisions to quarterly data are relatively small for the unemployment rate and CPI inflation and non-existent
for the T-bill rate, we simply use the historical time series available in the St. Louis Fed’s FRED database to measure the outcomes and corresponding forecast errors for these variables.

The replication set includes copies of the raw input files in the folder [`kensingtonDataUSSPF`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataUSSPF). Below we also describe code that transforms the input data before further processing by our main estimation routines.

In the ET paper, we only used data on real GDP, the GDP price index and the unemployment rate, as described in Section 2 therein. Additionally, we used the annual fixed-event probability forecasts for these variables from the [SPF page at the RTDRC](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/probability-variables)[^spfprob], and the RECESS variable (probability that real GNP/GDP will decline) available on the [SPF website](https://www.philadelphiafed.org/surveys-and-data/recess)[^spfrecess].

[^spfprob]: https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/probability-variables
[^spfrecess]: https://www.philadelphiafed.org/surveys-and-data/recess

### Code

All code used for this project was written in MATLAB. The code has been run on various, recent MATLAB versions (Versions 2019a through 2022b) as well as different operating systems (Linux, Windows and macOS) without the need for any particular adjustments across platforms. The code uses MATLAB’s Statistics and Machine Learning Toolbox toolbox as well as (optionally) the Parallel Computing Toolbox. 

The MATLAB code also creates LaTeX files collecting tables and figures produced by the MATLAB code. If a LaTeX installation is present (and if the “pdflatex” command is available on the command line via MATLAB’s “system” command), the LaTeX files will also be compiled into PDF files.

The replication code comes in the following subdirectories: 

- [`kensingtonDataUSSPF`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataUSSPF) contains raw data files obtained from the sources described above as well as MATLAB files for transforming the raw inputs into a set of `.mat` files (one for each of the variables). 
- [`kensingtonDataMatfiles`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataMatfiles) is a placeholder, to be populated by the `.mat` files generated by the scripts in `kensingtonDataUSSPF`. 
- [`kensingtonMCMC`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonMCMC) contains various scripts and functions to perform MCMC estimation of our baseline model as well as the various alternatives described in the  StateSpace and EntropicTilting papers and their appendices.
- [`kensingtonET`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonET) contains various scripts to perform entropic tilting using the model outputs generated in `kensingtonMCMC` (only relevant for the EntropicTilting paper).
- [`kensingtonTablesAndFigures`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonTablesAndFigures) contains various scripts to collect results (as generated by code provided in the `kensingtonMCMC` directory) and produce tables and figures. 
- [`kensingtonToolbox`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonToolbox) contains a few helper functions specific to this project.
- [`matlabtoolbox`](https://github.com/elmarmertens/em-matlabbox) contains a GitHub submodule with various folders providing different auxiliary `.m`-files (MATLAB scripts and functions) also available at https://github.com/elmarmertens/em-matlabbox, with further documentation found there. The scripts used in this repository load the appropriate toolbox folders automatically onto the MATLAB path as needed.

### To prepare input data files for estimation

The directory `kensingtonDataUSSPF` contains all of the raw data files obtained from  RTDRC and FRED described above. In addition, the directory contains several `.m`-file scripts to transform  the raw data input files for further processing by our main estimation routines contained in `kensingtonMCMC`. All `.m`-file scripts create `.mat` data files in MATLAB format.

To process raw data for RGDP and PGDP (both of which are matched against realized values collected from the RTDRC), and CPI, UNRATE, TBILL (all three of which are matched against realized values collected from FRED) please run `kensingtonCollectData.m`. For each variable, a data file is created and stored in MATLAB’s `.mat` format. (The resulting `.mat` should be copied into the `kensingtonDataMatfiles` directory for further use by the estimation routines described below.) 

In case of updating the data, please adapt the end of the sample in the appropriate scripts as indicated by comments therein (see variable `SPFsamstop` in the script `kensingtonCollectData.m` and its variants). 

Note that the EntropicTilting paper uses only data on RGDP, PGDP and UNRATE. Preparation of SPF histogram data, specific to the EntropicTilting paper, is described further below.

### To generate tables and figures

The directory [`kensingtonTablesAndFigures`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonTablesAndFigures) provides scripts to generate LaTeX tables and figures (as shown in our papers and the appendices). These scripts assume that MCMC results generated by the code in [`kensingtonMCMC`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonMCMC) and stored in `.mat` files are stored in a common directory; please edit the function `localresultsMCMC.m` to return a string that points to that directory.[^local]

To generate the tables and figures presented in the EntropicTilting paper, the results generated by the scripts in [`kensingtonET`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonET) and stored in `.mat` files are also needed, in addition to the MCMC results. Please edit the function `localresultsET.m` to return a string that points to the directory that contains the EntropicTilting results.[^local]

[^local]: The scripts  `localresultsMCMC.m` and  `localresultsET.m` are part of the matlabtoolbox that is provided as submodule to this repository.

## Specifics: StateSpace paper

Code for estimating the various model variants considered in our project is provided in [`kensingtonMCMC`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonMCMC). In addition,  the replication set assumes that [`kensingtonDataMatfiles`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataMatfiles) contains copies of the input data `.mat` files as created in [`kensingtonDataUSSPF`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataUSSPF). When updating the data, please copy the updated `.mat` files into [`kensingtonDataMatfiles`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataMatfiles).

### To estimate the models

- `kensingtonSTATEtrendgapSV.m` estimates the baseline SV model in real-time, as described in the StateSpace and EntropicTilting papers. The script loops over all five SPF variables.

- `kensingtonSTATEconst.m` estimates the alternative CONST model in real-time, as described in the  StateSpace and EntropicTilting papers. The script loops over all five SPF variables.

- `kensingtonVARSTATEtrendgapSV.m` estimates the SV model without the MDS assumption in real-time, as described in the StateSpace paper. The script loops over all five SPF variables.

- `kensingtonSTATEtrendgaplongSV.m` estimates the SV model using 10-year average forecasts (labeled 'SV-AVG10') in real-time, as described in the supplementary appendix of the StateSpace paper. The script loops over RGDP, CPI and TBILL.

The StateSpace paper's supplementary appendix also presents results of a model variant (labeled 'SV-AVG10') which, in addition to using the SPF forecasts up to maximum 3 years ahead, also utilizes long-term (10-year) SPF forecasts of the RGDP, CPI and TBILL variables (average growth rates for the former two variables, while average level in the case of the last variable). The corresponding raw SPF forecasts are processed by `kensingtonCollectDataRGDP10.m`, `kensingtonCollectDataCPI10.m` and `kensingtonCollectDataTBILL10.m`, respectively.

### To generate tables and figures

Note: tables and figures listed below without a prefix refer to the Cleveland Fed Working Paper version, and tables and figures with the prefix “S“ referring to the Supplementary materials.

- `tabulateRelativeEvalstats.m`: relative forecast accuracy of CONST and SV models, both using MDS assumption (Table 2), and SV model with and without MDS assumption (Table 3, Table S.1).
- `figuresETA4SV.m`: time-varying volatility estimates and absolute value of expectational updates of SV model (Figure 1, Figures S.3 to S.6).
- `figuresYtermstructure.m`: real-time term structures of expectations based on SV model (Figures 2 to 5). To plot the full-sample estimates, run `figuresYtermstructureFinal.m`.
- `figuresYuncertainty.m`: term structure of uncertainty based on SV model (Figure 6).
- `figuresRMSEVAReta.m`: time-varying bias and MSE decomposition of SPF-consistent expectations in SV model without MDS assumption (Figure 7, Figures S.13 to S.19).
- `figuresSEPfancharts.m`: SEP-style annual fan charts from SV model (Figure 8).
- `figuresSEPfanchartsUncertainty.m`: annual uncertainty measures from SEP and SV model (Figure 9).
- `goSTATEtrendgapSV.m`(in `kensingtonMCMC`): MA coefficients of univariate representation of $y_t$ (Figure S.1), and estimates of square root of SV factor (Figure S.2).
- `figuresYpredictivedensities.m`: forecast fan charts based on SV model (Figures S.7 to S.12).
- `figuresLongrunEstimates2.m`: shifting endpoints of term structures of expectations, SV and CONST models (Figure S.20).
- `figuresLongrunEstimatesAVG10.m`: shifting endpoints of term structures of expectations, SV and SV-AVG10 models (Figure S.21).
- `figuresYpredictivedensities2.m`: forecast fan charts based on SV and SV-AVG10 models (Figures S.22 and S.23).
- `tabulateRelativeEvalstatsAVG10.m`: relative forecast accuracy of SV and SV-AVG10 models (Tables S.2 and S.3).

## Specifics: EntropicTilting paper

### To prepare input data files for entropic tilting

In addition to preparing the input data for model estimation as described in the previous section, replication of results from the EntropicTilting paper requires some more additional data, detailed in this section.

The script `kensingtonCollectProbData.m` collects and processes the SPF probability forecasts for RGDP, PGDP and UNRATE (mnemonics: PRGDP, PRPGDP and PRUNEMP), saving the processed data in `kensingtonPROB.mat`, to be used in the entropic tilting exercise. It also computes the Discrete Ranked Probability Score (DRPS) associated with the SPF histograms - these are saved on a variable-by-variable basis in files named `kensingtonVARNAMEdataHISTOGRAMS.mat`, where VARNAME is RGDP, PGDP or UNRATE.

The script `kensingtonFitDistributionSPF.m` fits generalized beta and normal distributions to the SPF histograms, and saves the results in `kensingtonFitDistSPF.mat`.

### To perform entropic tilting

  Code for entropic tilting is provided in [`kensingtonET`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonET).

  Run the following MATLAB scripts:

- `kensington_ET.m` performs entropic tilting of the baseline SV and CONST models. The script loops over the three SPF variables (RGDP, PGDP, UNRATE) and the tilting variants considered in the paper and the supplementary appendix. The script assumes that the relevant `.mat` files generated earlier by the scripts in [`kensingtonDataUSSPF`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataUSSPF) are available in [`kensingtonDataMatfiles`](https://github.com/elmarmertens/ClarkGanicsMertensSPFfancharts/tree/main/kensingtonDataMatfiles), and the MCMC results should be stored in a common directory, whose name is provided by calls to the function `localresultsMCMC.m` (please edit this function as needed).[^local]
- `kensington_ET_eval.m` processes the results of entropic tilting, calculating statistics (e.g., forecast errors, CRPS, etc.).

### To generate tables and figures

Note: tables and figures listed below without a prefix refer to the Cleveland Fed Working Paper version, and tables and figures with the prefix “S“ referring to the Supplementary materials.

- `kensingtonTable1.m`: forecast accuracy comparison tables (SV and CONST), tilting to SPF histograms (Tables 2 and 3, Tables S.2 and S.14).
- `tabulateRelativeEvalstatsETSPFCDF.m`: forecast accuracy comparison table, SV and CONST models, tilting to SPF histograms vs. tilting to moments of fitted generalized beta distribution (Table 4, Tables S.7 to S.11 and S.16 to S.20).
- `SPFmeanRanges.m`: comparison of annual SPF point forecasts vs. range of means compatible with histograms (Figure 1, Figures S.1 and S.2).
- `figuresCDFybar4.m`: impact of entropic tilting on cumulative distribution functions (SV and CONST) along with cumulative SPF probabilities (Figures 2 and 3).
- `figuresYuncertaintyET.m`: term structures of uncertainty of SV model with and without entropic tilting (Figure 4, Figure S.14).
- `tabulateRelativeDRPS.m`: DRPS of models (SV and CONST) and SPF for annual forecasts (Figure 5).
- `tabulateRelativeDRPScurrent.m`: DRPS of models (SV and CONST) and SPF for annual forecasts, _including_ current-year predictions (Tables S.3, S.4 and S.15, Figures S.3, S.4 and S.15).
- `kensingtonTable1since2009.m`: forecast accuracy comparison tables (SV and CONST), tilting to SPF histograms since the Great Recession (Tables S.5 and S.6).
- `figuresYpredictivedensitiesET.m`: forecast fan charts of SV model with and without entropic tilting (Figures S.5 to S.7, S.12 and S.13).
- `figuresRECESS.m`: recession probabilities and average cumulative DRPS of SPF, and SV and CONST models with and without entropic tilting (Figures S.8 and S.9).
- `figuresFitCDFSPF.m`: skewness of generalized beta distribution fitted to SPF histograms (Figures S.10, S.11 and S.16).
- `kensingtonTableSVBins.m`: forecast accuracy comparison tables including CONST model tilted to probabilities assigned to SPF bins from SV model (Tables S.12, S.13 and S.21).
