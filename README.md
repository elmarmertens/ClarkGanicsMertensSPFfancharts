# REPLICATION FILES for “Constructing Fan Charts from the Ragged Edge of SPF Forecasts“ and “What Is the Predictive Value of SPF Point and Density Forecasts?“ by Clark, Ganics, and Mertens

# WORK IN PROGRESS

This readme file describes the set of replication files (“the replication set“) for two papers:

1. “Constructing Fan Charts from the Ragged Edge of SPF Forecasts“ (henceforth the FanCharts paper)

2. “What Is the Predictive Value of SPF Point and Density Forecasts?“ (henceforth the EntropicTilting paper)  

The replication set contains code as well as all of our input data in raw form as obtained from their original sources described further below. The data and code used by both papers overlap. Specifically, the EntropicTilting paper requires first estimating some of the models discussed in the FanCharts paper. 

### Authors

- Todd E. Clark (Federal Reserve Bank of Cleveland)
- Gergely Ganics (Banco de España and John von Neumann University)
- Elmar Mertens (Deutsche Bundesbank) [Corresponding author: [em@elmarmertens.com](mailto:em@elmarmertens.com)]

## Data

Data used for this project comprises SPF survey responses as well as realized values for up to five different macroeconomic variables. Additional details are also described in Sections 3 and 2 of the FanChart and EntropicTilting papers, respectively. All data were obtained from two, publicly available online sources: 

1. The [Real-Time Data Research Center (RTDRC)](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research) at the Federal Reserve Bank of Philadelphia.

2. The [FRED database](https://fred.stlouisfed.org) hosted by the Federal Reserve Bank of St. Louis.

#### Model estimation

From the RTDRC we obtained SPF mean responses for the following variables (with SPF mnemonics in parentheses as listed on their [SPF website](https://www.philadelphiafed.org/surveys-and-data/data-files)):

- Level of real GDP/GNP (RGDP)
- Level of the price index for GDP/GNP (PGDP)
- Civilian unemployment rate (UNEMP)
- CPI inflation rate (CPI)
- 3-month Treasury bill rate (TBILL)

In addition, for GDP growth and GDP inflation, we collected real-time measures of quarter $t-1$ data as these data were publicly available in quarter $t$ from the quarterly files of real-time data in the Philadelphia Fed's [Real-Time Data Set for Macroeconomists (RTDSM)](https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history) [https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/real-time-data-set-full-time-series-history].

From FRED, we collected realized values for CPI, UNRATE and TBILL (mnemonics: CPIAUCSL, UNRATE, TB3MS) using ”final” vintage data.

#### Forecast evaluation

For forecast evaluation, we measure the outcomes of GDP growth and GDP inflation with the RTDSM vintage published two quarters after the outcome date (that is, we use the quarterly vintage in $t + h + 2$ to evaluate forecasts for $t+h$ made in $t$; this is the second estimate available in the RTDSM’s vintages). Because revisions to quarterly data are relatively small for the unemployment rate and CPI inflation and non-existent
for the T-bill rate, we simply use the historical time series available in the St. Louis Fed’s FRED database to measure the outcomes and corresponding forecast errors for these variables.

The replication set includes copies of the raw input files in `kensingtonDataUSSPF`. Below we also describe code that transforms the input data before further processing by our main estimation routines.

In the ET paper, we only used data on real GDP, the GDP price index and the unemployment rate, as described in Section 2 therein. Additionally, we used the annual fixed-event probability forecasts for these variables from the [RTDRC website]([Probability Variables: Survey of Professional Forecasters].

## Code

All code used for this project was written in MATLAB. The code has been run on various, recent MATLAB versions (Versions 2019a through 2022b) as well as different operating systems (Linux, Windows and macOS) without need for any particular adjustments across platforms. The codes uses MATLAB’s Statistics and Machine Learning Toolbox toolbox as well as (optionally) the Parallel Computing Toolbox. 

The MATLAB code also creates LaTeX files collecting tables and figures produced by the MATLAB code. If a LaTeX installation is present (and if the “pdflatex” command is available on the command line via MATLAB’s “system” command), the LaTeX files will also be compiled into PDF files.

The replication code comes in the following subdirectories: 

- `kensingtonDataUSSPF` contains raw data files obtained from the sources described above as well as MATLAB files for transforming the raw inputs into a set of `.mat` files (one for each of the variables).  As part of the replication files, copies of these `.mat` files are also contained in the following subdirectory.
- `kensingtonDataMatfiles` contains the `.mat` files generated by the scripts in `kensingtonDataUSSPF`. 
- `kensingtonMCMC` contains various scripts and functions to perform MCMC estimation of our baseline model as well as the various alternatives described in the FanChart and EntropicTilting papers and their appendices.
- `kensingtonET` contains various scripts to perform entropic tilting using the model outputs generated in `kensingtonMCMC` (only relevant for the EntropicTilting paper).
- `kensingtonTablesAndFigures` contains various scripts to collect results (as generated by code provided in the `kensingtonMCMC` directory) and produce tables and figures. 
- `kensingtonToolbox` contains a few helper functions specific to this project
- `matlabtoolbox` contains a github submodule with various folders providing different auxiliary `.m`-files (MATLAB scripts and functions) used throughout.  The toolboxes are automatically loaded onto the MATLAB path upon invocation of any of the scripts contained in the previously described code directories. Please note that some toolbox files were obtained either from the [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/)[https://www.mathworks.com/matlabcentral/fileexchange/] or [James P. Le Sage's Econometrics Toolbox](https://www.spatial-econometrics.com) [https://www.spatial-econometrics.com], which are both freely available; please see the comment headers of the respective toolbox files for further attribution.


### To prepare input data files for estimation

The directory `kensingtonDataUSSPF` contains all of the raw data files obtained from the RTDRC and FRED described above. In addition, the directory contains several `.m`-file scripts to transform  the raw data input files for further processing by our main estimation routines contained in `kensingtonMCMC`. All `.m`-file scripts create `.mat` data files in MATLAB format.

To process raw data for RGDP and PGDP (both of which are matched against realized values collected from the RTDRC), and CPI, UNRATE, TBILL (all three of which are matched against realized values collected from FRED) please run `kensingtonCollectData.m`. For each variable, a data file is created and stored in MATLAB’s `.mat` format. (Resulting `.mat` files are also provided as part of the replication set and stored in the `kensingtonDataMatfiles` directory.) 

In case of updating the data, please adapt the end of the sample in the appropriate scripts as indicated by comments therein (see variable SPFsamstop in line 29). 

NOTES:

1. The EentropicTilting paper only uses data on RGDP, PGDP and UNRATE.

2. The FanChart paper's supplementary appendix presents results of a model variant (labeled 'SV-AVG10') which, in addition to using the SPF forecasts up to maximum 3 years ahead, also utilizes long-term (10-year) SPF forecasts of the RGDP, CPI and TBILL variables (average growth rates for the former two variables, while average level in the case of the last variable). The corresponding raw SPF forecasts are processed by `kensingtonCollectDataRGDP10.m`, `kensingtonCollectDataCPI10.m` and `kensingtonCollectDataTBILL10.m`, respectively.

### To prepare input data files for entropic tilting (only for EntropicTilting paper)

In addition to preparing the input data for model estimation as described in the previous section, the replicating the results in the EntropicTilting paper requires some more additional data, detailed below.

The script `kensingtonCollectProbData.m` collects and processes the SPF probability forecasts for RGDP, PGDP and UNRATE (mnemonics: PRGDP, PRPGDP and PRUNEMP), saving the processed data in `kensingtonPROB.mat`, to be used in the entropic tilting exercise. It also computes the Discrete Ranked Probability Score (DRPS) associated with the SPF histograms - these are saved on a variable-by-variable basis in files named `kensingtonVARNAMEdataHISTOGRAMS.mat`, where VARNAME is RGDP, PGDP or UNRATE.

The script `kensingtonFitDistributionSPF.m` fits generalized beta and normal distributions to the SPF histograms, and saves the results in `kensingtonFitDistSPF.mat`.

### To estimate the various models

Code for estimating the various model variants considered in the papers and appendices is provided in `kensingtonMCMC`. In addition, as part of the replication set, `kensingtonDataMatfiles` contains copies of the input data `.mat` files as created in `kensingtonDataUSSPF`. When updating the data, please copy the updated `.mat` files into `kensingtonDataUSSPF`.

Run the following MATLAB scripts:

- `kensingtonSTATEtrendgapSV.m` estimates the baseline SV model in real-time, as described in the FanChart and EntropicTilting papers. The script loops over all five SPF variables.

- `kensingtonSTATEconst.m` estimates the alternative CONST model in real-time, as described in the FanChart and EntropicTilting papers. The script loops over all five SPF variables.

- `kensingtonVARSTATEtrendgapSV.m` estimates the SV model without the MDS assumption in real-time, as described in the FanChart paper. The script loops over all five SPF variables.

- `kensingtonSTATEtrendgaplongSV.m` estimates the SV model using 10-year average forecasts (labeled 'SV-AVG10') in real-time, as described in the supplementary appendix of the FanChart paper. The script loops over RGDP, CPI and TBILL.
  
  ### To perform entropic tilting (EntropicTilting paper only)
  
  Code for entropic tilting is provided in `kensingtonET`.
  
  Run the following MATLAB scripts:
  
  - `kensington_ET.m` performs entropic tilting of the baseline SV and CONST models. The script loops over the three SPF variables (RGDP, PGDP, UNRATE) and the tilting variants considered in the paper and the supplementary appendix. The script assumes that the relevant mat files generated earlier by the scripts in
    - `kensingtonDataUSSPF` are available in `kensingtonDataMatfiles`,
    - `kensingtonMCMC` contains code to estimate the various state space models by MCMC
  - `kensington_ET_eval.m` processes the results of entropic tilting, calculating statistics (e.g., forecast errors, CRPS, etc.).
  
  ### To generate tables and figures
  
  The directory `kensingtonTablesAndFigures` provides scripts to generate LaTeX tables and figures (as shown in our papers and the appendices). These scripts assume that MCMC results generated by the code in `kensingtonMCMC` and stored in `.mat` files are stored in a common directory.
  
To generate the tables and figures presented in the EntropicTilting paper, the results generated by the scripts in `kensingtonET` and stored in `.mat` files are needed. 

