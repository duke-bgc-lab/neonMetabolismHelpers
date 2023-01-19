# neonMetabolismHelpers

an R package of functions intended to help users modeling stream and river metabolism from National Ecological Observatory Network (NEON) data.

At the coarsest description, these functions execute discrete steps to get data from NEON, prepare and model the data stream metabolism, and evaluate the outputs of the model. There are additional decisions that the user must make in these functions that will affect the metabolism estimates derived. Below, we describe the core of the steps as functions to use this package.

## Install and load

    devtools::install_github(duke-bgc-lab/neonMetabolismHelpers)
    library(neonMetabolismHelpers)

## Naming conventions

### Function names

All functions in this package use the prefix `nmh_`, which is a pseudo-acronym for neonMetabolismHelpers. Following this prefix is a verb (e.g. `nhm_get_`, `nhm_prep_`, `nhm_model_`, etc.) that specifies what is happening in that function. `nmh_get_` means data is being acquired from NEON or elsewhere to fill gaps in the NEON data. `nmh_prep_` prepares either input time series data or model designations.

### Modeling types

We allow the user to select three types of discharge data to serve as inputs to the metabolism model, specified as `q_type` throughout the package, data, and file structure.

-   'raw': If the user wants to use discharge direct from the NEON data product, we designate this as 'raw'. As NEON does not collect discharge at the 'TOMB' site, we exclude this site from the 'raw' designation

-   'qaqc': based on the results of [Rhea et al. (2022)](https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/), we recommend using their evaluation of the NEON discharge product, and the evaluation is imported using `get_neon_Q_eval()`.

-   'simulation': this uses discharge simulations prepared by the [MacroSheds](https://github.com/MacroSHEDS/macrosheds) project and can be accessed for each site using `nmh_get_neon_Q_sim()`. More details on these simulations can be found [add link here].

## Get data

The following functions obtain data requisite to model stream metabolism. The requisite data from NEON are found in the following data types and are downloaded using `nmh_get_neon_data()` . Users can specify which:

-   sites to download data from using the `sites = ...` argument

-   dates to download data from with the `date_range = ...` argument

The default for is to download all possible data, which may take quite a long time, so be careful! The raw data will be downloaded into a new directory called data/raw, and sub-directories for each Data Product Name (see below), and for each DP specified by site (e.g. data/raw/Continuous discharge/ARIK).

| Variable(s)                         | NEON Data Product Code                                                                                                                                                           | NEON Data Product Name                               |
|------------------|-----------------------------------|------------------|
| Discharge                           | [DP4.00130.001](https://data.neonscience.org/data-products/DP4.00130.001#:~:text=NEON%20calculates%20continuous%20stream%20discharge,Elevation%20of%20surface%20water%2C%20DP1.) | Continuous discharge                                 |
| Light                               | [DP1.20042.001](https://data.neonscience.org/data-products/DP1.20042.001)                                                                                                        | Photosynthetically active radiation at water surface |
| Barometric pressure                 | [DP1.00004.001](https://data.neonscience.org/data-products/DP1.00004.001)                                                                                                        | Barometric pressure                                  |
| Water temperature, Dissolved oxygen | [DP.202288.001](https://data.neonscience.org/data-products/DP1.20288.001)                                                                                                        | Water quality                                        |

Alternatively, `nmh_get_neon_data()` can be used to download any NEON data product indiviudally with the same site and date range arguments; this is effectively a wrapper around `neonUtilities::loadByProduct()` .

If the user needs to only load one or more of these data products, individual data products exist to get those products individually, if needed:

-   `nmh_get_neon_wq()`, `nmh_get_neon_discharge()`, `nmh_get_neon_bp()`, `nmh_get_neon_light()`

The remaining functions relate to how discharge is processed, discussed above:

-   `nmh_get_neon_Q_eval()` : downloads the results from Rhea et al. (2022) as a csv file
-   `nmh_get_neon_Q_sim()` : downloads the simulated discharge from the NEON sites from the MacroSheds portal

## Prepare time series data

After the data from NEON and/or MacroSheds has been downloaded, the data are prepared into sub-daily (e.g. 15 minute) time series files for each site using `nmh_prep_metab()` that will serve as inputs to `streamMetabolizer` . `streamMetabolizer` requires input time series with as a data frame with required column names: `solar.time`, `DO.obs`, `DO.sat`, `temp.water`, `discharge`, `depth`, `light` . Therefore, the goal of `nmh_prep_data()` is to perform the formatting, necessary calculations, and data access where needed and output csv files for each site.

The time-series starts with the dissolved oxygen and water temperature data, as these are the core data in modeling metabolism. In the data paper associated with this R package, we have taken additional steps to process these data before modeling. Data were downloaded from NEON and uploaded to the [StreamPULSE portal](https://data.streampulse.org/index), and were visually QAQC'd in addition to the quality flags (QF) that NEON populates in their data products. The QAQC'd data were accessed from the `StreamPULSE` R package and the following steps were applied.

From these data, saturation concentration of dissolved oxygen (`DO.sat`) is calculated using barometric pressure. If barometric pressure is available for a site and date from NEON (see above), we use those data; if not, we use an internal function in the `streamPULSE` R package which searches for the nearest NOAA weather station and accesses the barometric pressure at that station. With these data, we use `streamMetabolizer::calc_DO_sat(…, salinity = 0, model = 'garcia-benson')` to estimate `DO.sat`.

Next, discharge is added to the time series. Depending on the source of discharge, the handling of this data is different.

-   For `q_type = 'raw'` , the discharge time series is aggregated to mean 15-minute discharge from the 1-minute resolution data from NEON and converted from L s^-1^ to m^3^ s^-1^.

-   For `q_type = 'qaqc'` , the evaluation from Rhea et al. (2022) is applied and site-months that fail a rating curve threshold are excluded, resulting in a shorter time series from certain sites, but of known higher data quality. In addition, we include discharge data from the 'TOMB' site, as these data have passed USGS data checks. This is completed by running `nmh_get_tomb_q()`, which accesses USGS database using `dataRetrieval::readNWISuv('02459761', c('00060', '00065', startDate = '2015-01-01)` . These data are converted from cfs to m^3^ s^-1^ and appended to the time-series from the TOMB site

-   For `q_type = 'simulation'` , the time-series from the simulations are called into the function and appended to the time-series

Finally, light is added to the time-series. If light exists from NEON, those 1 minute data are added the 15 sum of PAR. In addition, the timestamp in UTC is converted to `solar.time` using `streamMetabolizer::convert_UTC_to_solartime(…, type = 'mean')`. If no light data exist, light is estimated using `streamMetabolizer::calc_light()` and the latitude + longitude of each site.

At this final step in building the time series, mean reach depth is estimated using `streamMetabolizer::calc_depth()`. This calculatation estimates the mean depth of the reach upstream of the sensor that contribute to the benthic area that has influence on the sensor. The coefficients, c and f, can be specific using the `z_method` argument. If `z_method == 'model'`, values of (`c = 0.409, f = 0.294`) are estimated in [Raymond et al. (2012)](https://aslopubs.onlinelibrary.wiley.com/doi/10.1215/21573689-1597669) for the continental US. Alternatively, if `z_method == 'meas'`, site-specific values are available through processing the width measurements in the [Reaeration field and lab collection](DP1.20190.001) and instantaneous depth, velocity, and discharge measurements in the [Discharge field collection](https://data.neonscience.org/data-products/DP1.20048.001) NEON data products. The site-specific coefficients for 24/27 sites can be found in `data/MacroSheds/NEON_site_scaling_coefs.csv` . For the sites where there are no measaurements (FLNT, TOMB, BLWA) we use the default. Additionally, where the R^2^ of the power law scaling is \<0.1, we use the default values.

As with the raw data, this function writes these output csv files into a directory. In the same data folder where the raw data are saved, each site csv are saved to `data/sm_input/`, and depending on the discharge data used, the next directory will be raw, qaqc, or simulation. Each file will be saved as `{site}_{qtype}_smReady.csv` (e.g. `data/sm_input/raw/ARIK_raw_smReady.csv`)

### Plot time-series \*this function is incomplete (1/12/23)

`nmh_plot_metab_inputs()` will plot the input time series and allow a comparison between different q_types included by the user. These files will be saved to a new directory `figures/sm_input/` .

## Define gas exchange-discharge relationship

We will first run the time-series model to estimate and define a gas exchange-discharge relationship for each site (K600 \~ Q). This is handled in the function `nmh_model_metab_mle()` .

The inputs for this function are the directory of the prepared time series. Within the function, we define the model name and specifications as a wrapper around `streamMetabolizer::mm_name()` and `streamMetabolizer::specs()`. More details about these functions can be found in [Appling et al. (2018)](https://www.nature.com/articles/sdata2018292) and the [streamMetabolizer repository](https://github.com/USGS-R/streamMetabolizer).

Models are run for each site-year (e.g. ARIK-2018, ARIK-2019, etc.). Output data from the model is stored in a new folder: `data/model_runs/{q_type}/`. The output data include the model fits, the model specs, and the predicted sub-daily DO.

There are internal checks and handling of each year time-series and model building that are populated in a csv output to `data/model_runs/q_type/MLE/neon_{q_type}_MLE_results.csv`. Each row has 5 columns: site, year, and Has Data (y/n), Fit Error (y/n), and Fit time (minutes).

The results of this process will become priors in the Bayesian model run. These priors are extracted using `nmh_eval_metab_mle()` . For each site-year a MLE model was fit, this function extracts the median K600 value, among other descriptive statistics, to use in the Bayesian model.

## Model metabolism

-   Before modeling, the following two functions should be run.

    -   `nmh_prep_bayes_model()` : checks for installation of `streamMetabolizer` v0.12.0 and `rstan`. If those packages are not installed, it will do so.

    -   `nmh_prep_site_year()`: Create a data frame with each site-year intended to model that is looped over in the model function. This function returns a data frame `site_years` directly to the environment

        -   `years = c('2016', '2021')` : specify the range of years intended to use. 2016 and 2021 are the defaults

        -   `site = 'all'` : this calls `get_neon_site_data()` to obtain all 27 sites, but can be subset to a string of the desired sites (e.g. `site = c('ARIK', 'FLNT')` )

The function that runs the model, `nmh_model_metab_bayes()`, is extensive and is broken into 5 steps, which we summarize here. All of these processes are handled within the function and rely on the suggested directory structure presented above.

-   Step 0: Tracking system

    -   For each site year, a tracker text file is written. This file records important meta-data and system information and tracks the results of the subsequent steps. Tracker files are written to `data/model_runs/{q_type}/Bayes/trackers`

-   Step 1: Data checking

    -   The models are run for each unique site-year in a `for` loop and referenced in `nmh_prep_site_year()`. For each site-year, the data for a given site are read in and filtered to the given year of interest.

    -   If there is no data for the specific site-year, the tracker file records there is no data and jumps to the next step in the loop. If there is data, a summary of the input data is run, calculating the % of NAs for each column in the input data file. This summary is printed to the tracker file

    -   If there are sufficiently low NAs, the data is cleared to enter the model

-   Step 2: Model name, priors, and specs

    -   The model name is defined outside the for loop using `streamMetabolizer::mm_name()` , though the specific arguments we recommend not changing.

    -   Gas exchange-discharge priors are defined by reading in the results of `nmh_model_metab_mle()` and filtering to the site-year of interest. Details on these priors can be found in the "Metabolism model configuration and application" section in [Appling et al. (2018)](10.1038/sdata.2018.292).

    -   Priors are passed to `streamMetabolizer::specs()` with additional arguments we recommend not changing. If there is an error in defining these specs, an error is flagged in the tracking file. Otherwise, the specs are ready to be input to the model

-   Step 3: Running the model

    -   The model is run by `streamMetabolizer::metab()` with inputs the time-series data (Step 1) and the specs (Step 2)

    -   If there is an error in fitting the model, a flag is recorded in the tracker file and internal log matrix. If the model run successfully, the time to run the model is recorded in the internal log matrix

-   Step 4: Saving the outputs

    -   The function records 6 outputs from `streamMetabolizer::metab()` : daily, overall, KQ_overall, specs, datadaily, mod_and_obs_DO. Each of these become a subdirectory in `data/model_runs/` and a csv for each site-year is written to each folder after successfully model run.

-   Step 5: Printing and saving the internal log matrix

    -   The internal log matrix, which serves as the data frame that is looped around, is printed with a brief summary of the information stored in the tracker files. This is written to `data/model_runs/{q_type}/Bayes/neon_bayes_row.csv` .

We highly recommend running the Bayesian models in parallel when possible. The models take a very long time to fit and converge and `nmh_model_metab_bayes()` takes advantage of parallelization. Below is example code to run the model in parallel that should be run each time, including and especially the `stopCluster()`, even if the model run was incomplete.

    # Before modeling

    nmh_prep_bayes_model()  # check installation of streamMetabolizer and rstan
    nmh_prep_site_years()   # create data.frame of site years to loop around

    library(foreach) 
    library(doParallel)

    n.cores <- parallel::detectCores()

    if(.Platform$OS.type == 'windows') {

    cl <- parallel::makeCluster(n.cores, type = 'PSOCK')
    } else {
      cl <- parallel::makeCluster(n.cores, type = 'FORK')
    }

    registerDoParallel(cl)
    foreach::getDoParWorkers()

    # Run the model
    nmh_model_metab_bayes()

    # After the model
    stopCluster(cl)
    unregister_dopar <- function() {
        env <- foreach:::.foreachGlobals
        rm(list = ls(env), pos = env)
    }
    unregister_dopar()

## Evaluate and prepare outputs

The following functions exist to parse and return an object of high confidence metabolism estimates.

-   `nmh_diagnose_metab_bayes()` extracts MCMC convergence of site-year models and the "real-world" results from each year (number of days GPP \> -0.5 g O2 m^-2^ d^-1^ or ER \< 0.5 g O2 m^-2^ d^-1^).

-   `nmh_estimate_metab_bayes()` returns "real-world" results with confidence bounds for GPP and ER, as well as model convergence stats, merged with relevant data from the input time series (e.g. DO min/max/range, mean depth, mean discharge, light inputs) and combines all sites into a single csv file.

-   `nmh_filter_metab_bayes()` reads in the output of `nmh_diagnose_metab_bayes()`, which identifies site-years with models that did not converge or had a frequency of "real-world" results lower than a 60% threshold, and removes those data with poor model fits from the output file returned by `nmh_estimate_metab_bayes()` .
