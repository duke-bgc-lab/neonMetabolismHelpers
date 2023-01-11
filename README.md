# neonMetabolismHelpers
an R package of functions intended to help users modeling stream and river metabolism from National Ecological Observatory Network (NEON) data.  

At the coarsest descritption, these functions execute discrete steps to get data from NEON, prepare and model the data stream metabolism, and evaluate the outputs of the model. There are additional decisions that the user must make in these functions that will affect the metabolism estimates derived. Below, we describe the core of the steps as functions to use this package.

## Naming convention
### Function names
All functions in this package use the prefix nmh_, which is a pseudo-acronym for neonMetabolismHelpers. Following this prefix is a verb (e.g. `get_`, `prep_`, `model_`, etc.) that specifies what is happening in that function. `nmh_get_` means data is being acquired from NEON or elsewhere to fill gaps in the NEON data. `nmh_prep_` prepares either input time series data or model designations.

### Modeling types
We allow the user to select various types of discharge data to serve as inputs to the metabolism model, specificed as `q_type` throughout the package and data. If the user wants to use discharge direct from the NEON data product, we designate this as 'raw'. Based on the results of [Rhea et al. (2022)] [https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/], we recommend using their evaluation of the NEON discharge product, and the evaluation is imported using `get_neon_Q_eval()`. A third option is to use discharge simulations prepared by the MacroSheds project and can be accessed for each site using `get_neon_Q_sim()`. 

## Get data
- specify sites, types, dates
- tracker files
- what are the requisite data to get from NEON (DO, temp, light, bp, Q)
`get_neon_data()`
`get_neon_Q_eval()`
`get_neon_Q_sim()`

## Prepare time series data

## Define gas exchange-discharge relationship

## Model metabolism

## Evaluate and prepare outputs
