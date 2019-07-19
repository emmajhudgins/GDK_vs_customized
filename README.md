# Comparing generalized to customized invasive pest dispersal models

This folder includes the data and code necessary to replicate Hudgins et al. (2019) "Comparing generalized to customized models for United States invasive forest pest spread". Ecological Applications. 

The code was written by Emma J. Hudgins, PhD candidate, McGill University

For any questions, e-mail me at emma.hudgins@mail.mcgill.ca

This code is licensed under the GPL v3.0. Please cite this GitHub folder and the associated publication upon re-use.

The folder contains 4 main R scripts along with associated functions called by each script.

1) *customized_models.R* - this fits the fully context specific models for each of the 3 focal species based on the best fit to year 2000. Validation is then performed on 2005.

2) *GDK_3spp.R* - this fits and/or validates the various levels of GDK correction (including GDKu for each of the 3 focal species)

3) *GDK_forecast.R* - this fits (no validation) the GDKic (and GDKic+c/GDKs+ic+c in a couple of cases) for the remaining species in the dataset. Species #3 (Asian longhorn beetle) should be removed from analysis based on the Corrigendum to Hudgins et al. (2017)

4) *plot_GDK_vs_custom_models.R* - this is the plotting script to generate predicted and observed distribution maps in 2005 (validation year) based on the fits of various models. Example input is given for HWA's best models (GDKu, GDKs+ic+c, customized).  

The first 2 scripts call r2_om_mse.R (which itself calls om_mse.R) during the validation phase. These scripts can optionally call random_nn_spp_2019.R in order to generate new null model draws for which to calculate R2om, but this is a time-consuming process, so output has been provided.

All data files are described in the script that calls them.
