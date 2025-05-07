# Comparison_of_PVK_and_Si_panel_temperatures

To run the present COMSOL Multiphysics simulations with MATLAB:

1. Import AM1.5 Global spectrum data e.g. from https://www.pveducation.org/pvcdrom/appendices/standard-solar-spectra and save it as AM15.mat

2. Create COMSOL model for the solar panels using 
	- SipanelTcomsolModel_v6.m and
	- PSCpanelTcomsolModel_v8.m

3. Set the boundary conditions (environmental conditions and heat generation) and run the simulation using 
	- main_PVKpanelT_VaryingConditions.m or
	- main_PVKpanelT_AllConditions.m


To determine Sandia, Faiman, PVsyst, Mattei, and TRNSYS temperature model parameters:

1. Store the simulated data in a format that only includes the information that is relevant for determining the model parameters using 
	- SaveResults.m

2. Determine the model parameters using
	- TemperatureModels.m


NOTES: 

- In main_PVKpanelT_VaryingConditions.m only one environmental parameter is varied at the time. Except the varied parameter, the ambient temperature is 20°C, wind speed is 1m/s, and irradiance is 800W/m^2. 
- In main_PVKpanelT_AllConditions.m the environmental conditions are varied simultaneously. 
- The temperature model parameters are determined using the data simulated with main_PVKpanelT_AllConditions.m
