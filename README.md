## Δ′17OBW Models

We generated three models for estimating δ17OBW, δ18OBW, and Δ′17OBW by modifying and adjusting the model and code developed by Steele et al. (2025) and Hu et al. (2023). The Hu et al. (2023) model was derived from the δ18OBW model developed by Kohn (1996), modifying this model to enable calculation of Δ′17OBW. For specific details of the differences between models included in this repository, please see the methods section of the manuscript. 

The first model is included in the R file 'E17O_isotope_scaling_waterflux_model'. This model estimates δ17OBW, δ18OBW, and Δ′17OBW using the following 32 parameters (in order of appearance in R code):

Metabolic pre-exponent (a in aM^b equation)

Metabolic exponent (b in aM^b equation)

Mass (kg)

Fecal H2O content (%)

Oxygen utilization fraction

Z-value (‰)

Animal body temperature (K)

Animal body temperature (C)

Meeh factor (proportionality constant - cm^2/kg^2/3) 

skin H2O evaporation rate (species specific constant - mg/cm^2/h/mm Hg)

Inhaled oxygen (mL)
Exhaled CO2 (mole)

Relative humidity (%)

Ambient Temperature (K)

Ambient Temperature (C)

δ18O of meteoric water (‰)

Δ′17O of meteoric water (per meg)

δ18O of atmospheric oxygen (‰)

Δ′17O of atmospheric oxygen (per meg)

Food consumed (kg)

Water intake (mL)

Relative Digestibility (%)

Drinker Nozzle Leakage (%)

Energy extraction efficiency (%)

Carbohydrate content of diet (% by mass)

Fat content of diet (% by mass)

Lipid content of diet (% by mass)

Food associated H2O content (%)

Fractional contribution from different food items

δ17O of different food items

δ18O of different food items

Fractional contribution of different plant components of diet

δ18O of different plant components of diet

While there are 30 parameters in total, the majority of the parameters will be fixed for a specific study animal (e.g., skin H2O evaporation rate, Carbohydrate content, Fecal H2O content, etc.) and environment (ambient Temperature, δ18O of meteoric water
Δ′17O of meteoric water, etc.). This leaves only a small amount of parameters that will change between individuals (e.g., animal body mass). In addition, the captive nature of our experiment resulted in several parameters being included that would otherwise be estimated using metabolic rate (e.g., food intake, inhaled oxygen, exhaled CO2) when studying free-ranging animals (see Hu et al. 2023).

Only lines 13-71 of the code that contain the parameters should be adjusted. The rest of the code simply incorporates these parameters to calculate fractional contributions from the different oxygen fluxes, and then uses fractionation factors and exponents and absolute ratios (i.e., of meteoric water & atmospheric oxygen) to estimate δ17OBW, δ18OBW, and Δ′17OBW. 

The data that is provided in the 'ExperimentData' csv can be used to recreate all our predicted δ17OBW, δ18OBW, and Δ′17OBW values using the Δ′17OBW  Model code. This csv should be used in combination with the information presented in Tables 1 & 3 in the manuscript that cover the fixed parameters for the species and captive setting, and some of the parameters that changed between experimental manipulations. 

The 'ExperimentData' consists of seven columns relevant to the different sampling periods of the experiment: 
MouseID - The ID for the specific mouse in the experiment (1-8) 
Diet - If the mouse was consuming the control (0.49% NaCl) or 4% NaCl (experimental) diet. The diet consumed is also relevant to the timeframe of the code (if the mouse was eating the 4%NaCl diet, then the turnover will be shorter [2 days], whereas if the animal was eating the 0.49% NaCl diet the turnover will be longer [4 days]. 
Temp - What temperature (C) the mouse was housed at during that sampling (either 25, 18, or 15C). 
O2_Inhaled - The amount of O2 inhaled (mL) during the timeframe 
CO2_Exhaled - The amount of CO2 exhaled (mL) during the timeframe. This should be converted to moles by converting to L and then dividing by 22.4 
Food_Consumed - The amount of food consumed (kg) during the timeframe 
Drinking_Water - The water intake (mL) during the timeframe 
Mass - The average body mass of the animal during the timeframe

If you have any questions when using this code and/or data, please email Zach Steele (ZacharyTaylorSteele@gmail.com).
