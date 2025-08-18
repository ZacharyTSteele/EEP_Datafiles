## Δ′17OBW Models

We generated three models for estimating δ17OBW, δ18OBW, and Δ′17OBW by modifying and adjusting the model and code developed by Steele et al. (2025) and Hu et al. (2023). The Hu et al. (2023) model was derived from the δ18OBW model developed by Kohn (1996), modifying this model to enable calculation of Δ′17OBW. For specific details regarding these models, please see the methods section of the manuscript. 

We have provided the R code for the 'Scaling and Water Flux' model in this repository because the 'Heart Rate and Water Flux' model only differs in the estimation of field metabolic rate, and because the 'Scaling and WEI' model code has been published elsewhere (see Hu et al. [2023]). If you would like access to the R code for the other two models, these can be requested by contacting Dr. Zachary Steele (zacharytaylorsteele@gmail.com). In addiiton, if you would like the Python code for any of the three models, please contact Dr. Steele as well.

The scaling and water flux model is included in the R file 'E17O_isotope_scaling_waterflux_model'. This model estimates δ17OBW, δ18OBW, and Δ′17OBW using the following 27 parameters (in order of appearance in R code):

Metabolic pre-exponent (a in aM^b equation)

Metabolic exponent (b in aM^b equation)

Mass (kg)

Percent of body mass that is water (%)

Percent of body water turnover per day (%)

Fecal H2O content (%)

Oxygen utilization fraction

Z-value (‰)

Animal body temperature (K)

Animal body temperature (C)

Meeh factor (proportionality constant - cm^2/kg^2/3) 

skin H2O evaporation rate (species specific constant - mg/cm^2/h/mm Hg)

Relative humidity (%)

Ambient Temperature (K)

Ambient Temperature (C)

δ18O of meteoric water (‰)

Δ′17O of meteoric water (per meg)

δ18O of atmospheric oxygen (‰)

Δ′17O of atmospheric oxygen (per meg)

NDVI (Normalized difference vegetation index - value betwee -1 to 1)

Relative Digestibility (%)

Energy extraction efficiency (%)

Carbohydrate content of diet (% by mass)

Fat content of diet (% by mass)

Lipid content of diet (% by mass)

Food associated H2O content (%)

Ratio of fractionating plant material (fractionating vs non-fractionating plant tissue)

While there are 27 parameters in total, the majority of the parameters will be fixed for a specific study animal (e.g., skin H2O evaporation rate, Carbohydrate content, Fecal H2O content, etc.) and environment (ambient Temperature, δ18O of meteoric water Δ′17O of meteoric water, etc.). This leaves only a small amount of parameters that will change between individuals (e.g., animal body mass).

Only lines 13-45 of the code that contain the parameters should be adjusted. The rest of the code simply incorporates these parameters to calculate fractional contributions from the different oxygen fluxes, and then uses fractionation factors and exponents and absolute ratios (i.e., of meteoric water & atmospheric oxygen) to estimate δ17OBW, δ18OBW, and Δ′17OBW. 

If you have any questions when using this code and/or data, please email Dr. Steele (ZacharyTaylorSteele@gmail.com).
