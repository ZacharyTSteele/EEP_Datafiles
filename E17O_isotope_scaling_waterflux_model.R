## alpha1 and alpha 2 are functions needed for calculating fractionation factors 
## please see Majoube 1971; Friedman & O’Neil 1977;O’Neil & Adami 1969
alpha1<-function(x){
  a1<-exp((1.137*10^6/(x^2)-4.156*10^2/x-2.0667)/1000)
  return(a1)}

exp((1.137*10^6/(288.15^2)-4.156*10^2/188.15-2.0667)/1000)

alpha2<-function (x) {
a2<-exp((1.66*10^4/x-15.69)/1000)
return(a2)}

# Animal-Specific Parameters
Met_pre_exp<- 2.96;                              # Metabolic pre-exponent (see Appendix C Kohn 1996)
Met_exp<- 0.73;                                  # Metabolic exponent (see Appendix C Kohn 1996)
Mass<- 59                                        # Animal Mass (kg) 
Perc_Water<-0.667                                # Percent of body mass that is water
Turnover_Fraction<-.167                          # Percent of body water turnover per day (see appendix of Nagy & Peterson 1988)
FecalH2O<- 0.5;                                  # Fecal H2O content (%)
Frac_used_O2<- 0.20;                             # Fraction of O2 used *oxygen utilization fractions (Johnsen et al. 1990; Joyce and Baxter 1964)
Z_value<- 10.5;                                  # Z-value *Kohn 1996 - see Zanconato et al. 1992
Body_T<- 37+273.15;                              # Body temperature (K)
Temp<-37                                         # Body temperature (C)
meeh_factor<-1100                                # Proportionality constant known as the Meeh factor (According to Dawson and Hulbert [1970] 
                                                 # k is often assumed to be ~ 1000, and has values of 1100-1200 for a variety of marsupials)
Skin_evap_constant<-0.027                        # transcutaneous water vapor loss (assumed value using Camelus spp.; Schmidt-Nielsen 1969)

# Environmental Parameters
rh<- 0.61;                                       # Relative humidity
T<- 273.15+ 7.04;                                # Environment temperature (K)
T_c<-7.04                                        # Environment temperature (C)
d18O_mw<- -12.3;                                 # Meteoric water [pre-formed water] d18O composition (per mil)
D17O_mw<- 0.029;                                 # Meteoric water D17O composition (per mil)
d18O_atmO2<- 24.046;                             # Atmospheric O2 d18O composition (Wostbrock et al. 2020)
D17O_atmO2<- -0.441;                             # Atmospheric O2 D17O composition (Wostbrock et al. 2020)
NDVI<- 0.323                                     # Normalized difference vegetation index

# Food Parameters
Digest<- 0.7;                                    # Relative Digestibility (see Kohn 1996 or Robbins 1983)
Energy_effi<- 0.9;                               # Energy extraction efficiency (see Robbins 1983)
Carb_content<- 0.85;                             # Food carbon content (Kohn 1996)
Fat_content<- 0.05;                              # Food fat content (Kohn 1996)
Protein_content<- 0.1;                           # Food protein content (Kohn 1996)
H2Oassoc_content<- 0.65;                         # Food associated H2O content (Kohn 1996)
ratio_frac_plant_material<- 0.75;                # Ratio of fractionating plant material - see Methods & Table 1 of manuscript

######## Do not Change Anything after this line ###############

#         Constants         #

lamda<- 0.528;              # Reference lamda value for water system
R_18_SMOW<- 0.0020052;      # O-18 isotope ratio, standard SMOW
R_17_SMOW<- 0.00038;        # O-17 isotope ratio, standard SMOW

#                Leaf Water Model                     #
#          By J. Roden, modified by C. Tai            #
#          Incorporated from model in Hu et al. 2023  #

# Inputs for leaf water
a_vap_mw<- 1/alpha1(T)       
# Fractionation factor between water vapor and meteoric water (Function:alpha1)
d18O_atm_vap<- a_vap_mw*(d18O_mw+1000)-1000;   
# Atmospheric humidity, d18O
RH<- rh*100;                  
# Relative humidity (%)
T_diff<- 0;                   
# Difference between leaf and air temperature (K)
baro_press<- 86.83;          
# Barometric pressure (kPa)
stoma_cond<- 0.3;             
# Stomatal conductance range
boud_cond<- 1;                
# Boundary layer conductance

# Constants for models
alpha_k<- 1.032;              
# Kinetic fracionation (Cappa et al. 2003)
alpha_kb<- 1.021;             
# Kinetic fractionation at boundary layer (Cappa et al. 2003)
eO_auto<- 27;                 
# Autotrophic fractionation (sucrose) (Yakir & DeNiro 1990; Luo & Sternberg 1992; Sternberg & DeNiro 1983)
eO_hetero<- 27;               
# Heterotrophic fractionation (Cellulose) (Yakir & DeNiron 1990; Luo & Sternberg 1992; Sternberg & DeNiro 1983)
fO_ex_mew<- 0.42;             
# Fraction of exchange with medium water (Roden & Ehleringer 1999)

# Calculations
Leaf_T<- T+T_diff;            
# Leaf temperature calculated
e_sat<- (101325*exp((((-0.1299*(1-(373.15/T))-0.6445)*(1-(373.15/T))-1.976)*(1-(373.15/T))+13.3185)*(1-(373.15/T))))/1000;
# Saturation vapor pressure (kPa)
e_a<- rh*e_sat;               
# Ambient vapor pressure (kPa)
e_i<- (101325*exp((((-0.1299*(1-(373.15/Leaf_T))-0.6445)*(1-(373.15/Leaf_T))-1.976)*(1-(373.15/Leaf_T))+13.3185)*(1-(373.15/Leaf_T))))/1000;
# Leaf vapor pressure (kPa)
g<- 1/(1/stoma_cond+1/boud_cond); 
# Total leaf conductance to water vapor (mol/m2/s)
E<- ((e_i-e_a)/baro_press)*g; 
# Leaf transpiration (mol/m2/s)
W_i<- e_i/baro_press;         
# Water vapor mole fraction
W_a<- e_a/baro_press;         
# Ambient water vapor mole fraction
W_s<- ((stoma_cond*W_i)-E*(1-W_i/2))/(stoma_cond-E/2);
# Leaf surface water vapor
e_s<- W_s*baro_press;         
# Vapor pressure at leaf surface (kPa)
alpha_equi<- alpha1(T);       
# Equilibrium fractionation (liquid-vapor) at leaf temperature (Function alpha1) 
R_source<- R_18_SMOW*(1+d18O_mw/1000);  
# O isotope ratio of source water
R_atm_vap<- R_18_SMOW*(1+d18O_atm_vap/1000);
# O isotope ratio of atmospheric water vapor
R_Leaf<- alpha_equi*(alpha_k*R_source*(e_i-e_s)/e_i+alpha_kb*R_source*(e_s-e_a)/e_i+R_atm_vap*e_a/e_i);
# O isotope ratio of leaf water

d18O_Leaf<- (R_Leaf/R_18_SMOW-1)*1000;
# Modeled leaf water d18O compositions (per mil)
d18O_medw<- 0.5*(d18O_mw+d18O_Leaf);
# Modeled medium water for cellulose synthesis (per mil) 
d18O_Cellulos<- fO_ex_mew*(d18O_medw+eO_hetero)+(1-fO_ex_mew)*(d18O_Leaf+eO_auto);
# Modeled cellulose d18O compositions (per mil)

# Energy Equation #

Energy<- (10^Met_pre_exp)*(Mass^Met_exp);
# Animal Total Metabolic Energy kJ

# Inputs #

# Ingested Food #
Food_energy<- (Carb_content*17.3+Fat_content*39.7+Protein_content*20.1)*Digest*Energy_effi*1000;
# Food metabolic energy kJ/kg
Food_O<- Carb_content*11.12
# Moles of O from per kg food (mole/kg)
Food_CO2<- Carb_content*11.12+Fat_content*2+Protein_content*3
# Moles of CO2 from per kg food (mole/kg)
Food_consumed<- (Energy/Food_energy)*(1+(NDVI-0.2));
# Total consumed food (kg) (different from actual ingested food)
Condensation_Rxn_H2O<- Food_consumed*Food_O*Digest*Energy_effi;                            
# Food O content ingested (mole) 
Food_CO2_content<- Food_consumed*Food_CO2*Digest*Energy_effi
# Food CO2 content ingested (mole) 
Decarboxylation_Ox<-Food_CO2_content*2;
#Food CO2 O content 
Food_H2<- Carb_content*30.9+Fat_content*60+Protein_content*11;
#  Moles of H2 from per kg food (mole/kg)
Food_H2_content<- Food_consumed*Food_H2*Digest*Energy_effi;
# Food H2 content ingested (mole)
Food_H2_O2<- Food_H2_content/2-Condensation_Rxn_H2O;                            
# Food H2/2- Food O2 (mole)
Food_H2O_assoc<- Food_consumed*55.56/(1-H2Oassoc_content)*H2Oassoc_content;
# Influx of unbound food H2O - 55.6 is the mole of ingested free food water
Food_H2O_assoc_frac<- Food_H2O_assoc*ratio_frac_plant_material; 
# Moles of leaf water ingested
Food_H2O_assoc_nonfrac<- Food_H2O_assoc-Food_H2O_assoc_frac;
# Moles of stem water ingested

# Oxidase Water #
O2_resp_2<-Energy*0.00216*22.4*1000
# Total amount of O2 consumed 
resp_Atoms_O<-O2_resp_2*(5.38*10^19)
# Total atoms of oxygen consumed in 24 hours (see Whiteman et al. 2019)
Oxidase_H2O<-resp_Atoms_O/(6.023*10^23)
# Total moles of oxygen consumed in 24 hours (see Whiteman et al. 2019)

# Atmospheric O2 #
O2_resp<- Oxidase_H2O/2
# O2 respired (mol)

# Atmospheric H2O #
Air_ex<- O2_resp/Frac_used_O2/0.21*22.4;
# Amount of air fluxed through the lungs (L)
H2O_sat_content<-10^(0.686+0.027*(T-273.15))/760/22.4;
# Saturation concentration of H2O in air at ambient temperature (mole/L)
Air_H2O<- rh*H2O_sat_content*Air_ex;
# Amount of H2O from atmosphere (mole)

# Drinking H2O #
Water_turnover<- ((Mass*1000)*Perc_Water*Turnover_Fraction)/18
# Calculate amount of water turnover
H2O_in_nodrink<-Condensation_Rxn_H2O+Food_H2O_assoc_frac+Food_H2O_assoc_nonfrac+Oxidase_H2O+Air_H2O+Decarboxylation_Ox
# Establish all other water inputs except drinking
if (Water_turnover-H2O_in_nodrink > 0)
  Drinking<-(Water_turnover-H2O_in_nodrink) else     # Animal is drinking water
    (Drinking<-0)                                    # Animal does not drink water
# Determine the amount of drinking water
T_in<-Condensation_Rxn_H2O+Food_H2O_assoc_frac+Food_H2O_assoc_nonfrac+Oxidase_H2O+Air_H2O+Decarboxylation_Ox+Drinking
# Total inputs (mole)

# Assign final water turnover value based on whether all other water inputs exceeds estimated water turnover or not
if (H2O_in_nodrink > Water_turnover)
  H2O_turnover<- (H2O_in_nodrink) else                
    (H2O_turnover<-Water_turnover)                   

#    Output    #

Breath_H2O<-Air_ex*10^(0.686+0.027*Temp)/760/22.4;
# Exhaled water vapor from breathing (mole)
Oral_H2O_breath<-Breath_H2O/2;      
# Oral water loss through breathing                           
M_Fecal<-Food_consumed*(1-Digest);
# Dry fecal output (kg)
Fecal_O2<-0;                 
# Fecal O2 content (mole)
Fecal_H2O<-M_Fecal/(1-FecalH2O)*FecalH2O*55.56;
# Fecal H2O content calculated from dry fecal output (mole)
Food_urea<- 6*Food_consumed*Protein_content*Digest*Energy_effi;  
# Mole of Urea or Uric acid (mole)
Urea_O2<-Food_urea/2;        
# Urea/uric acid O2 content (mole)
Nasal_H2O<- Breath_H2O/4;    
# Nasal water vapor loss (mole)
Vapor_Pressure_Deficit<-((610.78*exp((7.04/(7.04+237.3)*17.2694))/1000)*(1-(0.61)))
#kilopascals; For calculation of transcuatenous water vapor loss (Schmidt-Nielsen 1969)
Vapor_Pressure_Deficit_mmHg<-Vapor_Pressure_Deficit*(760/101.325) 
#convert Vapor_Pressure_Deficit to mmHg
Skin_Vap<-Vapor_Pressure_Deficit_mmHg*Skin_evap_constant
# Multiply by taxa specific skin H2O evaporation constant
Perm<-(meeh_factor*Skin_Vap/1000)*24/18 
# Skin permeability factor calculated from Skin_Vap and meeh factor
Skin_H2O<- Perm*Mass^(2/3); 
# Transcutaneous water vapor loss (mole)
CO2<- O2_resp-Food_H2_O2-Urea_O2;               
# CO2 loss, corrected for the H2_O2 difference and Urea or uric acid (mole)
RQ<- CO2/O2_resp;             
# CO2/O2 ratio  
Oral_H2O<-Oral_H2O_breath     
# Total oral water loss
Urine<- H2O_turnover-Fecal_H2O-Nasal_H2O-Skin_H2O-Oral_H2O-CO2*2-Urea_O2*2;
# Urinary water loss (mole) is assumed to be remaining residual
T_out<- Fecal_H2O+Urea_O2*2+Urine+Nasal_H2O+Skin_H2O+Oral_H2O+CO2*2;
# Total oxygen output (mole)

# Input fractions
ffo<- Condensation_Rxn_H2O/T_in;             # Fraction of condensation water
ffrac<- Food_H2O_assoc_frac/T_in;            # Fraction of fractionating plant tissue
fnfrac<- Food_H2O_assoc_nonfrac/T_in;        # Fraction of non-fractionating plant tissue
fO2<- Oxidase_H2O/T_in;                      # Fraction of oxidase water
fvap<- Air_H2O/T_in;                         # Fraction of inhaled water vapor
fdw<- Drinking/T_in;                         # Fraction of drinking water
fCO_in<-Decarboxylation_Ox/T_in              # Fraction of decarboxylation oxygen
F_in<- ffo+fO2+fvap+fdw+ffrac+fnfrac+fCO_in;     # Total input flux (=1)

# Output fractions
ffecO<- Fecal_O2/T_out;                      # Fraction of fecal O2
ffec<- Fecal_H2O/T_out;                      # Fraction of fecal water loss
furO<- 2*Urea_O2/T_out;                      # Fraction of oxygen loss in urea/uric acid
fur<- Urine/T_out;                           # Fraction of urinary water loss
fnas<- Nasal_H2O/T_out;                      # Fraction of nasal water vapor loss
ftrans<- Skin_H2O/T_out;                     # Fraction of transcutaenous water vapor loss
foral<- Oral_H2O/T_out;                      # Fraction of oral water vapor loss
fCO2<- CO2*2/T_out;                          # Fraction of exhaled CO2
F_out<- ffecO+ffec+furO+fur+fnas+ftrans+foral+fCO2;
# Total output flux (=1)

# Fractionation factors and exponents (modified from Hu et al. 2023) #

# Fractionation Factors #

# Input #
# Fractionation factors between atmospheric H2O vapor and meteoric water (mw)
a_vap_mw_18<-1/alpha1(T);                            
# alpha_18O
l_vap_mw<-0.529;                                    
# fractionation exponent between water vapor and liquid water (Barkan & Luz 2005)
a_vap_mw_17<-exp(l_vap_mw*log(a_vap_mw_18));         
# alpha_17O

# Fractionation factors between drinking water and mw                            
a_dw_mw_18<-1;                                       
# alpha_18O
a_dw_mw_17<-1;                                       
# alpha_17O

# Fractionation factors for atmospheric O2
a_in_O2_18<- (1000+(d18O_atmO2-Z_value)/(1-Frac_used_O2))/(1000+d18O_atmO2);   
# alpha_18O
l_in_O2<- 0.5179;                                    
# Fractionation exponent between lung-absorbed O2 and atm.O2, (Luz & Barkan 2005)
a_in_O2_17<- exp(l_in_O2*log(a_in_O2_18));           
# alpha_17O

# Fractionation factors between leaf water and mw (for fractionating plant tissue)                          
a_lw_mw_18<-(1000+d18O_Leaf)/(1000+d18O_mw);       
# alpha_18O
l_lw_mw<- -0.0078*rh+0.5216;                         
# Fractionation exponent between leaf H2O and meoric H2O (Landais et al. 2006)
a_lw_mw_17<- exp(l_lw_mw*log(a_lw_mw_18));       
# alpha_17O

# Fractionation factors between stem water and mw (for non-fractionating plant tissue)                           
a_stw_mw_18<-1;                                      
# alpha_18O
a_stw_mw_17<-1;                                      
# alpha_17O

# Fractionation factors between food cellulose and mw 
a_cell_mw_18<-(1000+d18O_Cellulos)/(1000+d18O_mw);   
# alpha_18O
l_cell_mw<- 0.525;                                  
# Fractionation exponent bewteen cellulose and meteoric water (lower end of Pack et al. 2013)
a_cell_mw_17<- exp(l_cell_mw*log(a_cell_mw_18)); 
# alpha_17O

# Output # 
# Fractionation factors between fecal water and body water (bw) 
a_fec_bw_18<-1;                                      
# alpha_18O
a_fec_bw_17<-1;                                      
# alpha_17O

# Fractionation factors between urinary water and bw 
a_ur_bw_18<-1;                                       
# alpha_18O
a_ur_bw_17<-1;                                       
# alpha_17O

# Fractionation factors between nasal water vapor and bw
a_nas_bw_18<-1/alpha1((T+Body_T)/2);                 
# alpha_18O
l_nas_bw<-0.529;                                     
# Fractionation exponent bewteen nasal water vapor and body water (Barkan & Luz 2005)
a_nas_bw_17<-exp(l_nas_bw*log(a_nas_bw_18));         
# alpha_17O

# Fractionation factors between transcutaneously-lost water vapor and bw
a_trans_bw_18<-1/1.018;                              
# alpha_18O, Kohn (1996)
l_trans_bw<-0.5235;                                  
# Fractionation exponent between transcutaneous water vapor loss and body water (Hu et al. 2023)
a_trans_bw_17<-exp(l_trans_bw*log(a_trans_bw_18));   
# alpha_17O

# Fractionation factors between exhaled CO2 and bw
frac_a_CO2_bw_18<-alpha2(Body_T);                    
# alpha_18O, Fuction alpha2, O'Neil & Adami (1969)
l_CO2_bw<-0.5248;                                    
# Fractionation exponent between CO2 and H2O (Cao & Liu 2010)
frac_a_CO2_bw_17<-exp(l_CO2_bw*log(frac_a_CO2_bw_18));    
# alpha_17O

# Fractionation factors between exhaled H2O and bw
a_brea_bw_18<-1/alpha1(Body_T);                      
# alpha_18O
l_brea_bw<-0.529;                                    
# Fractionation exponent between breath water vapor and body water (Barkan & Luz 2005)
a_brea_bw_17<-exp(l_brea_bw*log(a_brea_bw_18));  
# alpha_17O

# Calculations #

# Air O2 composition
d17O_atmO2<-(exp((D17O_atmO2+lamda*log(d18O_atmO2/1000+1)*1000)/1000)-1)*1000;
# d17O of atmospheric O2 based on input information
R_18_atmO2<-((d18O_atmO2/1000)+1)*R_18_SMOW;
# Calculated R_18O of atmospheric O2
R_17_atmO2<-((d17O_atmO2/1000)+1)*R_17_SMOW;
# Calculated R_17O of atmospheric O2

# Meteoric H2O composition
d17O_mw<-(exp((D17O_mw+lamda*log(d18O_mw/1000+1)*1000)/1000)-1)*1000;
# d17O of meteoric water based on input information
R_18_mw<-((d18O_mw/1000)+1)*R_18_SMOW;
# Calculated R_18O of meteoric H2O
R_17_mw<-((d17O_mw/1000)+1)*R_17_SMOW;
# Calculated R_17O of meteoric H2O

# Mass Balance Equation #

# 18R #
# Input
R_18_in<-R_18_mw*(fvap*a_vap_mw_18+fdw*a_dw_mw_18+ffrac*a_lw_mw_18+fnfrac*a_stw_mw_18+ffo*a_cell_mw_18+(fCO_in*a_cell_mw_18));
# Oxygen_18 input term exclude atmospheric O2
R_18_inatm<-R_18_atmO2*fO2*a_in_O2_18;
# Oxygen_18 input term from atmospheric O2

# Output
r_18_out<-ffec*a_fec_bw_18+furO*a_ur_bw_18+fur*a_ur_bw_18+fnas*a_nas_bw_18+ftrans*a_trans_bw_18+foral*a_brea_bw_18+(fCO2*frac_a_CO2_bw_18);
# Oxygen_18 Output term  

# 17R #
# Input
R_17_in<-R_17_mw*(fvap*a_vap_mw_17+fdw*a_dw_mw_17+ffrac*a_lw_mw_17+fnfrac*a_stw_mw_17+ffo*a_cell_mw_17+(fCO_in*a_cell_mw_17));
# Oxygen_17 input term exclude atmospheric O2
R_17_inatm<-R_17_atmO2*fO2*a_in_O2_17;
# Oxygen_17 input term from atmospheric O2

# Output
r_17_out<-ffec*a_fec_bw_17+furO*a_ur_bw_17+fur*a_ur_bw_17+fnas*a_nas_bw_17+ftrans*a_trans_bw_17+foral*a_brea_bw_17+(fCO2*frac_a_CO2_bw_17);
# Oxygen_17 Output term 

#    Results   #

R_18_bw<-(R_18_in+R_18_inatm)/r_18_out;
# Oxygen_18 mass balance
R_17_bw<-(R_17_in+R_17_inatm)/r_17_out;
# Oxygen_17 mass balance

d18O_bw_i<-(R_18_bw/R_18_SMOW-1)*1000
# d18O of body water
d17O_bw_i<-(R_17_bw/R_17_SMOW-1)*1000
# d17O of body water
D17O_bw_i<-(log(d17O_bw_i/1000+1)-lamda*log(d18O_bw_i/1000+1))*1000
# D17O of body water

#Prints d17O, d18O, and E17O#
AnimalBodyWater<-function(R_18_bw,R_18_SMOW, R_17_bw, R_17_SMOW, d17O_bw_i, d18O_bw_i, lamda) {
  d18O_bw_i<-(R_18_bw/R_18_SMOW-1)*1000
  d17O_bw_i<-(R_17_bw/R_17_SMOW-1)*1000
  D17O_bw_i<-(log(d17O_bw_i/1000+1)-lamda*log(d18O_bw_i/1000+1))*1000
  print(paste("d18O_bw=",d18O_bw_i))
  print(paste("d17O_bw=",d17O_bw_i))
  print(paste("E17O=",D17O_bw_i))}

AnimalBodyWater(R_18_bw,R_18_SMOW, R_17_bw, R_17_SMOW, d17O_bw_i, d18O_bw_i, lamda)