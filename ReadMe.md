# Information on files withing /Manuscript1

## ftms_process.r
R script for processing ICMB Ocean output and calculating key metrics
Contains some figures
### Related data
#### hotw_edata.csv
The matrix of sample names and relative intensities of masses
#### hotw_emeta.csv
Meta data of molecular fomulae found in samples
#### hotw_fdata.csv
Meta data on samples including factors like sampling trip and burned or unburned

## chem.r
R script with main statistical analysis, linear models, and graphing for manuscript visuals

## extra_analysis.r
R script with extra statistical analysis to describe data broadly and other figures

## transformations.r
R script for calculating putative biochemical transformations of sample DOM
Needs the following .rds as a reference list
  - Putative_Trannsform_Matrix_List_i.rds

## Datasets on DOM samples from ICBM Ocean
### ftms3.csv
### mfs.csv

## Datasets on general site data
Note: there are multiple because some version control over the two years of collecting data resulted in different copies of .csv files. Scripts are updated to use the apropriate file for the apropriate analysis and merge files when needed.

### water.csv
Includes the following variables:
#### Sample information
- sample_id : unique identifier for each sample
- watershed : unique watershed ID
- duplicate : yes/no value
- treatment : control or burn
- trip : month of trip as text
- date : date in mm/dd/yyyy format
- burn_mean : mean watershed areas of the Red Lake 2021 fire disturbance extracted from Google Earth Engine with the google earth engine script created by Holsinger et al, 2021. Improved fire severity mapping in the North American boreal forest using a hybrid composite method. - Original script : https://code.earthengine.google.com/9d63e1928205e238a5d49d6d9af1d1a7
- burn_sd : standard deviation of  watershed areas of the Red Lake 2021 fire disturbance extracted from Google Earth Engine with the google earth engine script created by Holsinger et al, 2021. Improved fire severity mapping in the North American boreal forest using a hybrid composite method. - Original script : https://code.earthengine.google.com/9d63e1928205e238a5d49d6d9af1d1a7
- year : year of collection
- area : watershed area in ha
- slope : ERIN TO FILL
- prop_wet_area  : ERIN TO FILL
- hw_burn_mean_haifls : stream Hydroweighted watershed areas of the Red Lake 2021 fire disturbance extracted from Google Earth Engine with the google earth engine script created by Holsinger et al, 2021. Improved fire severity mapping in the North American boreal forest using a hybrid composite method. - Original script : https://code.earthengine.google.com/9d63e1928205e238a5d49d6d9af1d1a7
#### Watershed Ecology Team variables
##### TOC-L
Data collected from the TOC-L in the Watershed Ecology lab at GLFC
- IC : IC  in mg/L. DL 0.5
- TOC.mgL : Dissolved Organic Carbon, filtered at 0.45 um in mg/L. DL 0.4
- TN.mgL : Dissolved Nitrogen, filtered at 0.45 um in mg/L. DL 0.2
##### DOM spectral characteristics collected on the Aqualog
Metrics were calculated with the staRdom package
- fi : Fluorescence index as in McKnight et al. (2001)
- hix : Humification index as in McKnight et al. (2001)
- mhix : Modified Humification index as in Ohno (2002)
- bix : Biological Fluorescence index as in Hueget et al. (2009)
- a254 : Absorbance coefficent at 254nm (m^-1)
- a300 : Absorbance coefficient at 300nm (m^-1)
- E2_E3 : Absorbance ratio 250nm/365 nm
- S275_295 : Ratio used for calculation of SR
- S300_400 : Ratio used for calculation of SR
- SR : Spectral Slope Ratio (S275_295:S350_400) as in Helms et al. (2008)
- Comp.1 to Comp.7 : Components from the internal Watershed Ecology Parafac model.
        Components relate to: •	C1: Humic A+M
        + C2: Humic A+C
        + C3: Humic
        + C4: Humic A
        + C5: Protein (TRP)
        + C6: Protein (TYR)
        + C7: Microbial Humic
- SUVA: the absorption of light at 254 nm per unit of carbon, has been shown to be a useful proxy for DOM aromatic content
##### Variables collected in the field using a YSI
- temp : temperature in Celcius
- do2_mg_l : dissolved oxygen in mg/L
- cond_us_cm : conductivity in uS/cm
- ph : pH
- orp_mv : oxidation reduction potential in mV
#### Variables from the GLFC water lab
- alk.meql: total alkalinity analyzed by auto-titration
- Cl.mg.L: total Chlorine as analyzed by auto-titration. DL 0.01
- SO4.mg.L : Total Sulphate as analyzed on an AA. DL 0.2
- NO2.NO3.mgL: dissolved nitrate and nitrite in mg/L analysed on an AA. DL 0.04
- NH4.mg.L: dissolved NH4 in mg/L analysed on an AA. DL 0.01
- DOC.mgL : Total dissolved organic carbon in mg/L analysed on an AA. DL 0.4
- DIC.mgL: Total dissolved inorganic Carbon in mg/L analysed on an AA. DL 0.5
- waterlab_TN.mgL: Total nitrogen in mg/L analysed on an AA. DL 0.2
- TP_testmark.mgL or TP.mgL_testmark: Total Phosphorus as analysed by testmark labs. DL 0.002
#### Variables from FTIR analysis
- GFE : ERIN TO FILL
- NOSC : ERIN TO FILL
- AI_mod : ERIN TO FILL
- DBE_1 : ERIN TO FILL
- h.c.wtavg : ERIN TO FILL
- o.c.wtavg : ERIN TO FILL
- transformations : ERIN TO FILL
- bc : ERIN TO FILL
  
### watersheds_characterized.csv
### watersheds_characterized1.csv 

