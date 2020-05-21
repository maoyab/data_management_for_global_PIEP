# Data Management For Global PIEP


**This repository contains code associated with data and result processing in:**

Bassiouni, M., S.P. Good, C.J. Still, and C.W. Higgins (2020), Plant water uptake thresholds inferred from satellite soil moisture. Geophysical Research Letters. https://doi.org/10.1029/2020GL087077 


Code to estimate ecohydrological parameters that best fit empirical soil saturation probability functions see:

*Probabilistic Inference of Ecohydrological Parameters (PIEP)*: http://doi.org/10.5281/zenodo.1257718 or https://github.com/maoyab/PIEP


Global dataset of pcohydrological parameters resulting from this study is available at https://doi.org/10.5281/zenodo.3351622


**File descriptions:**

*combine_L3_soil_moisture_data.py* - Combines SMAP L3 soil moisture data 

*combine_L4_climate_data.py* - Combines selected variables in SMAP L4 data

*compile_attributes.py* - Compiles global attributes for PIEP

*make_SMAP_PIEP_iterable.py* - Creates iterables for PIEP

*combine_results_0.py* - Combines PIEP results

*process_results.py* - Summarizes global PIEP results and creates .nc dataset of global ecoydrological parameters


**Datasets used**:

SMAP L3 soil moisture data: https://doi.org/10.5067/ZX7YX2Y2LHEB

SMAP L4 Geophysical data: https://doi.org/10.5067/KPJNN2GI1DQR

Global soil hydraulic parameters: https://doi.pangaea.de/10.1594/PANGAEA

IGBP landcover classification: https://doi.org/10.5067/KGLC3UH4TMAQ

