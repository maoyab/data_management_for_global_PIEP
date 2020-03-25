# Data Management For Global PIEP
This repository contains code associated with data and result processing in:

Bassiouni, M., S.P. Good, C.J. Still, and C.W. Higgins (2020), Plant water uptake thresholds inferred from satellite soil moisture. Geophysical Research Letters. https://doi.org/10.1029/2020GL087077 

Code for PIEP is available at https://doi.org/10.5281/zenodo.1257718 or https://github.com/maoyab/PIEP

Global Dataset of Ecohydrological Parameters resulting from this study is available at https://doi.org/10.5281/zenodo.3351622


Datasets used:

combine_L3_soil_moisture_data.py - Combine SMAP L3 soil moisture data (https://doi.org/10.5067/ZX7YX2Y2LHEB)

combine_L4_climate_data.py - Combine selected variables in SMAP L4 data (https://doi.org/10.5067/KPJNN2GI1DQR)

global soil hydraulic parameters (https://doi.pangaea.de/10.1594/PANGAEA)

IGBP landcover classification (https://doi.org/10.5067/KGLC3UH4TMAQ.)


Files:

compile_attributes.py - Compile global attributes for PIEP including

make_SMAP_PIEP_iterable.py - create iterables for PIEP

combine_results_0.py - combine PIEP results

process_results.py - summarize global PIEP results

