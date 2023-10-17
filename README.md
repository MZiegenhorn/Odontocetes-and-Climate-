# Odontocetes-and-Climate-
Code accompanying Ziegenhorn et al. (2023), "Odontocete detections are linked to oceanic conditions in the Hawaiian Archipelago", Nature Communications Earth and Environment

Please read on for a description of all files. 

1. allClimate.mat: a MATLAB file containing values of ENSO, PDO, and NPGO during the years 2009-2019 (monthly). For direct access to these datasets, please refer to the 'Data Availability' section of the related manuscript (see above).

2. climate_correlationAnalysis.mat: MATLAB code used to produce Supplemental Figs. 1-4 in the related manuscript.

3. final_timeseries_anom.mat: MATLAB code used to produce Figs. 2 and 4 in the related manuscript.

4. gmt_HAWAII.mat: MATLAB code (using Generic Mapping Tools) to produce Fig. 1 in the related manuscript.

5. Hawaii_allclasses_counts_finalData.mat; Manawai_allclasses_counts_finalData.mat: MATLAB files containing a table of odontocete detections and corresponding values of oceanographic variables (sea surface temperature, sea surface salinity, sea surface height, ENSO, PDO, NPGO) at the specified site. Daily values are given. More information can be found in the associated manuscript.

Note that for these files as well as their associated .csv files (which are used in the master modelling script), the following classes of odontocetes are considered: 

class1- Blainville's beaked whale, Mesoplodon densirostris

class2- Cuvier's beaked whale, Ziphius cavirostris

class3- false killer whale, Pseudorca crassidens

class4- Rough-toothed dolphin, Steno bredanensis

class5_6- short-finned pilot whale, Globicephala Macrorhyncus

class7_8- stenellid dolphins, mix of Stenella longirostris, Stenella attenuata, and Stenella coruleoalba

class9- common bottlenose dolphin/ melon-headed whale, Tursiops truncatus/Pepnocephala electra

class10- Kogia spp (pgymy and dwarf sperm whale), Kogia breviceps and Kogia sima


6. Hawaii_counts_finalData.csv; Manawai_counts_finalData.csv - these files contain the same data as in 5, but in CSV format.

7. masterModellingScript.R - An R script used for creating all models in the related manuscript.



Any questions about this code and dataset should be directed to Morgan A. Ziegenhorn, maziegenhorn36@gmail.com. Please feel free to reach out!






