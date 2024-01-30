# Postprocess_ASPECT
A matlab package that postprocess data generated from ASPECT geodynamics code
https://github.com/geodynamics/aspect

Current Aspect version used working with 
https://github.com/Minerallo/aspect/tree/main_compressible_density

--------------Instructions--------------

Use Postprocess_v1_2.m to set your parameters for processing and visualization and then run it.
The corresponding folder should be automatically added to the path. 
Be sure to read the comments. 
If you have any questions, write to me at 
ponsm@gfz-potsdam.de 
or on git

thank you !

----------------log----------------------

## January 30, 2024 Postprocess_ASPECTv1.2

Version 1.2 contains a new set of postprocesses, the so-called geofeatures, which in combination with
Paraview python script will be able to track and quantify the number of subduction zones
associated subduction zones and plumes. 
The models also create a Robinson projection of the data extracted from Paraview. 
-The data can be extracted from paraview by using the Global_extract_series.py. 
-Additional work has been done to quantify the distribution of the age of the oceanic seafloor relative to its surface
surface area, as in Colice, 2012. 
The postprocess can now calculate the Feret diameter of subduction zones and the VRMS of continents.

