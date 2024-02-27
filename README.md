# Postprocess_ASPECT

**Description:**
A Matlab package designed for postprocessing data generated from the [ASPECT geodynamics code](https://github.com/geodynamics/aspect). This version is specifically tailored to work with [Minerallo's compressible density branch](https://github.com/Minerallo/aspect/tree/main_compressible_density) of ASPECT.

## Instructions

1. Use `Postprocess_v1_2.m` to set your parameters for processing and visualization.
2. Run the script.
3. The corresponding folder should be automatically added to the path.
4. Read the comments for additional guidance.
5. For any questions, feel free to reach out at [ponsm@gfz-potsdam.de](mailto:ponsm@gfz-potsdam.de) or on [GitHub](https://github.com/geodynamics/aspect).

Thank you!

## Change Log

### February 27, 2024 - Postprocess_ASPECTv1.3

Version 1.3 introduce the automatization of Pvbatch with the matlab postprocess. To do so the user may : 

1) Set Python_paraview/build_autopostprocess_files, Postprocess_v1_2_auto.m (will be renameed later) and Postprocess.sh with the correct paths and make sure everything is correct.
2) Without running them, open your terminal and go to the Python_paraview repository, use chmod+x Postprocess.sh and then ./Postprocess.sh. 
This will change the model name in the Python paraview script and in the Matlab postprocess.
3) A new directory has been created. Change to this directory with your terminal and use chmod+x Postprocess_auto.sh and then ./Postprocess_auto.sh
This will copy the files to the server to run pvbatch and extract the model data. 
3b) Alternatively, you can simply restart without running pvbatch. Then you can use chmod+x Postprocess_auto_restart.sh and then ./Postprocess_auto_restart.sh.
This will check if new files should be downloaded and run the Matlab postprocess to create the figures.

### January 30, 2024 - Postprocess_ASPECTv1.2

Version 1.2 introduces a new set of postprocesses known as "geofeatures." These, in combination with a Paraview Python script, enable tracking and quantification of subduction zones, connected subduction zones, and plumes. Additionally:

- The models now create a Robinson projection of the data extracted from Paraview.
- Data extraction from Paraview is facilitated using `Global_extract_series.py`.
- There's an added focus on quantifying the distribution of the age of the oceanic seafloor relative to its surface area (as in Coltice, 2012).
- The postprocess can now calculate the Feret diameter of subduction zones and the VRMS of continents.
