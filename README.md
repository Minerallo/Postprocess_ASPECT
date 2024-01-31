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

### January 30, 2024 - Postprocess_ASPECTv1.2

Version 1.2 introduces a new set of postprocesses known as "geofeatures." These, in combination with a Paraview Python script, enable tracking and quantification of subduction zones, associated subduction zones, and plumes. Additionally:

- The models now create a Robinson projection of the data extracted from Paraview.
- Data extraction from Paraview is facilitated using `Global_extract_series.py`.
- There's an added focus on quantifying the distribution of the age of the oceanic seafloor relative to its surface area (as in Colice, 2012).
- The postprocess can now calculate the Feret diameter of subduction zones and the VRMS of continents.
