# MicrochipCriticalLiftOff

## MATLAB project for data analysis for the characterization of the critical lift-off of non-neutrally buoyant, square flat-plate microchip particles in microchannel flows
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17082194.svg)](https://doi.org/10.5281/zenodo.17082194)

This code was used in the following paper:
Yeung, R., Sainz, C., Mandala, J., Brisk, P., Grover, W. H., and Rodgers, V. G. J. “Characterization of the critical lift-off of a single flat-plate microchip particle in straight rectangular microchannel flows”. International Journal of Multiphase Flow, 193: 105355 (2025). doi: https://doi.org/10.1016/j.ijmultiphaseflow.2025.105355

This code was run using MATLAB R2024b (using the Java-based desktop and graphics system) on Windows 10/11. File paths need to be modified for Unix-based platforms (macOS, Linux). MATLAB versions R2025a and above, which uses a new environment based on HTML and JavaScript, may cause unintended graphical issues.

## How to setup the MATLAB Project
1. Open MATLAB (current version tested for MATLAB R2024b).
2. Navigate to the 'HOME' tab, and then click 'New' > 'Project', 'From Git'.
    1. Ensure 'Source control tool' uses Git
    2. For the 'Repository path', enter: "https://github.com/raykyeung/MicrochipCriticalLiftOff"
    3. For the 'Sandbox", select a suitable location for the local repository (ex. "C:\Users\USERNAME\LocalGitRepo\MicrochipCriticalLiftOff")
    4. Click 'Retrieve' and allow MATLAB to pull all files. Project dependencies and tracking are handled by the .prj file and 'resources' folder.

## Structure of the MATLAB Project
### Main script
**B2KCriticalLiftOff.m** - Main script to perform the analysis. Briefly, the script performs the following:
1. Defines parameters chosen for independent and dependent variables
2. Inputs all experimentally observed critical lift-off flow rates for all channel and solvent configurations
3. Inputs comparison data
4. Imports and parses data from computational simulations using COMSOL Multiphysics
5. Defines all parameter properties with their associated strings and dimensions
6. Evaluate all parameter properties for all configurations
7. Perform analysis including non-linear regression of parameters for the corrected generalized lift-off model using a remodified Archimedes number
8. Export all figures
9. Save log file containing dependency information and required program files

### Functions
**B2KAnnotateModelParams** - Function to output text annotation for model parameters with uncertainties to 1 significant figure\
**B2KFormatValuesToStrOneSigFigMaxUnc** - Function to rounds values based on 1 significant figure of the maximum uncertainty for all values and formats as strings\
**B2KTextAlignedToLineLogLogPlot** - Function to align text annotation to power-law linear line in log-log scale
### Classes
**B2KAllParametersLiftOff** - Class containing all equations for evaluating particulate flow properties including the area domain size, $HW/d^2$, and dimensionless numbers, $\mathrm{Re}$ and $\mathrm{Ar}$.
#### +B2KConstants
Classes for storing particle and fluid properties:
* @Acetonitrile
* @Ethanol
* @Glycerol
* @Isopropanol
* @Methanol
* @pChip
* @Water

## Copyright
Copyright (c) 2025 Raymond K. Yeung

This product includes software developed at the University of California, Riverside (https://www.ucr.edu/) for NSF CPS Grant Award #1740052 (https://www.nsf.gov/awardsearch/showAward?AWD_ID=1740052).

This product is licensed under the terms of the Apache 2.0 license.

The project contains the following code and derivatives of 3rd-party resources:
* Function suite, DataViz, Version 3.2.3 by Povilas Karvelis from Github (https://github.com/povilaskarvelis/DataViz) under MIT license, for plotting bar graphs
* Function suite, MATLABTools, Version 1.2.0.0 by Tyler Voskuilen from MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/46357-uncertainty-propagation-class-uc?s_tid=srchtitle) under BSD-3 License, for uncertainty propagation