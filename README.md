# MicrochipCriticalLiftOff

## MATLAB project for data analysis for the characterization of the critical lift-off of non-neutrally buoyant, square flat-plate microchip particles in microchannel flows
This code was used in the following paper:
Yeung, R., Sainz, C., Mandala, J., Brisk, P., Grover, W. H., and Rodgers, V. G. J. “Characterization of the critical lift-off of a single flat-plate microchip particle in straight rectangular microchannel flows”. International Journal of Multiphase Flow, 193: 105355 (2025). doi: https://doi.org/10.1016/j.ijmultiphaseflow.2025.105355

This code was run using MATLAB R2024b on Windows 10/11. File paths need to be modified for Unix-based platforms (macOS, Linux).

## Structure of the MATLAB Project
### Main script
**B2KCriticalLiftOff.m** - main script to perform the analysis.
### Functions
**B2KAnnotateModelParams** - function to output text annotation for model parameters with uncertainties to 1 significant figure.
**B2KFormatValuesToStrOneSigFigMaxUnc** - function to rounds values based on 1 significant figure of the maximum uncertainty for all values and formats as strings
**B2KTextAlignedToLineLogLogPlot** - function to align text annotation to power-law linear line in log-log scale
### Classes
**B2KAllParametersLiftOff** - class containing all equations for evaluating particulate flow properties including the area domain size, $HW/d^2$, and dimensionless numbers, $\mathrm{Re}$ and $\mathrm{Ar}$.
#### +B2KConstants
classes for storing particle and fluid properties
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