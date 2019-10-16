# Damage Allowance: Codes and data for Rasmussen et al. (2019)

README file last updated by DJ Rasmussen, dmr2-at-princeton-dot-edu, Tue Oct 15 14:55:29 PDT 2019

## Citation

This code is intended to accompany the results of

    D. J. Rasmussen, M. K. Buchanan, R. E. Kopp, and M. Oppenheimer (2019). A flood damage allowance framework for coastal protection with deep uncertainty in sea-level rise. arXiv:1908.02844

Please cite this work when using any results generated with this code.

## Overview

The R program and associated libraries are required to run the codes. In the included directories are:


###Data Sets
`manhattan_property_area_elevation.csv` 

* Manhattan property value and area summed at discrete elevations and cumulatively summed. Elevation is given in meters, area is given in square miles, and property value is given in billions of 2017 US$. Property value is building value only (excludes land).

`GPDfits_global_tidegauge_uhawaii.tsv`

* Generalized Pareto distribution parameters for estimating extreme sea level frequencies at a global network of tide gauges generated using tide gauge data from the University of Hawaii tide gauge data set. For codes that generated these parameters, see: <https://github.com/dmr2/hawaiiSL_process> and <https://github.com/dmr2/GPDfit>

`slr_du_rasmussen19.tgz`

* Sea level projections with multiple assumptions of future AIS melt behavior (for more details, see Rasmussen et al., 2019). Projections are given for the global mean, the Battery tide gauge in New York, and total AIS melt contributions to the global mean. Projections include local sea-level rise Monte Carlo samples needed for running `damage_allowance.R`. Projections for other locations could be generated using the MATLAB program LocalizeSL <https://github.com/bobkopp/LocalizeSL>


###Codes

* **LocalizeSL_DU** A version of the MATLAB program LocalizeSL that also explores deep uncertainty in the contribution of Antarctic ice melt to sea level rise if given two core files from ProjectSL <https://github.com/bobkopp/LocalizeSL/tree/master/scripts/Rasmussen2019>

* **damage_allowance.R** Codes to calculate flood damage allowances based on assumptions of Antarctic ice melt (see Rasmussen et al., 2019 for the exact parameters explored). Local sea-level rise projections must be generated before running (with the exception of New York City, which have been provided in `data/slr_du_rasmussen19.tgz`

* **routines** Various routines to support running `damage_allowance.R` 


----

    Copyright (C) 2019 by D.J. Rasmussen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.