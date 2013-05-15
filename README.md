StochasticPrecipitation
=======================

R code for generating multi-site stochastic precipitation from daily observations using the modified Wilks approach of Mhanna and Bauwens (2012).

In water resource management, climate change, hydrology and related disciplines long time series of precipitation/rainfall data is required. Since historical records are relatively short, typically 50 years or less, mathematical/statistical models are used to generate synthetic data of required length. This data generation process requires that the spatio-temporal statistics are preserved between the observations and the synthetic data. There are a number of procedures in the literature, largely beginning in the late 1990s with the work of Wilks (see references below). However, the code required to generate the synthetic data is typically not available. Only recently, two R packages have been made available (see links below).

This page provides the R code for the computationally and conceptually simple implementation of the modified Wilks approach of Mhanna and Bauwens (2012).

### Code
* "mod_wilks_lib.R" - is the library of R functions used by the code
* "mod_wilks_code.R" - this is the file you need to run
* "mod_wilks_plot.R" - for summary plots

### Input Data
The code on this page generates 1000 years of synthetic precipitation data for 15 stations (in the San Francisco Bay Area), only for the month of January. The input data is provided as the file "janData.txt". The format of this file is:
* year, day of month, rain @ stn 1, rain @ stn 2, ..., rain @ stn 15
* year, day of month, rain @ stn 1, rain @ stn 2, ..., rain @ stn 15
* ...

This code should work for any dail data set, provided the input format is consistent with the above specified format.


References
------------
* Mhanna and Bauwens (2012), A Stochastic space-time model for the generation of daily rainfall in the Gaza Strip, International Journal of Climatology, 32, 1098-1112.
* Brissette, Khalili and Leconte (2007), Efficient stochastic generation of multi-site synthetic precipitation data, Journal of Hydrology, 345, 121-133.

Links to R-based weather generators
------------
* HyetosR package from [University of Athens](http://itia.ntua.gr/en/docinfo/1200/)
* R Multi-site Auto-regressive Weather GENerator (RMAWGEN) package, [link](http://cran.r-project.org/web/packages/RMAWGEN/)
* Carey-Smith, Thomson and Sansom, 2012, [link](http://pagesperso.univ-brest.fr/~ailliot/roscoff/Carey.pdf)
* [Workshop on Stochastic Weather Generators](http://pagesperso.univ-brest.fr/~ailliot/SWGEN_workshop.html), 5/29-6/1, 2012, France

