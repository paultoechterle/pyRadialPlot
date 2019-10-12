# pyRadialPlot

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rbeucher/pyRadialplot/master)

A Python Package for Fission Track Radial Plots.

Romain Beucher, The University of Melbourne and The Australian National University

Inline-style: 
![RadialPlot](https://github.com/rbeucher/pyRadialPlot/blob/master/radialplot.png)


## What is a Radial plot?

The radial plot is a graphical method for displaying and comparing observations that have different precision. Invented by Rex Galbraith in 1988 it is commonly used in geochronology but as also a wide range of applications in business analytics or medical research. The observations are standardised and plotted against precision, with the precision defined as the reciprocal of the standard error. The original observations are given by slopes of lines through the origin and can be read using a circular scale. Radial plots provide a visual representation that can help to assess whether the estimates agree with a common value. They can also be used to identify outliers or groups of estimates differing in a systematic way because of some underlying factor or mixture of populations.    

## Why is it useful ?

In experimental sciences it is common to have measurements with different precision. This can arise from natural variations or from the experimental procedure. Geochronological methods such as fission track, $^{40}\text{Ar}/^{39}\text{Ar}$, U-Pb and Optically Stimulated Luminescence (OSL) dating produce age estimates and associated errors for each of several grains. The radial plot can be used to display and compare the age estimates and see how they agree or differ within standard statistical variation.

Another application of radial plot is in *meta-analysis*. Radial plots can be used to compare treatments effects from different clinical studies where the precision of the studies varies.

Learn more about Radial Plot here:
(Wikipedia Article on Galbraith Plots)[https://en.wikipedia.org/wiki/Galbraith_plot]

## Acknowledgments

Inspiration and design ideas have been borrowed from RadialPlotter developed by Peter Vermeesh at UCL.

## References

Campbell, I. H., Reiners, P. W., Allen, C. M., Nicolescu, S., Upadhyay, R., 2005: He-Pb double dating of detrital zircons from the Ganges and Indus Rivers: Implication for quantifying sediment recycling and provenance studies, Earth and Planetary Science Letters 237, 402-432.

Galbraith, R. F., 1988: Graphical display of estimates having differing standard errors, Technometrics, 30, 271-281.

Galbraith,R. F., 1990: The radial plot: graphical assessment of spread in ages, Nuclear Tracks and Radiation Measurements, 17, 207-214.

Galbraith, R. F. and Green, P. F., 1990: Estimating the component ages in a finite mixture, Nuclear Tracks and Radiation Measurements, 17, 197-206.

Galbraith, R. F., Roberts, R. G., Laslett, G. M., Yoshida, H. and Olley, J. M., 1999: Optical dating of single and multiple grains of quartz from Jinmium rock shelter, northern Australia: Part I, experimental design and statistical models, Archaeometry, 41, 339-364. 

Vermeesch, P., 2009, RadialPlotter: a Java application for fission track, luminescence and other radial plots, Radiation Measurements, 44, 4, 409-410
