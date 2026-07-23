# Cloud Plot for Spherical Data

## Overview

This repository contains the R code used to construct a Cloud Plot for spherical data and to identify potential outlying observations using a geodesic-distance-based procedure.

The repository includes:

1. An application to the eye dataset supplied with this repository.
2. An application to the gait dataset available through the `fda` R package.

The code implements the complete workflow used in the applications:

1. loading the data;
2. converting spherical coordinates to Cartesian coordinates on the unit sphere;
3. estimating the mean direction and concentration parameter;
4. calculating geodesic distances;
5. constructing the Cloud Plot;
6. identifying potential outlying observations; and
7. producing a three-dimensional visualization of the observations and the resulting Cloud Plot region.

## Repository contents

```
cloud_plot.R
eye_data.xlsx
README.md
```


`cloud_plot.R`

Main R script containing:

- data loading;
- spherical-to-Cartesian coordinate transformation;
- estimation of the concentration parameter;
- geodesic-distance calculations;
- Cloud Plot construction;
- outlier identification; and
- three-dimensional visualization.

`eye_data.xlsx`

The eye-movement dataset used in the first application.

`README.md`

This file provides the complete workflow for reproducing the analyses.




## Required R packages

The following packages are required:

```r
library(sphunif)
library(rgl)
library(geometry)
library(Directional)
library(readxl)
library(fda)
```

Install missing packages using:

```r
install.packages(c(
"sphunif",
"rgl",
"geometry",
"Directional",
"readxl",
"fda"
))
```

## Reproducibility

The analysis uses the following random seed:

```r
set.seed(123)
```

The R version and package versions can be checked using:

```r
sessionInfo()
```

## Analysis workflow

### Example 1: Eye dataset

The eye dataset is provided as:

```
eye_data.xlsx
```

To run the eye analysis:

1. Open `cloud_plot.R`.
2. Set:

```r
dataset <- "eye"
```

3. Run the complete script.

The script performs:

- loading the eye spherical coordinates,
- conversion from spherical coordinates to Cartesian coordinates,
- estimation of the mean direction,
- estimation of the concentration parameter,
- Cloud Plot outlier detection,
- spherical visualization,


### Example 2: Gait dataset

To analyze the gait dataset:

```r
dataset <- "gait"
```

The gait data are obtained from the `fda` R package.

The same workflow is applied:

- coordinate conversion,
- parameter estimation,
- outlier detection,
- Cloud Plot visualization,


## How to run

Clone this repository and open `cloud_plot.R` in RStudio.

Run the script after installing the required packages.
