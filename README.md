# Cloud Plots

## Overview

This repository contains an R implementation of a Cloud Plot method for spherical data. The method provides visualization of spherical observations and identifies potential outliers.

## Repository contents

```
cloud_plot.R
eye_data.xlsx
README.md
```

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
