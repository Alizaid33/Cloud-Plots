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




## Requirements

The analysis was developed in R.

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

The exact package versions used in the analysis can be inspected using:

```r
sessionInfo()
```
after installing the required packages.

## Reproducibility

The analysis uses the following random seed:

```r
set.seed(123)
```

This ensures reproducibility of any random-number generation used by the analysis.

The random seed should be set before running the analysis.

## General Analysis Workflow

The main function in `cloud_plot.R` is:

```r
analyze_and_visualize_spherical_data(DATA)
```
The function expects a data frame containing two columns:

```
theta
phi
```
where the angular measurements are expressed in radians.

The analysis proceeds as follows.

#### Step 1: Convert spherical coordinates to Cartesian coordinates

For observations with spherical coordinates (theta, phi), the corresponding Cartesian coordinates are calculated as:

```r
x = sin(theta) * cos(phi)
y = sin(theta) * sin(phi)
z = cos(theta)
```
This maps every observation to the unit sphere:

```r
x^2 + y^2 + z^2 = 1.
```

#### Step 2: Estimate the directional parameters

The code calculates:

- the mean direction;
- the mean resultant length; and
- the estimated concentration parameter.

The concentration parameter is estimated from the mean resultant length.

#### Step 3: Calculate geodesic distances

For two unit vectors p1 and p2, the geodesic distance is calculated as:

```r
acos(sum(p1 * p2))
```

with numerical protection against small floating-point errors.

#### Step 4: Construct the Cloud Plot

The Cloud Plot is constructed using the spherical median and the distribution of geodesic distances from the observations to the spherical median.

The method calculates distance quantiles and constructs an upper boundary for identifying observations that lie unusually far from the central directional structure.

#### Step 5: Identify potential outliers

Observations whose geodesic distance from the spherical median exceeds the upper Cloud Plot boundary are identified as potential outliers.

#### Step 6: Visualize the data

The observations are displayed on the unit sphere using the rgl package.

The visualization includes:

the unit sphere;
spherical coordinate grid lines;
the observed spherical data points;
the spherical median;
observations identified as potential outliers; and
the Cloud Plot region.



## Application 1: Eye Data

The eye data are provided in the repository as:

```
eye_data.xlsx
```

The data contain the spherical coordinates required by the analysis.
### Running the eye-data analysis

Open: `cloud_plot.R`.

Set:

```r
dataset <- "eye"
```
The relevant part of the script is:

```r
if (dataset == "eye") {

  my_data_df <- read_excel("eye_data.xlsx")

}
```
Then run:

```r
analyze_and_visualize_spherical_data(my_data_df)
```

3. Run the complete script.

The script performs:

- loading the `eye_data.xlsx` spherical coordinates,
- read the spherical coordinates;
- convert the observations to Cartesian coordinates;
- estimate the directional parameters;
- calculate geodesic distances;
- construct the Cloud Plot;
- identify potential outliers; and
- produce the three-dimensional visualization.


### Application 2: Gait Data

The second application uses the `gait` dataset from the `fda` R package.

The dataset contains angular measurements describing the hip and knee positions of 39 children over a complete gait cycle.

For each child, the gait cycle is represented at 20 equally spaced time points.

The two angular measurements are transformed into points on the unit sphere using the spherical-coordinate transformation described above.

#### Full gait-trajectory analysis

The gait data are stored as a three-dimensional array containing:

```
subject × time point × angular variable
```

The complete gait trajectory for a subject can be extracted as:

```r
subject_data <- gait[subject_index, , ]
```
where:

```r
subject_index
```
specifies the child.

For example:

```r
subject_index <- 1

subject_data <- gait[subject_index, , ]

my_data_df <- data.frame(
  theta = subject_data[, 1] * pi / 180,
  phi   = subject_data[, 2] * pi / 180
)
```
This produces 20 spherical observations corresponding to the 20 equally spaced time points of the complete gait cycle for the selected child.

The resulting data frame has the form:
```
theta     phi
------    ------
t1        t1
t2        t2
...       ...
t20       t20
```
The Cloud Plot analysis is then applied to the complete trajectory:

```r
analyze_and_visualize_spherical_data(my_data_df)
```

#### Running the Gait Example

To run the gait example, set:

```r
dataset <- "gait"
```
The corresponding code is:
```r
library(fda)

subject_index <- 1

subject_data <- gait[subject_index, , ]

my_data_df <- data.frame(
  theta = subject_data[, 1] * pi / 180,
  phi   = subject_data[, 2] * pi / 180
)

analyze_and_visualize_spherical_data(my_data_df)
```

#### The analysis then:

- loads the gait dataset from the fda package;
- selects one child's complete gait trajectory;
- extracts the 20 equally spaced gait-cycle observations;
- converts the hip and knee angles from degrees to radians;
- maps the angular observations to the unit sphere;
- estimates the directional parameters;
- calculates geodesic distances;
- constructs the Cloud Plot; and
- produces the three-dimensional visualization.


#### Changing the Gait Subject

The subject analyzed can be changed by modifying:

```r
subject_index <- 1
```
For example:

```r
subject_index <- 10
```
will analyze the complete 20-point gait trajectory for the tenth child.

The same analysis workflow is then applied:
```r
subject_data <- gait[subject_index, , ]

my_data_df <- data.frame(
  theta = subject_data[, 1] * pi / 180,
  phi   = subject_data[, 2] * pi / 180
)
```

### Complete Minimal Example

The following example reproduces the gait analysis:

```r
library(fda)

# Select one child
subject_index <- 1

# Extract the complete gait trajectory
subject_data <- gait[subject_index, , ]

# Convert angular measurements from degrees to radians
my_data_df <- data.frame(
  theta = subject_data[, 1] * pi / 180,
  phi   = subject_data[, 2] * pi / 180
)

# Run the Cloud Plot analysis
analyze_and_visualize_spherical_data(my_data_df)
```


For the eye dataset:

```r
library(readxl)

# Load the data
my_data_df <- read_excel("eye_data.xlsx")

# Run the Cloud Plot analysis
analyze_and_visualize_spherical_data(my_data_df)
```


# Reproducing the Results

To reproduce the analyses:

- Install R.
- Install the required packages.
- Clone or download this repository.
- Open cloud_plot.R in R or RStudio.
- Select either:
```r
dataset <- "eye"
```
or:
```r
dataset <- "gait"
```
- Run the complete script.

The code will perform the full analysis and generate the corresponding three-dimensional visualization.

## Notes on the Gait Application

The gait application uses the complete trajectory for one child at a time.

The 20 observations correspond to 20 equally spaced time points over one complete gait cycle. Each pair of hip and knee angular measurements is mapped to a point on the unit sphere.

Consequently, the analysis is performed on the spherical representation of the complete gait trajectory rather than on a single cross-sectional time point.

## Citation

If you use this code, please cite the associated research article.

The gait dataset is obtained from the fda R package.

## Contact

For questions concerning the implementation or the Cloud Plot methodology, please refer to the associated research article or contact the authors.

