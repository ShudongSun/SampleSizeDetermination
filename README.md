# SampleSizeDetermination
A transparent simulator that generate different sample size data from pilot data based on multivariate guassian distribution. Then it will calculate the corresponding classification error/ARI/AMI and draw the plot which has the same trend as the true data. People can determine the sample size according the plot we draw.

## Installation
```r
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("ShudongSun/SampleSizeDetermination")
```
If dependency R packages are updated when installing the `devtools` R package, try closing and restarting the R session before proceeding.

## Tutorials
This [tutorial](https://github.com/ShudongSun/SampleSizeDetermination) introduces how to use the SSD R package.
