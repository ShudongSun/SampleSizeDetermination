# SampleSizeDetermination
A transparent simulator that generate different sample size data from pilot data based on multivariate guassian distribution. Then it will calculate the corresponding classification error/ARI/AMI and draw the plot which has the same trend as the true data. People can determine the sample size according the plot we draw.

## Installation
```r
if(!require(devtools)) install.packages("devtools")
library(devtools)

devtools::install_github("ShudongSun/SampleSizeDetermination")
```
