# BaFTA: Bayesian Fertility Trajectory Analysis (Beta).

`BaSTA` is an R package for parametric Bayesian estimation of age-specific fertility for aggregated and individual level data.

## What's in BaFTA?

- Ability to analyze aggregated or individual level data for seasonal or continuous breeding with variable inter-birth intervals;
- Testing different models of age-specific fertility and compared their performance by means of DIC or posterior predictive check;
- Plots of traces to visually inspect convergence, of parameter posterior densities, of predicted age-specific fertility, and predictive plots for goodness of fit;

## How to install BaFTA?
To install BaFTA from GitHub, type the following lines of code on the R console:

```R
# Install and load 'devtools':
install.packages("devtools")
library(devtools)

# Install BaSTA2.0:
install_git("https://github.com/fercol/BaFTA", subdir = "pkg/")
```
