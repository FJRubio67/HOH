# On Harmonic Oscillator Hazard Functions

## Overview

This repository contains the R code and data used to reproduce the results in:

> Christen, J.A. and Rubio, F.J. (2025). On harmonic oscillator hazard functions.
> *Statistics & Probability Letters* **217**: 110304.
> [DOI: 10.1016/j.spl.2024.110304](https://doi.org/10.1016/j.spl.2024.110304) |
> [Preprint (arXiv)](https://arxiv.org/abs/2408.15964)

The **harmonic oscillator hazard (HOH)** model is a parametric survival model
constructed by enforcing positivity in the solution to the damped harmonic
oscillator differential equation. This produces a flexible hazard function capable
of capturing a wide range of shapes, including bathtub, unimodal, and
monotone, that arise naturally in applications such as cancer epidemiology and
reliability analysis, while retaining a mechanistic physical interpretation.

## Repository structure

```
HOH/
├── routinesHO.R     # Core functions required by all other scripts
├── appHO.R          # Shiny app to interactively visualise HOH shapes
├── Example.R        # Real data application from Christen & Rubio (2025)
├── Example.RData    # Saved output of Example.R
├── Example.Rmd      # R Markdown source for the real data application
├── Example.html     # Rendered tutorial with code and outputs
├── bound.R          # Upper bounds from the proof of Proposition 2
└── references.bib   # Bibliography
```

## Requirements

The following R packages are required. Install them before running any script:

```r
install.packages(c("shiny", "survival", "numDeriv"))
```

Source the core routines before running the example or the Shiny app:

```r
source("routinesHO.R")
```

## Shiny app

`appHO.R` launches an interactive Shiny app for exploring the shape of the HOH
hazard function across different parameter values:

```r
library(shiny)
source("appHO.R")
```

## Tutorial

- [On harmonic oscillator hazard functions](https://rpubs.com/FJRubio/HOH) — 
  illustrative RPubs tutorial accompanying the paper
- `Example.html` — rendered R Markdown tutorial available directly in this
  repository

## Related resources

- [ODESurv](https://github.com/FJRubio67/ODESurv) — broader framework for
  survival modelling via ordinary differential equations
- [SurvMODE](https://github.com/FJRubio67/SurvMODE) — hazard-based
  distributional regression via ODEs (R and Julia)
- [HazReg](https://github.com/FJRubio67/HazReg) — parametric hazard-based
  regression models (R package)

## Citation

If you use this code, please cite:

```bibtex
@article{christen:2025,
  author  = {Christen, J.A. and Rubio, F.J.},
  title   = {On harmonic oscillator hazard functions},
  journal = {Statistics \& Probability Letters},
  volume  = {217},
  pages   = {110304},
  year    = {2025},
  doi     = {10.1016/j.spl.2024.110304}
}
```
