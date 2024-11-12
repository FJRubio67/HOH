# On harmonic oscillator hazard functions

This repository contains information about the project `HOH`, which concerns the construction of a parametric hazard model obtained by enforcing positivity in the damped harmonic oscillator. In particular, this repository contains real data applications and R code for the survival models in:

> Christen, J. A., and Rubio, F. J. (2024+). On harmonic oscillator hazard functions. *Statistics & Probability Letters*, in press.https://doi.org/10.1016/j.spl.2024.110304. Preprint: https://arxiv.org/abs/2408.15964

- The code `appHO.R` produces a shiny app to visualise the shapes of the harmonic oscillator hazard function for different values of the parameters.
- The code `routinesHO.R` contains the functions required in other codes presented in this repository.
- The code `Example.R` implements the real data example in Christen and Rubio (2024) (the output of this file is stored in `Example.RData`). The same example is also available in R Markdown in the files `Example.Rmd` and `Example.html`.
- The code `bound.R` implements the upper bounds presented in the proof of Proposition 2.

See also:
- [On harmonic oscillator hazard functions](https://rpubs.com/FJRubio/HOH)
- [ODESurv](https://github.com/FJRubio67/ODESurv)
