# CompSurv - an R package for competing risk and misc survival modelling

## Aim

_CompSurv_ is an R package that aims to offer tools for R users to handle competing risk analyses in survival modelling and more comprehensively handle also survival modelling in general. As such, the package is closely related to e.g. Fine-Gray models (competing risks), Cox models, and Kaplan-Meiers.
Some key functionality includes model diagnostics (such as inspecting proportional hazards assumptions or colliearity of variables) and visualizations (such as 

## Dependencies

CompSurv's intentionally kept lean, with only very few dependencies beyond ```base``` R, such as ```survival```. Plots aim to be compatible with the base R graphics.

## Author(s)

Teemu Daniel Laajala (contact teelaa@utu.fi or daniel.laajala@helsinki.fi)

## Acknowledgements / citations / etc

CompSurv's been used (at least) in the following publications, for visualizations, diagnostics, etc:

```
@article{brannback2025competing,
  title={Competing risk of death in patients with low, intermediate and high risk of recurrence after radical surgery for clear cell renal cell carcinoma},
  author={Br{\"a}nnb{\"a}ck, Anna and Mustonen, Ivan and Laajala, Teemu D and Vainio, Paula and Lindskog, Magnus and Kjellman, Anders and Lundgren, Per-Olof and Jaakkola, Panu M and Mattila, Kalle E},
  journal={BJUI compass},
  volume={6},
  number={7},
  pages={e70047},
  year={2025},
  publisher={Wiley Online Library}
}
``` 

To cite the package itself, use standard R package citation style (such as URL of this git repo and version of the package you used).