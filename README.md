# Dependent Dirichlet processes via thinning


[![arXiv:2506.18223](https://img.shields.io/badge/arXiv-2506.18223-b31b1b.svg)](https://arxiv.org/abs/2506.18223)
![Last Commit](https://img.shields.io/github/last-commit/laura-dangelo/thinned_DDP)

This repository contains the code to replicate the results in the paper "_Dependent Dirichlet processes via thinning_" by D'Angelo, Nipoti, and Ongaro (available on [arXiv](https://arxiv.org/abs/2506.18223)). It also contains the R package `thinnedDDP`, which implements the Gibbs sampler for posterior inference (see [Installation](#installation)).


## Structure of the Repository

Main structure:
<ul>
  <li> 01_Simulation_study : folder that contains the scripts and outputs of the simulation study; </li>
  <li> 02_Application : folder that contains the script and outputs of the real data analysis;  </li>
  <li> 03_Sensitivity_study : folder that contains the scripts and outputs of the sensitivity study;  </li>
  <li> thinnedDDP : R package implementing the Gibbs sampler. </li>
  <li> thinnedDDP_project : R project to be used as main directory for running the code.</li>
</ul>


Large files (typically, outputs of the Gibbs sampler algorithm and extracted quantities) are excluded from the repository and are available in the Google Drive [folder](https://drive.google.com/drive/folders/1VVzDpo1fLS5QG8ByCHO6oXbZSlVjNahu?usp=sharing). The name of each subfolder corresponds to the path of the files in the repository. Note that these files are not necessary to produce the plots in the article, which can be obtained starting from the RDS files available in the folders named `output_RDS`.  


