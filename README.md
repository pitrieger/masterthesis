# Measurement Invariance in Confirmatory Factor Analysis: Methods for Detecting Non-invariant Items
The most recent draft version of the thesis is available [here](https://www.overleaf.com/read/ggznrkmtkwxp). Feel free to [contact me](mailto:prieger@ethz.ch?subject=[Transparent%20Master%20Thesis]%20Inquiry) with questions, suggestions, etc. An interactive version of the simulation study in my thesis is available as a shiny app [here](https://prieger.shinyapps.io/miapp/).

## Abstract
Violations of measurement invariance (MI) of a given confirmatory factor analysis (CFA) model can arise as a result of non-invariant items and pose a significant threat to the validity of latent variable comparisons across subgroups of a study population. While methods for detecting such items under partial MI exist, there hasn't been a systematic study to compare their performance. This thesis makes three contributions. First, by means of a simulation study, the performance of six detection methods is assessed. Second, two versions of a novel detection approach are introduced and included in the simulation study. The advantage of the novel approach is that it is arguably much easier to interpret than existing methods. Instead of relying on likelihood inference, it builds on residuals and only requires a basic understanding of linear regression, thus being much more accessible to a broad audience of applied researchers. Finally, the detection methods are applied to different CFA models for measuring populist attitudes using survey data, demonstrating that they can easily be generalized to fairly complex measurement models. 
In terms of performance, the findings indicate that one of the existing methods and one of the novel methods can reliably detect non-invariant items. Other existing approaches perform the task they were developed for so poorly that they cannot be recommended under any circumstance. For the exemplary application, the results corroborate findings of significant issues with respect to cross-cultural validity at the model level, but also provide a starting point for model improvement to be taken up by further research.

## Navigating this Repo
This repo contains all the resources to replicate and analyze the simulation study and application with the exception of data from Castanho et al.'s (2020) paper, which can be downloaded from the authors' own repo [here](https://github.com/bcastanho/PRQ2019). All relevant scripts are included in the folder [/Rscripts](https://github.com/pitrieger/masterthesis/tree/main/Rscripts). The implementations of the detection methods for single- and multi-factor models are contained in the sub-folders [/simulation](https://github.com/pitrieger/masterthesis/tree/main/Rscripts/simulation) and [/application](https://github.com/pitrieger/masterthesis/tree/main/Rscripts/application), respectively, which corresponds with where they are used in the paper. Within these folders, the main scripts for running and analyzing the simulation study are [Simulation_final.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/simulation/Simulation_final.R) and [Analysis_final.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/simulation/Analysis_final.R), respectively. The main script for running the application is [application_main.R](https://github.com/pitrieger/masterthesis/blob/main/Rscripts/application/application_main.R). Simulation data is available in [/publicdata](https://github.com/pitrieger/masterthesis/blob/main/publicdata).

## Current work
Doing some final polishing and proofreading.

## Open Questions


## Log
- 16/08/2021 - Finished initial proposal
- 15/09/2021 - Registered thesis with study administration
- 05/10/2021 - Officially started writing
- 11/2021 - Finished draft for introduction of (confirmatory) factor analysis
- 11/2021 - Finished draft for introduction of measurement invariance
- 11/2021 - Implemented detection methods (existing and original contribution) for single-factor models
- 11/2021 - Implemented preliminary simulation study
- 12/2021 - Implemented detection methods for more general CFA models
- 12/2021 - Finished draft for introduction of detection methods
- 12/2021 - Finished draft for simulation study
- 01/2022 - Finished draft introduction
- 01/2022 - Finished draft for application
- 01/2022 - Programmed Shiny app for interactive version of simulation study
- 01/2022 - Formalized idea of novel method (to be added to draft)
- 02/2022 - Ran final version of simulation study
- 02/2022 - Added formalized novel method introduction


<!--- ![70%](https://progress-bar.dev/70) --->


