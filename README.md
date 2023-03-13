
# Readme <a href='https://osf.io/zcvbs/'><img src='worcs_icon.png' align="right" height="139" /></a>

<!-- Please add a brief introduction to explain what the project is about    -->

## Where do I start?

You can load this project in RStudio by opening the file called 'ERV_simulation.Rproj'.
Afterwards, follow the "Reproducibility" section to reproduce the simulation and reanalysis results.

## Project structure

<!--  You can add rows to this table, using "|" to separate columns.         -->
Folder| Description                | Usage         
------------------------- | -------------------------- | --------------
manuscript | Folder that holds manuscript pdfs (empty at the moment)      | Human editable
RSD                 | R scripts that calculate [relative variability (Mestdagh et al., 2018)](https://ppw.kuleuven.be/okp/software/relative_variability/) | Read only     
codebook| Codebooks for [dataset 1 - 3 (Blanke et al., 2020)](https://osf.io/mxjfh/)      | Read only



File                      | Description                | Usage         
------------------------- | -------------------------- | --------------
README.md                 | Description of project     | Human editable
ERV_simulation.Rproj      | Project file               | Loads project 
LICENSE                   | User permissions           | Read only     
prepare_data.R            | Script to process raw data | Run to reproduce  results
sim1_VAR(1).R                 | Simulation 1: data generation with VAR(1) | Run to reproduce  results
sim1_performance.R                 | Simulation 1: evaluate performance | Run to reproduce  results
sim2_lorenz.R                 | Simulation 2: data generation with Lorenz system | Run to reproduce  results
reanalysis_main.R                 | Reproducible R environment | Run to reproduce  results
.worcs                    | WORCS metadata YAML        | Read only     
renv.lock                 | Reproducible R environment | Read only     
dfraw1.csv                | [Dataset 1 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     
dfraw2.csv                | [Dataset 2 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     
dfraw3.csv                | [Dataset 3 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     
metric_functions.R                 | Reproducible R environment | Human editable
reanalysis_MLM.R                 | Reproducible R environment | Human editable
reanalysis_calculateERV.R                 | Reproducible R environment | Human editable

<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

# Reproducibility

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to
ensure transparency and reproducibility. The workflow is designed to meet the
principles of Open Science throughout a research project. 

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles,
read the preprint at https://osf.io/zcvbs/


## WORCS: Advice for readers

Please refer to the vignette on [reproducing a WORCS project]() for step by step advice.
<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->
