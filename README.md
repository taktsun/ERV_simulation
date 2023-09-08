



# Readme <a href='https://osf.io/zcvbs/'><img src='worcs_icon.png' align="right" height="139" /></a>

This is the github of my manuscript ["A Theory-Informed Emotion Regulation Variability Index: Bray-Curtis Dissimilarity"](https://psyarxiv.com/twk9m/download?format=pdf). You can reproduce the analysis results in the manuscript following this readme. 

## Where do I start?

You can load this project in RStudio by opening the file called 'ERV_simulation.Rproj'.
Afterwards, follow the "Reproducibility" section to reproduce the simulation and reanalysis results.

## Project structure

Tables of R scripts, other files, and folders.

### R scripts

File                      | Description                | Usage         
------------------------- | -------------------------- | --------------
prepare_data.R            | Reanalysis: load raw data from OSF | Optional because data are already in .csv
sim1_VAR(1).R                 | Simulation 1: data generation with VAR(1) | Run to reproduce results
sim2_lorenz.R                 | Simulation 2: data generation with Lorenz system | Run to reproduce results
reanalysis_main.R                 | Reanalysis: descriptive statistics and main analyses | Run to reproduce results 
RMSE_bootstrap.R                 | Produce bootstrapped RMSEs following reanalysis | Run to reproduce results
sim1_VAR(1)_measurement.R                | Extra analyses for supplemental materials| Run to reproduce results
sim2_lorenz_analysis_measurement.R                | Extra analyses for supplemental materials| Run to reproduce results
func_indices.R                 | Dissimilarity indices for simulations & reanalyses | Required; read only
func_rean_desStat.R                 | Reanalysis 1: descriptive statistics | Required; read only
func_rean_calculateERV.R                 | Reanalysis 2: calculate dissimilarity indices | Required; read only
func_rean_MLM.R                 | Reanalysis 3: multilevel models | Required; read only
packages_needed.R                 | States package dependency for renv | Required; read only



### Other Files

File                      | Description                | Usage         
------------------------- | -------------------------- | --------------
ERV_simulation.Rproj      | Project file               | Loads project 
sim1_input.RData      | Simulation 1 parameter input with seed               | Optional; load with ReadRDS in step 1 for replication
README.md                 | Description of project     | Read only
LICENSE                   | User permissions           | Read only     
.worcs                    | WORCS metadata YAML        | Read only     
renv.lock                 | Reproducible R environment | Read only     
dfraw1.csv                | [Dataset 1 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     
dfraw2.csv                | [Dataset 2 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     
dfraw3.csv                | [Dataset 3 (Blanke et al., 2020)](https://osf.io/mxjfh/) for reanalysis| Read only     


### Folders
Folder| Description                | Usage         
------------------------- | -------------------------- | --------------
manuscript | Folder that holds manuscript pdfs (empty at the moment)      | Read only
RSD                 | R scripts that calculate [relative variability (Mestdagh et al., 2018)](https://ppw.kuleuven.be/okp/software/relative_variability/) | Required; Read only     
codebook| Codebooks for [dataset 1 - 3 (Blanke et al., 2020)](https://osf.io/mxjfh/)      | Read only



<!--  You can consider adding the following to this file:                    -->
<!--  * A citation reference for your project                                -->
<!--  * Contact information for questions/comments                           -->
<!--  * How people can offer to contribute to the project                    -->
<!--  * A contributor code of conduct, https://www.contributor-covenant.org/ -->

# Reproducibility
Reproduce the results by these 5 steps.

 1. Install RStudio and R
 2. Install WORCS dependencies
		
		install.packages("worcs", dependencies = TRUE)
		tinytex::install_tinytex()
		renv::consent(provided = TRUE)
		
 3. [Clone](https://resources.github.com/github-and-rstudio/#:~:text=Clone%20the%20repository%20with%20RStudio&text=On%20GitHub%2C%20navigate%20to%20the,RStudio%20on%20your%20local%20environment.) this repo (https://github.com/taktsun/ERV_simulation) to your RStudio
 4. Restore the package dependencies
	

	    renv::restore()

 5. Run 3 R scripts to reproduce the results. Start new R session (Ctrl+Shift+F10 in Windows) before you run each simulation or reanalysis.
 
	- sim1_VAR(1).R
	- sim2_lorenz.R
	- reanalysis_main.R, then RMSE_bootstrap.R

Step 1 to 4 are detailed in the vignette on [reproducing a WORCS project](https://cjvanlissa.github.io/worcs/articles/reproduce.html).

<!-- If your project deviates from the steps outlined in the vignette on     -->
<!-- reproducing a WORCS project, please provide your own advice for         -->
<!-- readers here.                                                           -->

### Adherence to WORCS

This project uses the Workflow for Open Reproducible Code in Science (WORCS) to ensure transparency and reproducibility. The workflow is designed to meet the principles of Open Science throughout a research project. We used WORCS for the simulation studies and data analysis, but we did not use WORCS for manuscript preparation.

### More about WORCS

To learn how WORCS helps researchers meet the TOP-guidelines and FAIR principles, read the preprint at https://osf.io/zcvbs/

