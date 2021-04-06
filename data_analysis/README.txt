This readme.txt file was generated on 2020-11-30 by Hugo Salinas

GENERAL INFORMATION

1. Title of Dataset: Importance  of elongation and organogenesis on the rhizome length of Zostera marina in an individual-­based simulation model

2. Author Information

	A. Principal Investigator Contact Information
		Name: Hugo Salinas
		Email: hugosal@comunidad.unam.mx

	B. Co-investigator Contact Information
		Name: Elena Solana-Arellano
		Institution: Centro de Investigacion Cientifica y de Educacion Superior de Ensenada (CICESE)
		Address: Carretera Tijuana-Ensenada 3918, Ensenada, Mexico.
		Email: esolana@cicese.mx

3. Date of data collection: Fortnightly samples during 2000 and 2018

4. Geographic location of data collection: Punta Banda Estuary, Mexico. 31°40.9' N to 31°48.9' N, 116°37.9' W to 116°40.9' W.

OVERVIEW

This dataset contains the R scripts and files to perform the statistical analysis of the research paper "Importance of elongation and organogenesis on the rhizome length of Zostera marina in an individual­based simulation model."
The scripts include the methods used for the parametrization of an individual-based model, its validation, and sensitivity analysis of its parameters.
The model was written in the Python programming language and can be downloaded from https://github.com/hugosal/Zostera-marina-IBM
The model output files used for the validation and sensitivity analysis are included in this data set and can be recreated using the model.
The R scripts were implemented in R version 3.6.1, in Windows. 
The package that are used is the BayesFactor package (version 0.9.12-4.2), and the JAGS program (version 4.3.0) with the runjags package (version 2.0.4.6).

1. File List: 

environmental_time_series_plot_2000_Salinas_Solana.R
Fig4_Salinas_Solana
Fig5_Salinas_Solana
full_report_R.Rnw
internodes_and_branches_2000.csv
internode_length_rhizome_and_lateral.csv
internodes_and_branches_2000_environment.csv
internodes_and_branches_2018_environment.csv
latin_hipercube.csv
length_of_new_internodes_JAGS_sampling_Salinas_Solana.R
model_length_new_internodes_JAGS_output.csv
model_select_length_phyt.R
model_selection_n_phyt.R
new_internode_lenght_and_environmental.csv
new_internode_lenght_and_environmental_2018.csv
number_new_branches_JAGS_output.csv
number_new_branches_JAGS_sampling_Salinas_Solana.R
number_new_leaves_JAGS_output.csv
number_new_leaves_JAGS_sampling_Salinas_Solana.R
observed_rhizomes_2000.csv
observed_rhizomes_2018.csv
outputs_for_latin.csv
outs2000\
outs2018\
phytomer_age_mortality_plot_Salinas_Solana.R
phytomer_mortality_JAGS_output.csv
phytomer_mortality_JAGS_sampling_Salinas_Solana.R
README.txt
sensitivity_analysis_Salinas_Solana.R
shrink_factor_JAGS_output.csv
shrink_factor_JAGS_sampling_Salinas_Solana.R
validation_ANOVA_BF_2000_and_2018_Salinas_Solana.R
validation_n_phtomers_ANOVA_BF_2000_and_2018_Salinas_Solana.R

ANALYSIS STRUCTURE:
1) Subodel selection
The submodels of length of new internodes and number of new phytomers were selected after a model selection process.
The scripts model_select_length_phyt.R and model_selection_n_phyt.R can be used to re-create the submodel analysis used in the selection.
These scripts parametrize eight submodels and of each extract the DIC and the 0.95 credible interval of the distribution of its parameters.

2) Model parametrization 
The parametrization consisted of using the JAGS program thorugh the package. 
The model is decomposed into five submodels: number of new internodes in a fortnight, the length of new internodes, the mortality of internodes, the probability of branching, and the internode length shrink factor. 
The posterior distribution of the parameters of each submodel was approximated using the R scripts: 
	number_new_branches_JAGS_sampling_Salinas_Solana.R
	length_of_new_internodes_JAGS_sampling_Salinas_Solana.R
	phytomer_mortality_JAGS_output.csv
	number_new_leaves_JAGS_sampling_Salinas_Solana.R
	shrink_factor_JAGS_sampling_Salinas_Solana.R
Each script creates a CSV output file with a summary of the posterior distribution of the parameters of the submodel.
The parametrization was performed using data of observed plants during the year 2000 in the study site.
The values of these parameters were used for the construction of the individual-based model in Python.

3) Model validation
The model's validation consisted of comparing several simulations of the model with the observed plants during the years 2000 and 2018. 
The R script validation_ANOVA_BF_2000_and_2018_Salinas_Solana.R was used for this comparison using the BayesFactor package.
The R script uses the output files in the directories outs2000 and outs2018 to compare the simulated with the observed plants (observed_rhizomes_2000 and observed_rhizomes_2018.csv).
Figure 3 shows one comparison of observed and simulated plants through time and was created with the script Fig3_Salinas_Solana.R

4) Sensitivity analysis
The sensivity analysis consisted of running simulations using different values of the model parameters and measuring its effect on the model output. 
The Latin Hypercube method was used to generate a set of parameter values (sensitivity_sampling.py, in the model source code)
The R script sensitivity_analysis_Salinas_Solana.R was used to compute the correlation of the parameter values (latin_hipercube.csv) and the output measure (outputs_for_latin.csv).

The file full_report.Rnw contains all the analysis and creates a PDF and TEX report of the results.


DATA-SPECIFIC INFORMATION FOR: internodes_and_branches_2000.csv

1. Number of variables: 8

2. Number of cases/rows: 443

date_stamp: Date stamp representing the date of the beginning of fortnight, date stamp from 1970/01/01
number_of_sampled_rhizome: Identifier of each rhizome to which an internode belonged
number_internodes_in_rhizome: Number of internodes that comprised a rhizome
new_leaves_number: Number of leaves that were added to a rhizome during the fortnight, estimated using leaf marking.
number_final_branches: Number of branches of an individual at the end of the fortnight, including the rhizome and lateral branches
number_initial_branches: Number of branches of an individual at the beginning of the fortnight, including the rhizome and lateral branches
new_branches_per_branch: Number of new branches added to the individual per initial branch; this is the number of new branches divided by the number of initial branches
Date: Date of the beginning of the fortnight

DATA-SPECIFIC INFORMATION FOR: internode_length_rhizome_and_lateral.csv

1. Number of variables: 3

2. Number of cases/rows: 65

Individual: the identifier number of each individual
rhizome_0_else_lateral_1: this variable indicates if an internode belongs to a rhizome (0), or to a lateral branch (1)
internode_length: Length of the internode (mm)

DATA-SPECIFIC INFORMATION FOR: internodes_and_branches_2000_environment 

1. Number of variables: 12

2. Number of cases/rows: 268

date_stamp: Date stamp representing the date of the beginning of fortnight, date stamp from 1970/01/01
number_of_sampled_rhizome: Identifier of each rhizome to which an internode belonged
number_internodes_in_rhizome: Number of internodes that comprised a rhizome
new_leaves_number: Number of leaves that were added to a rhizome during the fortnight, estimated using leaf marking.
number_final_branches: Number of branches of an individual at the end of the fortnight, including the rhizome and lateral branches
number_initial_branches: Number of branches of an individual at the beginning of the fortnight, including the rhizome and lateral branches
new_branches_per_branch: Number of new branches added to the individual per initial branch; this is the number of new branches divided by the number of initial branches
SST: Sea surface temperature (°C)
SSTA Sea Surface Temperature Anomaly (°C)
Irradiance: Irradiance (kW fortnight^-1 m^-2)
exposition_hours: Hours of air at a 0.5 m depth (h)
Date: Date of the beginning of the fortnight

DATA-SPECIFIC INFORMATION FOR: internodes_and_branches_2018_environment

1. Number of variables: 12

2. Number of cases/rows: 103

date_stamp: Date stamp representing the date of the beginning of fortnight, date stamp from 1970/01/01
number_of_sampled_rhizome: Identifier of each rhizome to which an internode belonged
number_internodes_in_rhizome: Number of internodes that comprised a rhizome
new_leaves_number: Number of leaves that were added to a rhizome during the fortnight, estimated using leaf marking.
number_final_branches: Number of branches of an individual at the end of the fortnight, including the rhizome and lateral branches
number_initial_branches: Number of branches of an individual at the beginning of the fortnight, including the rhizome and lateral branches
new_branches_per_branch: Number of new branches added to the individual per initial branch; this is the number of new branches divided by the number of initial branches
SST: Sea surface temperature (°C)
SSTA Sea Surface Temperature Anomaly (°C)
Irradiance: Irradiance (kW fortnight^-1 m^-2)
exposition_hours: Hours of air at a 0.5 m depth (h)
Date: Date of the beginning of the fortnight

DATA-SPECIFIC INFORMATION FOR: latin_hipercube.csv

1. Number of variables: 99

2. Number of cases/rows: 16

Each of the colums correspond to a set of values of a parameter of the model.
the rows correspond to the parameters of the model:
1: n_lambda
2: l_sigma_alpha
3: l_sigma_beta
4: l_mu_T
5: l_mu_T^2
6: l_mu_A^2
7: l_mu_I
8: l_mu_H
9: a_shrink
10: b_shrink
11: a_branchloss
12: b_branchloss
13: a_branching
14: b_branching
15: m_alpha
16: m_beta



DATA-SPECIFIC INFORMATION FOR: model_length_new_internodes_JAGS_output.csv, number_new_branches_JAGS_output.csv, number_new_leaves_JAGS_output.csv, phytomer_mortality_JAGS_output.csv, shrink_factor_JAGS_output.csv

1. Number of variables: 12

2. Number of cases/rows: one for each parameter in the submodel

The variables are the summary values from the runjags package:
Lower95: The lower confidence limit for the highest posterior density
Median: The median value
Upper95: The upper confidence limit for the highest posterior density
Mean: The mean value
SD: The sample standard deviation
Mode: The mode of the variable
MCerr: The Monte Carlo standard error associated with this variable
MC.ofSD: The Monte Carlo standard error expressed as a percentage of the standard deviation of the variable
SSeff: The effective sample size
AC.100: The autocorrelation of the sample
psrf: The potential scale reduction factor


DATA-SPECIFIC INFORMATION FOR: new_internode_lenght_and_environmental.csv

1. Number of variables: 9

2. Number of cases/rows: 1157

date_stamp: Date stamp representing the date of the beginning of fortnight, date stamp from 1970/01/01
rhizome_number: Identifier number of each rhizome to which an internode belonged
Numbrer_of_internodes_in_rhizome: Number of internodes that comprised the sampled rhizome
length_of_internode: Length of the internode (mm)
SST: Sea surface temperature (°C)
SSTA Sea Surface Temperature Anomaly (°C)
Irradiance: Irradiance (kW fortnight^-1 m^-2)
exposition_hours: Hours of air at a 0.5 m depth (h)
Date: Date of the beginning of the fortnight



DATA-SPECIFIC INFORMATION FOR: observed_rhizomes_2000.csv, observed_rhizomes_2018.csv

1. Number of variables: 3

2. Number of cases/rows: 5158, 1862

Date: Date of the beginning of the fortnight
Rhizome:  Identifier of each rhizome to which an internode belonged
Internode_length: Length of the internode (mm)



DATA-SPECIFIC INFORMATION FOR: outputs_for_latin.csv

1. Number of variables: 1

2. Number of cases/rows: 100

Each row in this file correspond to the final mean rhizome length (cm) of a simulation



DATA-SPECIFIC INFORMATION FOR: files in directories outs2000 and outs2018

Each of the files in these directories correspond to the output of a simulation

Date: Date of the beginning of the fortnight
Rhizome: Identifier of each rhizome to which an internode belonged
Internode number: Position number of the internode in the rhizome
Internode length: Length of the internode (mm)
