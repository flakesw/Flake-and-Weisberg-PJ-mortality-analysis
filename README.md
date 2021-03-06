# Flake-and-Weisberg-PJ-mortality-analysis

<a href="https://zenodo.org/badge/latestdoi/111577553"><img src="https://zenodo.org/badge/111577553.svg" alt="DOI"></a>

Analysis by Sam Flake, swflake@ncsu.edu

Abstract: 

Severe drought has resulted in widespread tree die-off events in forests and woodlands globally and is forecast to become more frequent in coming decades. Tree mortality is a complex process influenced by climate, soils, characteristics of individual trees, interactions between trees, and the dynamics of pests and pathogens. The role of stand structure and stand density in mediating the resistance of trees to drought remains poorly understood, especially in semiarid woodlands, which are expected to be highly susceptible to future severe drought. We sampled permanent plots in central Nevada woodlands dominated by single-leaf pinyon pine and Utah juniper before and after a severe multi-year drought (2013-2015) to investigate the importance of climate, tree attributes, and local-neighborhood stand structure on tree mortality and canopy dieback at the level of individual trees and 0.1-ha plots.
	
  We observed widespread tree mortality of pinyon at approximately eight times the reported background mortality rate, and substantial canopy dieback in both pinyon and juniper. Both species were more prone to mortality and dieback in hotter, drier sites. Canopy dieback was associated with both long-term average climate and the severity of recent drought, with elevated mortality on sites with higher water deficits, average summer temperatures, and vapor pressure deficits. Soils also played a role in tree dieback, with greater mortality on deeper soils. 
  
  While mortality was driven largely by climate at coarse scales, fine-scale stand structure interacted with climate to mediate mortality and dieback. Neighborhood statistics showed that trees were susceptible to competitive influence, and pinyon trees were especially sensitive to neighborhood density on drier sites. Mortality and dieback were associated with diverse, co-occurring insect and parasitic plant mortality agents. Canopy dieback prior to the drought was strongly associated with tree mortality during the drought, implying that current widespread defoliation caused by these agents may foreshadow future elevated woodland decline. Fine-scale influences such as stand structure and soil characteristics play a key role in the long-term dynamics of semiarid woodlands, and these factors should be considered in predictive models of forest and woodland susceptibility to drought.

Description of files:

This code takes raw data, found in the folder "Raw Data," processes it into clean analyzable data, and generates all the figures used in the manuscript Flake, S.W., P.J. Weisberg. Accepted. Fine-scale stand structure mediates drought-induced tree mortality in pinyon-juniper woodlands. Ecological Applications.

The first script that should be run is pj_mortality_data_prep_112117.R, which processes the raw data and populates the Clean Data folder. This script merges climate data, plot-level data, and tree-level data. Some plot-level variables (for example, basal area) are aggregated from tree-level data. 

After generating clean data, plot_level_analysis_for_paper_061718.R and ind_tree_analysis_for_paper_061718.R can be run. These scripts run all of the analysis for the paper, including variable selection, model fitting, bootstrapping, and validation. It takes processed data from the Clean Data folder and outputs several large files to the model output folder, including saved model objects as .Rds files, which are used for generating figures. This script will take AT LEAST 1 Day to run on a run-of-the-mill laptop, so if you are interested in the analysis, I'd recommend skipping the model dredging and bootstrapping and looking at the exported objects instead. 

For each response variable (juniper dieback, pinyon dieback, pinyon mortality), three objects are saved as .Rds files: the model selection object from MuMIn, the output from the bootstrapping process (from package boot), and the final model. The summarized results from model selection (tables in the appendix) and bootstrapping (with confidence intervals) are saved as a .csv file. 

The final script generates the effects plots and other figures (Figures 1-6 and Figures S1 and S2). It requires the .rds files for the fitted models as well as the bootstrapped parameter estimates. It also needs the raw data and the clean data generated by the data prep script. Figures are output to the ./plots folder. 

The code here is not very well documented and a little spaghetti-stringy, so please contact me if you have any questions. 
