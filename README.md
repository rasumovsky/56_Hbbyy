# A search for di-Higgs production with ATLAS

### Introduction
This package implements a search for non-resonant and resonant di-Higgs
production in ATLAS data at the LHC using analytic statistical models. The two
Higgs bosons are tagged using the diphoton and bb final states, to provide 
a combination of good resolution and clean triggers and large cross-section. 

The code has been structured so that all of the general analysis settings are 
stored in the configuration file in the data/ directory.

If you have questions about the code, please contact the author, who will be happy to assist you in getting the code running.

##### Input files

The analysis requires text files with the mass points from data (or a suitable MC sample) in each category. 

##### Config files

The configuration files for the analysis are included in the data directory. Currently, the di-Higgs analyses use the files below:

- settingsDH_nonresonant.cfg (for the non-resonant analysis),
- settingsDH_resonant.cfg (for the fit-based resonant analysis),
- settingsDH_countresonant.cfg (for the counting resonant analysis).

The config files allow the user to specify everything from the luminosity, to the blinding status, to the systematic uncertainty sizes, without need to recompile or edit C++ code.

##### Running the code

Every program in this package can be run via the "DHMaster" main method. 

 - Workspace (make workspaces for the analysis)
 - TossPseudoExp (create a pseudoexperiment ensemble for mu=0 and mu=1)
 - PlotPseudoExp (plot the results of the step above)
 - CLScanSubmitToys (submit toys for the CL scan, if toy result desired).
 - CLScanAnalysis (get the 95% CL value from a CL scan (toys or asymptotics!)
