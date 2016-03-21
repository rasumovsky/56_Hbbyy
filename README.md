# A search for di-Higgs production with ATLAS

### Introduction
This package implements a search for non-resonant and resonant di-Higgs
production in ATLAS data at the LHC using analytic statistical models. The two
Higgs bosons are tagged using the diphoton and bb final states, to provide 
a combination of good resolution and clean triggers and large cross-section. 

The code has been structured so that all of the general analysis settings are 
stored in the configuration file in the data/ directory.

If you have questions about the code, please contact the author, who will be 
happy to assist you in getting the code running.

##### Input files

The analysis requires text files with the mass points from data (or a suitable 
MC sample) in each category. For ATLAS users with access to afs, the files 
corresponding to 3.2 fb-1 of 13 TeV data collected in 2015 can be found in the 
directory: `/afs/cern.ch/user/a/ahard/public/ForDiHiggs/DHMassPoints/`.

##### Config files

The configuration files for the analysis are included in the data directory. 
Currently, the di-Higgs analyses use the files below:
 - settingsDH_nonresonant_no1b.cfg (non-resonant analysis),
 - settingsDH_countresonant_no1b.cfg (counting resonant analysis).


The config files allow the user to specify everything from the luminosity, to 
the blinding status, to the systematic uncertainty sizes, without any need to 
recompile or edit C++ code.

##### Compiling and Running the Package

There are several executables that must be compiled prior to running. This 
package currently uses a custom makefile (makefile). Before compiling, set up
the ROOT and GCC versions:
```
source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh 
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh
```

Alternatively, one can just run the setup script:
```
source scripts/package_setup.sh
```

To compile the executables:
```
make bin/DHMaster
make bin/DHPseudoExp
make bin/DHCLScan
make bin/DHNuisanceParameters
make bin/DHPlotCLvsMX
make bin/DHNLLScan
```

Every program in this package can be run via the "DHMaster" main method. The 
user only needs to specify the step of the analysis to run and the config file
with all of the analysis settings of interest. To produce the workspace for the
non-resonant analysis, simply execute the following command:
```
./bin/DHMaster Workspace data/settingsDH_nonresonant_no1b.cfg
```

The analysis steps are listed and explained below.
 - Systematics: Load, store, combine, parameterize systematic uncertainties.
 - Workspace: Make workspaces for the analysis (resonant or non-resonant).
 - TossPseudoExp: Create a pseudoexperiment ensemble for mu=0 and mu=1.
 - PlotPseudoExp: Plot the results of fits to a pseudo experiment ensemble.
 - TestStat: Get asymptotic p0 and CL statistical results.
 - CLScanSubmitToys: Submit toys for the CL scan, if toy result desired.
 - CLScanAnalysis: Get the 95% CL value from a CL scan (toys or asymptotics!).
 - PlotCLVsMX: Plot the resonant analysis limits vs. resonant mass.
 - ScanNLL: Scan the NLL for the cross-section in order to compare errors.
 - RankNP: Rank and plot the nuisance parameters in the analysis.