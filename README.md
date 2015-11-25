# A search for di-Higgs production with ATLAS

### Introduction
This package implements a search for non-resonant and resonant di-Higgs
production in ATLAS data at the LHC using analytic statistical models. The two
Higgs bosons are tagged using the diphoton and bb final states, to provide 
a combination of good resolution and clean triggers and large cross-section. 

The code has been structured so that all of the general analysis settings are 
stored in the configuration file in the data/ directory.

NOTE: This package has recently undergone significant changes. The primary difference is that the statistical model for fitting is now described entirely within the config file. The nuisance parameters, global observables, and constraint terms for systematic uncertainties are handled automatically based on specifications in the config file.