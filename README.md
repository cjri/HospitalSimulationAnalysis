# HospitalSimulationAnalysis
Code for analysing simulated data describing nosocomial infection in hospitals

This repository contains code to analyse simulated data describing nosocomial transmission in hospitals, modelling an intervention in which UDCA was given to patients on a ward following the detection of a case of infection.  The code is written in C++ and can be compiled with a make command.

In compiling this code note that:
  The code makes use of the GSL library.  You will need to have this installed on your machine.
  The Makefile may need editing to point it towards an appropriate C compiler, and the GSL library.

Example scripts to run the code are located in Code/Calculations/Output.  These scripts assume that:
  You have the GNU parallel package installed on your machine.
  Your machine has access to at least 12 computational cores.  This assumption can be altered by editing the script, replacing --jobs 12 with --jobs X where X is the number of cores you want to use to run the calculations.
   You have available files describing simulated data describing outbreaks in hospitals.  These are available from Stephanie Evans by request.

In brief, the script iterates over running the code on simulations 1 to 180, denoted by the {} marker in each script.

Scripts run the code for an intervention length of between 6 and 14 days e.g. run_calculations12.sh runs for an intervention length of 12 days from the detection of a case.

The directory Intervention10 directory contains the outputs of running the code for a 10-day intervention.  The 5 .gz files contain the outputs from each of 5 replica runs for the 180 simulated hospital datasets.  Key files used in our analysis were:

Cluster_sizes.dat: Contains the number of individuals in each SARS-CoV-2 outbreak, as a distribution, with the first column containing the number of individuals, and the second column containing the number of instances.  An outbreak of size 1 is where an individual with SARS-CoV-2 did not infect anybody else.

RunX.out: The output of the model.  The critical lines here are the final two, which detail
#Nosocomial patient infections (original) 	#Nosocomial patient infections (intervention) #Percentage remaining
#Nosocomial HCW infections (original) 	#Nosocomial HCW infections (intervention) #Percentage remaining

The files in the Intervention10 directory contains the outputs of five replicate runs of the code for simulated data given a 10-day intervention.


The folder /Data contains information of where in each hospital the intervention took place.  File names are in the format:

Ward_detection_windows_patient_onlyX_Y.dat

Here X is the number of the simulation, from 1 to 180.  Y is the number of days modelled in the intervention.  The format of each file is:

Ward number  Intervention start (day)  Intervention end (day)


The file FindInterventionWindows.nb is a Mathematica worksheet that was used to generate the ward detection windows provided in the Data directory.

Full simulation data will be made available upon the publication of a manuscript accompanying this code.
