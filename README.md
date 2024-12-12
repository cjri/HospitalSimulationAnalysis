# HospitalSimulationAnalysis
Code for analysing simulated data describing nosocomial infection in hospitals

This repository contains code to analyse simulated data describing nosocomial transmission in hospitals, modelling an intervention in which UDCA was given to patients on a ward following the detection of a case of infection.  The code is written in C++ and can be compiled with a make command.

In compiling this code note that:
  The code makes use of the GSL library.  You will need to have this installed on your machine.
  The Makefile may need editing to point it towards an appropriate C compiler, and the GSL library.

Example scripts to run the code are located in Code/Calculations/Output.  These scripts assume that:
  You have the GNU parallel package installed on your machine.
  Your machine has access to at least 12 computational cores.  This assumption can be altered by editing the script.

In brief, the script iterates over running the code on sumulations 1 to 180, denoted by the {} marker in each script.

Scripts run the code for an intervention length of between 6 and 14 days.

The files in the Intervention10 directory contains the outputs of five replicate runs of the code for simulated data given a 10-day intervention.


The folder /Data contains information of where in each hospital the intervention took place.  The format is

Ward number  Intervention start (day)  Intervention end (day)


The file FindInterventionWindows.nb is a Mathematica worksheet that was used to generate the ward detection windows provided in the Data directory.

Full simulation data will be made available upon the publication of a manuscript accompanying this code.
