using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>


#include "process_simulation.h"
#include "io.h"
#include "utilities.h"

int main(int argc, const char **argv) {
    
    int seed=(int) time(NULL);
    gsl_rng_env_setup();
    gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rgen, seed);
    
    run_params p;
    GetOptions(p,argc,argv);
    
    //Set up conditional distribution of level of exposure and baseline ACE2 expression.  We want the quantile of each one and the distributions.
    
    //First set up basic distributions
    vector<double> exposures;
    vector<double> ace_basic;
    vector<double> ace12;
    vector<double> ace3;
    SetupDistributions (exposures,ace_basic,ace12,ace3);
    
    //Find conditional distribution
    vector< vector<double> > condo;
    MakeConditionalDistribution (exposures,ace_basic,condo);

    //Find partial probabilities
    vector<double> coord1;
    vector< vector<double> > coord2;
    FindPartialProbabilities (exposures,ace_basic,condo,coord1,coord2);
    
    int choice1;
    int choice2;
    SampleProbabilities (choice1,choice2,exposures,coord1,coord2,rgen);

    //Calculate reduction matrices
    vector< vector<double> > reductions12;
    vector< vector<double> > reductions3;
    GetReductions (exposures,ace_basic,ace12,ace3,reductions12,reductions3);

  
    //Import patient data
    cout << "Import patient data\n";
    vector<precord> pat_records;
    ReadPatientData(p,pat_records);
    cout << "Import HCW data\n";
    vector<precord> hcw_records;
    ReadHCWData(p,hcw_records);
    cout << "Import intervention data\n";
    vector< vector<int> > interventions;
    ReadInterventionData (p,interventions);
    
    ofstream inv_file;
    inv_file.open("Interventions.dat");
    cout << "Interventions\n";
    for (int i=0;i<interventions.size();i++) {
        inv_file << interventions[i][0] << " " << interventions[i][1] << " " << interventions[i][2] << "\n";
    }
    inv_file.close();
    
    //Sort out dates.  Convert to days since 1/1/2020
    GenerateYMD(pat_records);
    GenerateYMD(hcw_records);
    
    //Find unique individuals
    vector<int> patients;
    vector<int> hcws;
    FindUniqueIndividuals(pat_records,patients);
    FindUniqueIndividuals(hcw_records,hcws);
    
    //Make set of records corresponding to each individual
    vector< vector<precord> > indiv_pat_records;
    vector< vector<precord> > indiv_hcw_records;
    SplitRecordsByIndividual (patients,pat_records,indiv_pat_records);
    SplitRecordsByIndividual (hcws,hcw_records,indiv_hcw_records);
    
    //cout << "Size " << indiv_pat_records.size() << " " << indiv_hcw_records.size() << "\n";
    
    
    //Next step is to take the split records by individual, and extract the data in the indiv structure
    vector<indiv> all_indiv;
    CompileAllIndiv(indiv_pat_records,indiv_hcw_records,all_indiv);
    
    if (p.verb==1) {
        WriteIndividualDetails(all_indiv);
    }
    
    cout << "Generate clusters\n";
    //Find networks.  First step - group individuals
    vector< vector<int> > clusters;
    GenerateClusters(all_indiv,clusters);
    if (p.verb==1) {
        WriteClustersRaw(clusters);
    }
    
    vector< vector<int> > clusters_index;
    GenerateClustersIndex(all_indiv,clusters,clusters_index);
    if (p.verb==1) {
        WriteClustersIndexed(clusters);
    }

    cout << "Match to intervention details\n";

    MatchTreatmentNetworks(p,clusters_index,coord1,coord2,interventions,exposures,reductions12,reductions3,all_indiv,rgen);

    cout << "Link transmission clusters\n";

    //Find the links within clusters and probabilties.  We want to work out the cumulative probabilties.
    LinkTransmissionClusters (all_indiv,clusters_index);
    cout << "Write processed information\n";

    //Print cluster information
    if (p.verb==1) {
        WriteClustersProcessed (all_indiv,clusters_index);
    }
    
    //Find number of nosocomial infections among patients and HCWs.
    FindReductionInCases (all_indiv,clusters_index);
    
    
    return 0;
    //Import hospital worker data
}
	
