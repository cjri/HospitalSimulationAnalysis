#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void SetupDistributions (vector<double>& exposures, vector<double>& ace_basic, vector<double>& ace12, vector<double>& ace3);
void NormaliseDistributions (double mo, vector<double>& ace_basic, vector<double>& ace12, vector<double>& ace3);
void MakeConditionalDistribution (const vector<double>& exposures, const vector<double>& ace_basic, vector< vector<double> >& condo);
void FindPartialProbabilities (const vector<double>& exposures, const vector<double>& ace_basic, const vector< vector<double> >& condo, vector<double>& coord1, vector< vector<double> >& coord2);
void SampleProbabilities (int& choice1, int& choice2, const vector<double>& exposures, const vector<double>& coord1, const vector< vector<double> >& coord2, gsl_rng *rgen);
void GetReductions (const vector<double>& exposures, const vector<double>& ace_basic, const vector<double>& ace12, const vector<double>& ace3, vector< vector<double> >& reductions12, vector< vector<double> >& reductions3);


void FindUniqueIndividuals (const vector<precord>& pat_records, vector<int>& patients);
void SplitRecordsByIndividual (const vector<int>& patients, const vector<precord>& pat_records, vector< vector<precord> >& indiv_pat_records);
void CompileAllIndiv(const vector< vector<precord> >& indiv_pat_records, vector< vector<precord> >& indiv_hcw_records, vector<indiv>& all_indiv);
void GenerateClusters(const vector<indiv>& all_indiv, vector< vector<int> >& clusters);
void GenerateClustersIndex(const vector<indiv>& all_indiv, const vector< vector<int> >& clusters, vector< vector<int> >& clusters_index);
void MatchTreatmentNetworks (run_params& p, const vector< vector<int> >& clusters_index, const vector<double>& coord1, const vector< vector<double> >& coord2, const vector< vector<int> >& interventions, const vector<double>& exposures, const vector< vector<double> > reductions12, const vector< vector<double> > reductions3, vector<indiv>& all_indiv, gsl_rng *rgen);
void LinkTransmissionClusters (vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index);
void FindReductionInCases (const vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index);
