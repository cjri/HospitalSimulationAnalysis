#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>


void GetOptions (run_params& p, int argc, const char **argv);
void ReadPatientData (run_params& p, vector<precord>& pat_records);
void ReadHCWData (run_params& p, vector<precord>& hcw_records);
void GenerateYMD (vector<precord>& pat_records);
void MakeYMD (const string subs, char delim, vector<int>& ymd);
int ConvertYMD (vector<int>& ymd);

void RemovePunc(string& str);
void SplitCommas(const string str, vector<string>& subs);
void ReadInterventionData (run_params& p, vector< vector<int> >& interventions);

void WriteIndividualDetails(const vector<indiv>& all_indiv);
void WriteClustersRaw (const vector< vector<int> >& clusters);
void WriteClustersIndexed (const vector< vector<int> >& clusters_index);
void WriteClustersProcessed (const vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index);

