using namespace std;
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

struct run_params {
	string method; //Describe the calculation the simulation will perform.  This will take text input and will be essential for any calculation.  Methods are as follows:
    int sim; //Number of simulation
    int verb;
    int intervention_length; //In days.  Default is 10 days
};


struct precord {
    string date;
    int idee;
    int wardno;
    string status;
    string infected_by_status;
    int infected_by_number;
    int day_no;  //Days since 1-1-2020
    int detected;
};

struct indiv {
    int hcw; //Flag if HCW
	int index; //Number increases for each person
    int name;  //From Steph's simulation data
	int infected_by; //Who infected them
    double p_infected; //Nominal probability of being infected
    double p_case; //Probability of being infected accounting for network structure
	int day_i; //Day infected
    int day_d; //Day detected if relevant
    vector<int> ward;
    vector<int> ward_time;
};

struct indivh {  //To encompass differential ACE2 levels.  Don't need sequence model here
	int index;
	int day_i; //Day infected
	int day_symp; //Day symptomatic
	//vector<double> ace2; //Time-dependent expression depending upon treatment
};
