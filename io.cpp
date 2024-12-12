#include "process_simulation.h"
#include "io.h"
#include "utilities.h"
#include <string>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
    p.sim=1;
	p.verb=0;
    p.intervention_length=10;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--sim")==0) {
			x++;
			p.sim=atoi(argv[x]);
		} else if (p_switch.compare("--verb")==0) {
			x++;
			p.verb=atoi(argv[x]);
        } else if (p_switch.compare("--intervention")==0) {
            x++;
            p.intervention_length=atoi(argv[x]);

        } else {
			cout << "Incorrect usage " << argv[x] << "\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ReadPatientData (run_params& p, vector<precord>& pat_records) {
    ostringstream convert;
    convert << p.sim;
    string temp=convert.str();
    string pat_file_name = "/Users/christopher.illingworth/Documents/Coronavirus/Simulation/StephData/patient_data/pat_data_full_stay"+temp+".csv";
    ifstream pat_file;
    pat_file.open(pat_file_name.c_str());
    cout << pat_file_name.c_str() << "\n";
    string str;
    int i=-1;
    while (getline(pat_file,str)) {
        i++;
        if (i>0) {
            precord pp;
            RemovePunc(str);
            vector<string> subs;
            SplitCommas(str,subs);
            pp.date=subs[1];
            pp.idee=stoi(subs[2],nullptr,10);
            pp.status=subs[3];
            pp.wardno=stoi(subs[5],nullptr,10);
            pp.infected_by_status=subs[9];
            if (subs[10].compare("ED")!=0) {  //Assume we can't diagnose these...
                pp.infected_by_number=-1;
            }
            if (subs[10].compare("NA")!=0) {
                pp.infected_by_number=stoi(subs[10],nullptr,10);
            } else {
                pp.infected_by_number=-1;
            }
            if (subs[11].compare(0,4,"TRUE")==0) {
                pp.detected=1;
            } else {
                pp.detected=0;
            }
            /*if (pp.idee==48666) {
                cout << pp.date << " " << pp.idee << " " << pp.status << " " << pp.wardno << " " << pp.infected_by_status << " " << pp.infected_by_number << "\n";
            }*/
            pat_records.push_back(pp);
        }
    }
}

void ReadHCWData (run_params& p, vector<precord>& hcw_records) {
    ostringstream convert;
    convert << p.sim;
    string temp=convert.str();
    string hcw_file_name = "/Users/christopher.illingworth/Documents/Coronavirus/Simulation/StephData/hcw_data/hcw_data_inf_30b10a"+temp+".csv";
    ifstream hcw_file;
    hcw_file.open(hcw_file_name.c_str());
    cout << hcw_file_name.c_str() << "\n";
    string str;
    int i=-1;
    while (getline(hcw_file,str)) {
        i++;
        if (i>0) {
            precord pp;
            RemovePunc(str);
            vector<string> subs;
            SplitCommas(str,subs);
            pp.date=subs[1];
            pp.idee=stoi(subs[2],nullptr,10)+1000000;
            pp.status=subs[3];
            pp.wardno=-1;
            pp.infected_by_status=subs[6];
            if (subs[7].compare("ED")!=0) {  //Assume we can't diagnose these...
                pp.infected_by_number=-1;
            }
            if (subs[7].compare("Comm")!=0) {  //Assume we can't diagnose these...
                pp.infected_by_number=-1;
            }
            if (subs[7].compare("NA")!=0&&subs[7].compare("")!=0) {
                pp.infected_by_number=stoi(subs[7],nullptr,10);
            } else {
                pp.infected_by_number=-1;
            }
            if (subs[8].compare(0,4,"TRUE")==0) {
                pp.detected=1;
            } else {
                pp.detected=0;
            }
            /*if (pp.idee == 1004529) {
                cout << str << "\n";
                cout << pp.date << " " << pp.idee << " " << pp.status << " " << pp.wardno << " " << pp.infected_by_status << " " << pp.infected_by_number << "\n";
            }*/
            hcw_records.push_back(pp);
        }
    }

}

    
void RemovePunc(string& str) {
    //Edit string to remove "
    vector<int> rem;
    for (unsigned int j=0;j<str.size();j++) {
        if (str[j]=='"') {
            rem.push_back(j);
        }
    }
    reverse(rem.begin(),rem.end());
    for (unsigned int j=0;j<rem.size();j++) {
        str.erase(str.begin()+rem[j]);
    }
}

void SplitCommas(const string str, vector<string>& subs) {
    stringstream ss(str);
    while (ss.good() ) {
        string sr;
        getline(ss,sr,',');
        subs.push_back(sr);
    }
}

void GenerateYMD (vector<precord>& pat_records) {
    char delim='-';
    vector<int> ymd;
    for (int i=0;i<pat_records.size();i++) {
        MakeYMD(pat_records[i].date,delim,ymd);
        pat_records[i].day_no=ConvertYMD(ymd);
        //cout << ymd[2] << " " << ymd[1] << " " << ymd[0] << " " << pat_records[i].day_no << "\n";
    }
    for (int i=0;i<pat_records.size();i++) {
        
    }
}

void MakeYMD (const string subs, char delim, vector<int>& ymd) {
    ymd.clear();
    stringstream sss(subs);
    while (sss.good()) {
        string sr;
        getline(sss,sr,delim);
        ymd.push_back(atoi(sr.c_str()));
    }
}

int ConvertYMD (vector<int>& ymd) {
    int x=0;
    if (ymd[0]>1999) {
        ymd[0]=ymd[0]-2000;
    }
    int years=ymd[0]-20;
    x=x+(years*365);
    //Account for previous leap years
    if (years>0) {
        int prev_lp=1+((years-1)/4);
        x=x+prev_lp;
    }
    //Previous months of this year
    if (ymd[1]>1) {
        x=x+31;
    }
    if (ymd[1]>2) {
        x=x+28;
        if (years%4==0) {
            x=x+1;
        }
    }
    if (ymd[1]>3) {
        x=x+31;
    }
    if (ymd[1]>4) {
        x=x+30;
    }
    if (ymd[1]>5) {
        x=x+31;
    }
    if (ymd[1]>6) {
        x=x+30;
    }
    if (ymd[1]>7) {
        x=x+31;
    }
    if (ymd[1]>8) {
        x=x+31;
    }
    if (ymd[1]>9) {
        x=x+30;
    }
    if (ymd[1]>10) {
        x=x+31;
    }
    if (ymd[1]>11) {
        x=x+30;
    }
    //Account for day
    x=x+ymd[2]-1;
    return(x);
}


void ReadInterventionData (run_params& p, vector< vector<int> >& interventions) {
    ostringstream convert;
    convert << p.sim;
    string temp=convert.str();
    ostringstream intervene;
    intervene << p.intervention_length;
    string temps=intervene.str();
    string int_file_name = "/Users/christopher.illingworth/Documents/Coronavirus/Simulation/UDCA_Project/WriteupCI/Data/Ward_detection_windows_patient_only"+temp+"_"+temps+".dat";
    cout << int_file_name << "\n";
    ifstream int_file;
    int_file.open(int_file_name.c_str());
    int n;
    for (int i=0;i<100000000;i++) {
        vector<int> ini;
        if (!(int_file >> n)) break;
        ini.push_back(n);
        if (!(int_file >> n)) break;
        ini.push_back(n);
        if (!(int_file >> n)) break;
        ini.push_back(n);
        interventions.push_back(ini);
    }
}

void WriteIndividualDetails(const vector<indiv>& all_indiv) {
    ofstream indiv_file;
    indiv_file.open("Individual_details.dat");
    for (int i=0;i<all_indiv.size();i++) {
        indiv_file << all_indiv[i].index << " " << all_indiv[i].hcw << " " << all_indiv[i].name << " " << all_indiv[i].hcw << " " << all_indiv[i].day_i << " " << all_indiv[i].day_d << " " << all_indiv[i].infected_by << "\n";
        /*for (int j=0;j<all_indiv[i].ward.size();j++) {
         cout << all_indiv[i].ward[j] << " ";
         }*/
        
        //cout << "\n";
    }
    indiv_file.close();
}

void WriteClustersRaw (const vector< vector<int> >& clusters) {
    ofstream cl_file;
    cl_file.open("Clusters.dat");
    for (int c=0;c<clusters.size();c++) {
        cl_file << c << " ";
        for (int d=0;d<clusters[c].size();d++) {
            cl_file << clusters[c][d] << " ";
        }
        cl_file << "\n";
    }
    cl_file.close();
}

void WriteClustersIndexed (const vector< vector<int> >& clusters_index) {
    ofstream cl2_file;
    cl2_file.open("ClustersIndex.dat");
    for (int c=0;c<clusters_index.size();c++) {
        for (int d=0;d<clusters_index[c].size();d++) {
            cl2_file << clusters_index[c][d] << " ";
        }
        cl2_file << "\n";
    }
    cl2_file.close();
}

void WriteClustersProcessed (const vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index) {
    ofstream cl3_file;
    cl3_file.open("Clusters_Processed.dat");
    for (int c=0;c<clusters_index.size();c++) {
        cl3_file << "Cluster " << c << "\n";
        for (int d=0;d<clusters_index[c].size();d++) {
            cl3_file << "Indiv " << clusters_index[c][d] << " " << all_indiv[clusters_index[c][d]].name << " ";
            if (all_indiv[clusters_index[c][d]].hcw==1) {
                cl3_file << "HCW ";
            } else {
                cl3_file << "PAT ";
            }
            cl3_file << all_indiv[clusters_index[c][d]].infected_by << " " << all_indiv[clusters_index[c][d]].p_infected << " " << all_indiv[clusters_index[c][d]].p_case << "\n";
        }
    }
}
