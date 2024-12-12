#include "process_simulation.h"
#include "io.h"
#include "utilities.h"
#include <string>

void SetupDistributions (vector<double>& exposures, vector<double>& ace_basic, vector<double>& ace12, vector<double>& ace3) {
    //Exposure values
    for (int i=0;i<1000;i++) {
        double p=(i+1.)/1001.;
        double v=gsl_cdf_gamma_Pinv(p,0.15580249308940242,2.7972762585308675);
        exposures.push_back(v);
    }
    //ACE2 relative expression
    double mo=0;
    for (int i=0;i<1000;i++) {
        double p=(i+1.)/1001.;
        double v=gsl_cdf_gamma_Pinv(p,3.3081986097097356,0.03436092435733757);
        ace_basic.push_back(v);
        mo=mo+v;
    }
    mo=mo/1000;
    for (int i=0;i<1000;i++) {
        double p=(i+1.)/1001.;
        double v=gsl_cdf_gamma_Pinv(p,2.233325253871696,0.022225496524666207);
        ace12.push_back(v);
    }
    for (int i=0;i<1000;i++) {
        double p=(i+1.)/1001.;
        double v=gsl_cdf_gamma_Pinv(p,1.081675307659759,0.031041321249140563);
        ace3.push_back(v);
    }
    NormaliseDistributions (mo,ace_basic,ace12,ace3);
}

void NormaliseDistributions (double mo, vector<double>& ace_basic, vector<double>& ace12, vector<double>& ace3) {
    for (int i=0;i<ace_basic.size();i++) {
        ace_basic[i]=ace_basic[i]/mo;
        ace12[i]=ace12[i]/mo;
        ace3[i]=ace3[i]/mo;
    }
}

void MakeConditionalDistribution (const vector<double>& exposures, const vector<double>& ace_basic, vector< vector<double> >& condo) {
    double ctot=0;
    for (int i=0;i<exposures.size();i++) {
        vector<double> c;
        for (int j=0;j<ace_basic.size();j++) {
            double cv=1-exp(-exposures[i]*ace_basic[j]);
            c.push_back(cv);
            ctot=ctot+cv;
        }
        condo.push_back(c);
    }
    for (int i=0;i<exposures.size();i++) {
        for (int j=0;j<ace_basic.size();j++) {
            condo[i][j]=condo[i][j]/ctot;
        }
    }
}

void FindPartialProbabilities (const vector<double>& exposures, const vector<double>& ace_basic, const vector< vector<double> >& condo, vector<double>& coord1, vector< vector<double> >& coord2){
    for (int i=0;i<exposures.size();i++) {
        double etot=0;
        for (int j=0;j<ace_basic.size();j++) {
            etot=etot+condo[i][j];
        }
        coord1.push_back(etot);
    }
    coord2=condo;
    for (int i=0;i<exposures.size();i++) {
        double rtot=0;
        for (int j=0;j<ace_basic.size();j++) {
            rtot=rtot+condo[i][j];
        }
        for (int j=0;j<ace_basic.size();j++) {
            coord2[i][j]=coord2[i][j]/rtot;
        }
    }
}

void SampleProbabilities (int& choice1, int& choice2, const vector<double>& exposures, const vector<double>& coord1, const vector< vector<double> >& coord2, gsl_rng *rgen) {
    double pexposure[coord1.size()];
    for (int i=0;i<coord1.size();i++) {
        pexposure[i]=coord1[i];
    }
    unsigned int N=1;
    unsigned int draw[exposures.size()];
    gsl_ran_multinomial(rgen,coord1.size(),N,pexposure,draw);
        for (int i=0;i<coord1.size();i++) {
        if (draw[i]==1) {
            choice1=i;
        }
    }
    
    double pexpression[coord2[choice1].size()];
    for (int i=0;i<coord2[choice1].size();i++) {
        pexpression[i]=coord2[choice1][i];
    }
    unsigned int draw2[exposures.size()];
    gsl_ran_multinomial(rgen,coord2.size(),N,pexpression,draw2);
    for (int i=0;i<coord2[choice1].size();i++) {
        if (draw2[i]==1) {
            choice2=i;
        }
    }
}

void GetReductions (const vector<double>& exposures, const vector<double>& ace_basic, const vector<double>& ace12, const vector<double>& ace3, vector< vector<double> >& reductions12, vector< vector<double> >& reductions3) {
    for (int i=0;i<exposures.size();i++) {
        vector<double> r12;
        vector<double> r3;
        for (int j=0;j<ace_basic.size();j++) {
            double val12=(1-exp(-exposures[i]*ace12[j]))/(1-exp(-exposures[i]*ace_basic[j]));
            double val3=(1-exp(-exposures[i]*ace3[j]))/(1-exp(-exposures[i]*ace_basic[j]));
            r12.push_back(val12);
            r3.push_back(val3);
        }
        reductions12.push_back(r12);
        reductions3.push_back(r3);
    }
}

void FindUniqueIndividuals (const vector<precord>& pat_records, vector<int>& patients) {
    for (int i=0;i<pat_records.size();i++) {
        patients.push_back(pat_records[i].idee);
    }
    sort(patients.begin(),patients.end());
    patients.erase(unique(patients.begin(),patients.end()),patients.end());
}

void SplitRecordsByIndividual (const vector<int>& patients, const vector<precord>& pat_records, vector< vector<precord> >& indiv_pat_records) {
    for (int i=0;i<patients.size();i++) {
        vector<precord> ipat;
        for (int j=0;j<pat_records.size();j++) {
            if (pat_records[j].idee==patients[i]) {
                ipat.push_back(pat_records[j]);
            }
        }
        indiv_pat_records.push_back(ipat);
    }
}

void CompileAllIndiv(const vector< vector<precord> >& indiv_pat_records, vector< vector<precord> >& indiv_hcw_records, vector<indiv>& all_indiv) {
    int index=0;
    for (int i=0;i<indiv_pat_records.size();i++) {
        indiv ii;
        ii.hcw=0;
        ii.index=index;
        index++;
        ii.day_i=-1;
        ii.day_d=-1;
        ii.name=indiv_pat_records[i][0].idee;
        ii.p_infected=0;
        ii.infected_by=-1;
        for (int j=0;j<indiv_pat_records[i].size();j++) {
            /*if (ii.name==48666) {
                cout << "Read " << indiv_pat_records[i][j].status << " " << indiv_pat_records[i][j].infected_by_status << "\n";
            }*/

            ii.ward.push_back(indiv_pat_records[i][j].wardno);
            ii.ward_time.push_back(indiv_pat_records[i][j].day_no);
            if (ii.day_i==-1&&indiv_pat_records[i][j].status.compare(0,8,"INFECTED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_pat_records[i][j].day_no;
                ii.infected_by=indiv_pat_records[i][j].infected_by_number;
            }
            if (ii.day_i==-1&&indiv_pat_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_pat_records[i][j].infected_by_number!=-1) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_pat_records[i][j].day_no;
                ii.infected_by=indiv_pat_records[i][j].infected_by_number;
            }
            if (ii.day_i==-1&&indiv_pat_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_pat_records[i][j].infected_by_status.compare(0,2,"ED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_pat_records[i][j].day_no;
            }
            if (ii.day_i==-1&&indiv_pat_records[i][j].status.compare(0,9,"RECOVERED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_pat_records[i][j].day_no;
                ii.infected_by=indiv_pat_records[i][j].infected_by_number;
            }
            if (ii.day_d==-1&&indiv_pat_records[i][j].detected==1) {
                ii.day_d=indiv_pat_records[i][j].day_no;
            }
        }
        if (ii.day_d==-1) {  //Some cases exposed but source not recorded.
            for (int j=0;j<indiv_pat_records[i].size();j++) {
                if (indiv_pat_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_pat_records[i][j].wardno!=-1) {
                    ii.p_infected=1;
                    ii.p_case=1;
                    ii.day_i=indiv_pat_records[i][j].day_no;
                }
            }
        }

        /*if (ii.name==48666) {
            cout << "Add case " << ii.p_infected << " " << ii.infected_by  << "\n";
        }*/

        all_indiv.push_back(ii);
    }

    for (int i=0;i<indiv_hcw_records.size();i++) {
        indiv ii;
        ii.hcw=1;
        ii.index=index;
        index++;
        ii.day_i=-1;
        ii.day_d=-1;
        ii.name=indiv_hcw_records[i][0].idee;
        ii.p_infected=0;
        ii.infected_by=-1;
        for (int j=0;j<indiv_hcw_records[i].size();j++) {
            /*if (ii.name==1004529) {
                cout << "Read " << indiv_hcw_records[i][j].status << " " << indiv_hcw_records[i][j].infected_by_status << " " << indiv_hcw_records[i][j].day_no << "\n";
            }*/

            ii.ward.push_back(indiv_hcw_records[i][j].wardno);
            ii.ward_time.push_back(indiv_hcw_records[i][j].day_no);
            if (ii.day_i==-1&&indiv_hcw_records[i][j].status.compare(0,8,"INFECTED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_hcw_records[i][j].day_no;
                ii.infected_by=indiv_hcw_records[i][j].infected_by_number;
            }
            if (ii.day_i==-1&&indiv_hcw_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_hcw_records[i][j].infected_by_number!=-1) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_hcw_records[i][j].day_no;
                ii.infected_by=indiv_hcw_records[i][j].infected_by_number;
            }
            if (ii.day_i==-1&&indiv_hcw_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_hcw_records[i][j].infected_by_status.compare(0,2,"ED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_hcw_records[i][j].day_no;
            }
            if (ii.day_i==-1&&indiv_hcw_records[i][j].status.compare(0,7,"EXPOSED")==0&&indiv_hcw_records[i][j].infected_by_status.compare(0,4,"Comm")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_hcw_records[i][j].day_no;
            }
            if (ii.day_i==-1&&indiv_hcw_records[i][j].status.compare(0,9,"RECOVERED")==0) {
                ii.p_infected=1;
                ii.p_case=1;
                ii.day_i=indiv_hcw_records[i][j].day_no;
                ii.infected_by=indiv_hcw_records[i][j].infected_by_number;
            }
            if (ii.day_d==-1&&indiv_hcw_records[i][j].detected==1) {
                ii.day_d=indiv_hcw_records[i][j].day_no;
            }
        }
        all_indiv.push_back(ii);
    }
}



void GenerateClusters(const vector<indiv>& all_indiv, vector< vector<int> >& clusters) {
    vector<int> networked;
    for (int i=0;i<all_indiv.size();i++) {
        networked.push_back(0);
    }
    
    for (int i=0;i<all_indiv.size();i++) {  //Find introductions
        if (networked[i]==0&&all_indiv[i].infected_by==-1) {
            vector<int> c;
            /*if (i==194) {
                cout << "Start with " << i << " " << all_indiv[i].name << "\n";
            }*/
            c.push_back(all_indiv[i].name);
            /*if (all_indiv[i].name==6934) {
                cout << "Add 6934 " << i << " " << all_indiv[i].hcw << "\n";
            }*/
            networked[i]=1;
            clusters.push_back(c);
        }
    }
    
    for (int c=0;c<clusters.size();c++) {  //Build clusters in turn
        int add=1;
        while (add>0) {
            //Go through the cluster.  Find anyone infected
            add=0;
            for (int i=0;i<all_indiv.size();i++) {
                if (networked[i]==0) {
                    /*if (i==194) {
                        cout << "Check 194\n";
                    }*/
                    for (int d=0;d<clusters[c].size();d++) {
                        if (all_indiv[i].infected_by==clusters[c][d]) {
                            add=1;
                            /*if (all_indiv[i].name==6934) {
                                cout << "Add 6934 " << i << " " << all_indiv[i].hcw << "\n";
                            }*/
                            clusters[c].push_back(all_indiv[i].name);
                            networked[i]=1;
                        }
                    }
                }
            }
        }
    }
}


void GenerateClustersIndex(const vector<indiv>& all_indiv, const vector< vector<int> >& clusters, vector< vector<int> >& clusters_index) {
    clusters_index=clusters;
    for (int i=0;i<clusters.size();i++) {
        for (int j=0;j<clusters[i].size();j++) {
            clusters_index[i][j]=-1;
            for (int k=0;k<all_indiv.size();k++) {
                if (all_indiv[k].name==clusters[i][j]) {
                    clusters_index[i][j]=k;
                    break;
                }
            }
            if (clusters_index[i][j]==-1) {
                cout << "Error " << i << " " << j << " " << all_indiv.size() << " " << clusters[i][j] << "\n";
            }
        }
    }
}
    
    
void MatchTreatmentNetworks (run_params& p, const vector< vector<int> >& clusters_index, const vector<double>& coord1, const vector< vector<double> >& coord2, const vector< vector<int> >& interventions, const vector<double>& exposures, const vector< vector<double> > reductions12, const vector< vector<double> > reductions3, vector<indiv>& all_indiv, gsl_rng *rgen) {
    int choice1;
    int choice2;
    for (int c=0;c<clusters_index.size();c++) {
        //cout << "Cluster " << c << "\n";
        for (int d=0;d<clusters_index[c].size();d++) {
            //cout << "Indiv " << clusters_index[c][d] << " " << all_indiv[clusters_index[c][d]].name << " " << all_indiv[clusters_index[c][d]].hcw << "\n";
            int match=0;
            vector<int> interventions_match;
            if (all_indiv[clusters_index[c][d]].hcw==0) {
                //Is there a day when the ward/day information falls into an intervention?
                for (int i=0;i<all_indiv[clusters_index[c][d]].ward.size();i++) {
                    for (int j=0;j<interventions.size();j++) {
                        match=0;
                        for (int k=0;k<interventions_match.size();k++) {
                            if (interventions_match[k]==j) {
                                match=1;
                            }
                        }
                        if (all_indiv[clusters_index[c][d]].ward[i]==interventions[j][0]) {
                            if (match==0&&all_indiv[clusters_index[c][d]].ward_time[i]>=interventions[j][1]&&all_indiv[clusters_index[c][d]].ward_time[i]<=interventions[j][2]) {
                                interventions_match.push_back(j);
                                //Patient was on a ward with an intervention.
                                
                                //Therefore, find the first date on which the drug was given.
                                /*cout << "Match " << all_indiv[clusters_index[c][d]].ward[i] << " " << all_indiv[clusters_index[c][d]].ward_time[i] << " ";
                                for (int k=0;k<3;k++) {
                                    cout << interventions[j][k] << " ";
                                }
                                cout << "\n";*/
                                //Check for being infected rather than simply community infection
                                if (all_indiv[clusters_index[c][d]].infected_by!=-1 && all_indiv[clusters_index[c][d]].day_i != -1) {//Infected by someone else i.e. not community
                                    if (c==685) {
                                        cout << "Day of infection " << all_indiv[clusters_index[c][d]].day_i << " by " << all_indiv[clusters_index[c][d]].infected_by << " baseline " << all_indiv[clusters_index[c][d]].ward_time[i] << "\n";
                                    }
                                    int tti=all_indiv[clusters_index[c][d]].day_i - all_indiv[clusters_index[c][d]].ward_time[i];
                                    //cout << "Time treatment to infection " << tti << "\n";
                                    //Calculate reduction in the probability of infection if the time of infection was within the first day of
                                    //being on the ward during the intervention
                                    if (tti>0&&tti<=2) {
                                        SampleProbabilities (choice1,choice2,exposures,coord1,coord2,rgen);
                                        //cout << "Reduction " << reductions12[choice1][choice2] << "\n";
                                        all_indiv[clusters_index[c][d]].p_infected=reductions12[choice1][choice2];
                                        
                                        //Note: Could alter the reduced value here...
                                        
                                    }
                                    if (tti>2&&tti<=p.intervention_length) {
                                        SampleProbabilities (choice1,choice2,exposures,coord1,coord2,rgen);
                                        //cout << "Reduction " << reductions3[choice1][choice2] << "\n";
                                        all_indiv[clusters_index[c][d]].p_infected=reductions3[choice1][choice2];

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void LinkTransmissionClusters (vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index) {
    for (int c=0;c<clusters_index.size();c++) {
        //cout << "Cluster " << c << " " << clusters_index.size() << "\n";
        vector<int> considered;
        vector<int> current;
        vector<int> current_index;
        vector<int> next;
        vector<int> next_index;
        int add=1;
        for (int d=0;d<clusters_index[c].size();d++) {
            considered.push_back(0);
            if (all_indiv[clusters_index[c][d]].infected_by==-1) {
                current.push_back(all_indiv[clusters_index[c][d]].name);
                current_index.push_back(d);
                considered[d]=1;
            }
        }
        while (add==1) {
            add=0;
            for (int d=0;d<clusters_index[c].size();d++) {
                for (int e=0;e<current.size();e++) {
                    if (considered[d]==0 && considered[e]==1 && all_indiv[clusters_index[c][d]].infected_by==all_indiv[clusters_index[c][current_index[e]]].name) {
                        if (c==685) {
                            cout << "Join " << d << " " << all_indiv[clusters_index[c][d]].name << " infected by " << all_indiv[clusters_index[c][d]].infected_by << " to " << e << " " << all_indiv[clusters_index[c][current_index[e]]].name << " with p(i) " << all_indiv[clusters_index[c][current_index[e]]].p_case << "\n";
                            cout << "From " << all_indiv[clusters_index[c][d]].p_case << " to ";
                        }
                        all_indiv[clusters_index[c][d]].p_case=all_indiv[clusters_index[c][d]].p_case*all_indiv[clusters_index[c][current_index[e]]].p_case;
                        if (c==685) {
                            cout << all_indiv[clusters_index[c][d]].p_case;
                        }
                        all_indiv[clusters_index[c][d]].p_case=all_indiv[clusters_index[c][d]].p_case*all_indiv[clusters_index[c][d]].p_infected;
                        if (c==685) {
                            cout << " via " << all_indiv[clusters_index[c][d]].p_infected << " " << all_indiv[clusters_index[c][d]].p_case << "\n";
                        }
                        next.push_back(all_indiv[clusters_index[c][d]].name);
                        next_index.push_back(d);
                        considered[d]=1;
                        add=1;
                    }
                }
            }
            current=next;
            current_index=next_index;
            next.clear();
            next_index.clear();
        }
        for (int d=0;d<clusters_index[c].size();d++) {
            //all_indiv[clusters_index[c][d]].p_case=all_indiv[clusters_index[c][d]].p_case*all_indiv[clusters_index[c][d]].p_infected;
        }
    }
}
    
void FindReductionInCases (const vector<indiv>& all_indiv, const vector< vector<int> >& clusters_index) {
    double p_infections=0;
    double h_infections=0;
    double p_infections_int=0;
    double h_infections_int=0;
    for (int c=0;c<clusters_index.size();c++) {
        for (int d=0;d<clusters_index[c].size();d++) {
            if (all_indiv[clusters_index[c][d]].infected_by!=-1) {
                if (all_indiv[clusters_index[c][d]].hcw==0) {
                    p_infections++;
                    p_infections_int=p_infections_int+all_indiv[clusters_index[c][d]].p_case;
                } else {
                    h_infections++;
                    h_infections_int=h_infections_int+all_indiv[clusters_index[c][d]].p_case;
                }
            }
        }
    }
    cout << "Patients " << p_infections << " " << p_infections_int << " " << p_infections_int/p_infections << "\n";
    cout << "HCWs " << h_infections << " " << h_infections_int << " " << h_infections_int/h_infections << "\n";
}
