#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

typedef struct {
	
	string type;
    string sel;
	double F;
	double uF;
	double diffF;
	double udiffF;
	double eb;
	double ueb;
	double diffeb;
    double nomdiffeb;
    double unomdiffeb;
	double udiffeb;
	double xs;
	double uxs;
	double diffxs;
	double udiffxs;
    double SF;
    double diffSF;

} systematic;

double max (double a, double b) {
	if (fabs(a)>fabs(b)) return fabs(a);
	return fabs(b);
}

double maxPDF (double a, double b) {
	if (a>0 && a>b) return a;
	else if (b>0 && a<b) return b;
	return 0;
}
	
double round(double depth, double n)
{
    double d;
    int    i;
	
    /* rescale 123.45678 to 12345.678 */ 
    d = fabs(n * depth);
    /* round off: 12345.678 + 0.5 = 12346.178 -> 12346 */ 
    i = d + 0.5;
    /* restore to its original scale: 12346 -> 123.46 */
    d = (double)i / depth;
	
    if (n < 0)
        return -d;
    return d;
}

void printSystematics(std::string type, std::map<std::string, std::vector<double> > systematics, std::map<std::string,double> bias, int depth) {
    cout <<  "----------------- Systematics on "<< type <<" -----------------" << endl;
	double sumSyst=0;
	
	for (std::map<std::string, std::vector<double> >::iterator it = systematics.begin(); it != systematics.end(); ++it) {
		
		if (it->second.size() > 2)
			cout << "Error, " << it->first << " has more than two values. This is unexpected!" << endl;
		
		else {
            
            double sys=max(fabs(it->second[0]),fabs(it->second[1]));
            
            //cout << max(it->second[0],it->second[1]) << endl;
            
            //if (bias.find(it->first) != bias.end()) {
            //cout << sys << " " << bias[it->first] << endl;
            //    sys = sqrt(pow(sys,2)+pow(bias[it->first],2));
            //}
            
            //if (sys > 1./(depth*10) ) {
            if (it->second.size() ==1) (it->second).push_back(0);
			
            cout << "+> Systematic " << it->first << " = " << round(depth,sys) << endl;
			
            sumSyst+=pow(sys,2);
            
			//}
		}	
		
	}
	
	cout << endl << "+> TOTAL: " << round(depth,sqrt(sumSyst)) << endl;
    
    ofstream f(("systematics/tmp/tmp_total_"+type).c_str(), ios::out | ios::trunc );
    
    f << round(1000,sqrt(sumSyst)/100);
    
    f.close();
    
}

void printAssymSystematics(std::string type, std::map<std::string, std::vector<double> > systematics, std::map<std::string,double> bias, int depth) {
    cout <<  "----------------- Asymmetric Systematics on "<< type <<" -----------------" << endl;
	double sumSyst_up=0;
	double sumSyst_down=0;
    	
	for (std::map<std::string, std::vector<double> >::iterator it = systematics.begin(); it != systematics.end(); ++it) {
		
        //if (it->second[0] == 0 && it->second[1] == 0) continue;
        
		if (it->second.size() > 2)
			cout << "Error, " << it->first << " has more than two values. This is unexpected!" << endl;
		
		else {
            
            double sys_up = 0;
            double sys_down = 0;
            
            //if (it->second[0] > 0 && it->second[1] > 0) cout << it->first << " ERROR both > 0 " << it->second[0] << " " << it->second[1] << endl;
            //if (it->second[0] < 0 && it->second[1] < 0) cout << it->first << " ERROR both < 0 " << it->second[0] << " " << it->second[1] << endl;
            
            /*if (it->second[0] > 0 && it->second[1] < 0) {
                sys_up = it->second[0];
                sys_down = it->second[1];
            }

            else if (it->second[0] < 0 && it->second[1] > 0) {
                sys_up = it->second[1];
                sys_down = it->second[0];
            }
            
            else if (fabs(it->second[0]) != 0 && fabs(it->second[1]) != 0) {
                cout << it->first << " ERROR we don't have a + and a - " << it->second[0] << " " << it->second[1] << endl;
                if (fabs(it->second[0]) < fabs(it->second[1])) {
                    sys_up=fabs(it->second[0]);
                    sys_down=it->second[1];
                } else {
                    sys_up=it->second[0];
                    sys_down=fabs(it->second[1]);
                }
            }*/
            
            sys_up=it->second[0];
            sys_down=it->second[1];
            //cout << max(it->second[0],it->second[1]) << endl;
            
            //if (bias.find(it->first) != bias.end()) {
            //cout << sys << " " << bias[it->first] << endl;
            //    sys = sqrt(pow(sys,2)+pow(bias[it->first],2));
            //}
            
            //if (sys > 1./(depth*10) ) {
            //if (it->second.size() ==1) (it->second).push_back(0);
			
            cout << "+> Systematic " << it->first << " = " << round(depth,sys_up) << " (up) " << round(depth,sys_down) << " (down)" << endl;
            //cout << "+> Systematic " << it->first << " = " << it->second[0] << " (up) " << it->second[1] << " (down)" << endl;
			
            sumSyst_up+=pow(sys_up,2);
            sumSyst_down+=pow(sys_down,2);
            
			//}
		}	
		
	}
	
	cout << endl << "+> TOTAL: " << round(depth,sqrt(sumSyst_up)) << " (up) " << round(depth,sqrt(sumSyst_down)) << " (down) " <<  endl;
}

int main (int argc, char* argv[]) {
	
	std::map<std::string, std::vector<double> > systematics_f;
	std::map<std::string, std::vector<double> > systematics_eb;
	std::map<std::string, std::vector<double> > systematics_xs;
	std::map<std::string, std::vector<double> > systematics_SF;
    std::map<std::string, double> bias_eb;
    
    std::vector<double> pdf_up;
    std::vector<double> pdf_down;
    
    
    std::map<std::string, std::vector<double> > systematics_xs_ass;
    
	if (argc < 2) exit(1);
    
    string chan="";
    
    if (argc > 2) chan=argv[2];
	
	string fileName = (string)argv[1];
	
	fstream inFile(fileName.c_str(), ios::in);
	
	vector<systematic> systematics;
	
	cout << "Reading file " << fileName << endl;
    
    double nomdiffeb = -1;
    double unomdiffeb = -1;
    
    double nomxs = -1;
	
	while (!inFile.eof()) {
		
		systematic sys;
		
		inFile >> sys.type >> sys.sel >> sys.F >> sys.uF >> sys.diffF >> sys.udiffF >> sys.eb >> sys.ueb >> sys.diffeb >> sys.udiffeb >> sys.xs >> sys.uxs >> sys.diffxs >> sys.udiffxs >> sys.SF >> sys.diffSF;
        
        if (sys.type.find("nominal") != string::npos) {
            nomdiffeb=sys.eb;
            unomdiffeb=sys.ueb;
            //cout << sys.type << " " << nomdiffeb << endl;
            if (nomxs == -1)
                nomxs = sys.xs;
        } else if (sys.type.find("Data") == string::npos){
            
            //cout << sys.type << endl;
            
            /*sys.diffF = sys.diffF;
            sys.diffeb = sys.diffeb;
            sys.diffxs = sys.diffxs;
            */
            sys.SF = sys.SF*100;
            sys.diffSF = sys.diffSF*100;
            
            sys.nomdiffeb=nomdiffeb;
            sys.unomdiffeb=unomdiffeb;
            
            //cout << sys.type << " " << sys.nomdiffeb << endl;

            systematics.push_back(sys);
        }
 		
	}
	
	////////////////////
	// calcing syst unc
	////////////////////
	
	//for PDF
	double sumplusF = 0;
	double summinF = 0;
	double sumplusE = 0;
	double summinE = 0;
	double sumplusXS = 0;
	double summinXS = 0;
    double summinSF = 0;
    double sumplusSF = 0;
	
	//for Backgrounds
    std::map<std::string, double> map_scaleXS_eb;
    std::map<std::string, double> map_scaleXS_SF;
    std::map<std::string, double> map_scaleXS_F;
    std::map<std::string, double> map_scaleXS_XS;

    std::map<std::string, vector<double> > map_scaleXS_XS_ass;

	//for matching
    std::map<std::string, double> map_match_eb;
    std::map<std::string, double> map_match_SF;
    std::map<std::string, double> map_match_F;
    std::map<std::string, double> map_match_XS;
	
	//for scale
    std::map<std::string, double> map_scale_eb;
    std::map<std::string, double> map_scale_SF;
    std::map<std::string, double> map_scale_F;
    std::map<std::string, double> map_scale_XS;
    
	int nPdf=0;
	for (unsigned int s=0;s<systematics.size();s++) {
	
		systematic sys = systematics[s];
        
        //sys.diffeb = (sys.diffeb/sys.eb)*100;
        
        //if (sys.type.find("JES") == string::npos) continue;
        
        /*string title = "Method bias";
        if (systematics_eb.find(title) == systematics_eb.end()) {
            if (sys.nomdiffeb >= sys.unomdiffeb) {
                systematics_eb[title].push_back(sys.nomdiffeb);
                systematics_SF[title].push_back(sys.nomdiffeb);
            }
            else {
                systematics_eb[title].push_back(sys.unomdiffeb);
                systematics_SF[title].push_back(sys.unomdiffeb);                
            }
            //cout << "lol" << endl;
        }*/
        
		//pdfunc
		if (sys.type.find("PDF-Uncert") != string::npos) {
		
			nPdf++;
			/*if (sys.diffF > 0) {
				sumplusF+=sys.diffF*sys.diffF;
			}	
			if (sys.diffeb < 0) {
				summinF+=sys.diffF*sys.diffF;
			}	

			if (sys.diffeb > 0) {
				sumplusE+=sys.diffeb*sys.diffeb;
			}	
			if (sys.diffeb < 0) {
				summinE+=sys.diffeb*sys.diffeb;
			}	
			
			if (sys.diffxs > 0) {
				sumplusXS+=sys.diffxs*sys.diffxs;
			}	
			if (sys.diffxs < 0) {
				summinXS+=sys.diffxs*sys.diffxs;
			}	*/
			if (nPdf%2 == 0) {
			
				//cout << nPdf << "  down" << endl;
				sumplusF+=sys.diffF*sys.diffF;
				sumplusE+=sys.diffeb*sys.diffeb;
				sumplusXS+=sys.diffxs*sys.diffxs;
				sumplusSF+=sys.diffSF*sys.diffSF;
                
                pdf_up.push_back(sys.diffxs);

			} else {
			
				//cout << nPdf << "  up" << endl;
				summinF+=sys.diffF*sys.diffF;
				summinE+=sys.diffeb*sys.diffeb;
				summinXS+=sys.diffxs*sys.diffxs;
				summinSF+=sys.diffSF*sys.diffSF;
                
                pdf_down.push_back(sys.diffxs);

			}
            
            bias_eb["PDF Uncertainties"]=sys.nomdiffeb;
		}
		// the simple systematics
		else if (sys.type.find("JES") != string::npos) {
			systematics_f["Jet Energy Scale"].push_back(sys.diffF);
			systematics_eb["Jet Energy Scale"].push_back(sys.diffeb);
			systematics_xs["Jet Energy Scale"].push_back(sys.diffxs);
			systematics_SF["Jet Energy Scale"].push_back(sys.diffSF);
            bias_eb["Jet Energy Scale"]=sys.nomdiffeb;
            
            string s="Jet Energy Scale";
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("+") != string::npos || sys.type.find("Up") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;
            else if (sys.type.find("-") != string::npos || sys.type.find("Down") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
        }
        else if (sys.type.find("JER") != string::npos) {
			systematics_f["Jet Energy Resolution"].push_back(sys.diffF);
			systematics_eb["Jet Energy Resolution"].push_back(sys.diffeb);
			systematics_xs["Jet Energy Resolution"].push_back(sys.diffxs);
			systematics_SF["Jet Energy Resolution"].push_back(sys.diffSF);
            bias_eb["Jet Energy Resolution"]=sys.nomdiffeb;
            
            string s="Jet Energy Resolution";
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("+") != string::npos || sys.type.find("Up") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
            else if (sys.type.find("-") != string::npos || sys.type.find("Down") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;

		}
		else if (sys.type.find("PileUp") != string::npos) {
			systematics_f["PileUp"].push_back(sys.diffF);
			systematics_eb["PileUp"].push_back(sys.diffeb);
			systematics_xs["PileUp"].push_back(sys.diffxs);
			systematics_SF["PileUp"].push_back(sys.diffSF);
            bias_eb["PileUp"]=sys.nomdiffeb;
            
            string s="PileUp";
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("+") != string::npos ) {
                systematics_xs_ass[s][0]=sys.diffxs;                // cout << "a " << sys.type << endl;
            }
            else if (sys.type.find("-") != string::npos ) {
                systematics_xs_ass[s][1]=sys.diffxs;             //    cout << "b " <<  sys.type << endl;
            }
            //if (sys.type.find("+") != string::npos || sys.type.find("Up") != string::npos)
            //    cout << "lalalalal" << sys.diffxs << systematics_xs_ass[s][0]<<   endl;
		}
		else if (sys.type.find("UE") != string::npos) {
			systematics_f["Underlying event"].push_back(sys.diffF);
			systematics_eb["Underlying event"].push_back(sys.diffeb);
			systematics_xs["Underlying event"].push_back(sys.diffxs);
			systematics_SF["Underlying event"].push_back(sys.diffSF);
            bias_eb["Underlying event"]=sys.nomdiffeb;
            
            systematics_xs_ass["Underlying event"].push_back(sys.diffxs);
			systematics_xs_ass["Underlying event"].push_back(-sys.diffxs);

		}
		else if (sys.type.find("ISRFSR") != string::npos) {
			systematics_f["ISRFSR"].push_back(sys.diffF);
			systematics_eb["ISRFSR"].push_back(sys.diffeb);
			systematics_xs["ISRFSR"].push_back(sys.diffxs);
			systematics_SF["ISRFSR"].push_back(sys.diffSF);
            bias_eb["ISRFSR"]=sys.nomdiffeb;
		}
		
		else if (sys.type.find("RightReg") != string::npos) {
			systematics_f["Right Region definition"].push_back(sys.diffF);
			systematics_eb["Right Region definition"].push_back(sys.diffeb);
			systematics_xs["Right Region definition"].push_back(sys.diffxs);
			systematics_SF["Right Region definition"].push_back(sys.diffSF);
            bias_eb["Right Region definition"]=sys.nomdiffeb;
           
            string s="Right Region definition";
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("bigger") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
            else if (sys.type.find("smaller") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;
		}

        else if (sys.type.find("TopMass") != string::npos) { // divide by 9 to get the 1GeV/c2 effect
			systematics_f["Top Quark Mass"].push_back(sys.diffF/9.);
			systematics_eb["Top Quark Mass"].push_back(sys.diffeb/9.);
			systematics_xs["Top Quark Mass"].push_back(sys.diffxs/9.);
			systematics_SF["Top Quark Mass"].push_back(sys.diffSF/9.);
            bias_eb["Top Quark Mass"]=sys.nomdiffeb;
            
            string s="Top Quark Mass";
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("+") != string::npos || sys.type.find("181") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs/9.;
            else if (sys.type.find("-") != string::npos || sys.type.find("163") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs/9.;
		}

		else if (sys.type.find("SF") != string::npos) {
			systematics_xs["Muon ID/Eff and Trigger eff SF"].push_back(sys.diffeb);
            //bias_eb["M"]=sys.nomdiffeb;
            systematics_xs_ass["Muon ID/Eff and Trigger eff SF"].push_back(sys.diffeb);
            systematics_xs_ass["Muon ID/Eff and Trigger eff SF"].push_back(sys.diffeb);

		}
		
		
		// background scaling
		else if (sys.type.find("XSDown") != string::npos || sys.type.find("XSUp") != string::npos) {
            
            if (map_scaleXS_eb.find(sys.sel) == map_scaleXS_eb.end()) {
                map_scaleXS_eb[sys.sel]=sys.diffeb;
                map_scaleXS_SF[sys.sel]=sys.diffSF;
                map_scaleXS_F[sys.sel]=sys.diffF;
                map_scaleXS_XS[sys.sel]=sys.diffxs;
            } else {
                map_scaleXS_eb[sys.sel]=max(map_scaleXS_eb[sys.sel],sys.diffeb);
                map_scaleXS_SF[sys.sel]=max(map_scaleXS_SF[sys.sel],sys.diffSF);
                map_scaleXS_F[sys.sel]=max(map_scaleXS_F[sys.sel],sys.diffF);
                map_scaleXS_XS[sys.sel]=max(map_scaleXS_XS[sys.sel],sys.diffxs);
                
            }
            
            //systematics_xs_ass["Cross Section "+sys.sel].push_back(sys.diffxs);
            
            bias_eb["Background Composition"]=sys.nomdiffeb;
            
            string s="Cross Section "+sys.sel;
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("Up") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
            else if (sys.type.find("Down") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;

            
		}
        
        // Q2 scaling
		else if (sys.type.find("ScaleDown") != string::npos || sys.type.find("ScaleUp") != string::npos) {
            
            if (map_scale_eb.find(sys.sel) == map_scale_eb.end()) {
                map_scale_eb[sys.sel]=sys.diffeb;
                map_scale_SF[sys.sel]=sys.diffSF;
                map_scale_F[sys.sel]=sys.diffF;
                map_scale_XS[sys.sel]=sys.diffxs;
            } else {
                map_scale_eb[sys.sel]=max(map_scale_eb[sys.sel],sys.diffeb);
                map_scale_SF[sys.sel]=max(map_scale_SF[sys.sel],sys.diffSF);
                map_scale_F[sys.sel]=max(map_scale_F[sys.sel],sys.diffF);
                map_scale_XS[sys.sel]=max(map_scale_XS[sys.sel],sys.diffxs);
                
            }
            
            /*if (sys.type.find("ScaleDown") != string::npos && sys.sel == "TTJets") {
                systematics_xs_ass["Factorisation Scale "+sys.sel].push_back(fabs(sys.diffxs));
            } else if (sys.type.find("ScaleUp") != string::npos && sys.sel == "TTJets") {
                systematics_xs_ass["Factorisation Scale "+sys.sel].push_back(-fabs(sys.diffxs));
            } else {
                systematics_xs_ass["Factorisation Scale "+sys.sel].push_back(sys.diffxs);
            }*/
            
            string s="Factorisation Scale "+sys.sel;
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 

            if (sys.type.find("Up") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
            else if (sys.type.find("Down") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;
            
            bias_eb["Factorisation Scale"]=sys.nomdiffeb;
            
		}
        
        // ME-PS matching
		else if (sys.type.find("MatchingDown") != string::npos || sys.type.find("MatchingUp") != string::npos) {
            
            if (map_match_eb.find(sys.sel) == map_match_eb.end()) {
                map_match_eb[sys.sel]=sys.diffeb;
                map_match_SF[sys.sel]=sys.diffSF;
                map_match_F[sys.sel]=sys.diffF;
                map_match_XS[sys.sel]=sys.diffxs;
            } else {
                map_match_eb[sys.sel]=max(map_match_eb[sys.sel],sys.diffeb);
                map_match_SF[sys.sel]=max(map_match_SF[sys.sel],sys.diffSF);
                map_match_F[sys.sel]=max(map_match_F[sys.sel],sys.diffF);
                map_match_XS[sys.sel]=max(map_match_XS[sys.sel],sys.diffxs);
                
            }
            
            //systematics_xs_ass["ME-PS Matching Threshold "+sys.sel].push_back(sys.diffxs);

            bias_eb["ME-PS Matching Threshold"]=sys.nomdiffeb;
            
            string s="ME-PS Matching Thresholde "+sys.sel;
            if (systematics_xs_ass[s].size() < 2) {
                systematics_xs_ass[s].push_back(0);
                systematics_xs_ass[s].push_back(0);
            } 
            if (sys.type.find("Up") != string::npos)
                systematics_xs_ass[s][0]=sys.diffxs;
            else if (sys.type.find("Down") != string::npos)
                systematics_xs_ass[s][1]=sys.diffxs;
		}
    }
    
    
    // finish the map for assym xs errors -- OLD
    
    /*for (std::map<std::string, vector<double> >::iterator it = systematics_xs.begin(); it != systematics_xs.end(); ++it) {
        systematics_xs_ass[it->first] = it->second;
        if (it->first == "Underlying event")
            systematics_xs_ass[it->first].push_back(-systematics_xs_ass[it->first][0]);
    }*/
    
    // BackGround composition
    
    double bgComp[4]={-1,-1,-1,-1};
    for (std::map<std::string, double>::iterator it = map_scaleXS_eb.begin(); it != map_scaleXS_eb.end(); ++it) {
        cout << "BGComp - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["Background Composition"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["Background Composition"],2)) << endl;
        cout << "BGComp - SF - " << it->first << " " << map_scaleXS_SF[it->first] << endl;
        cout << endl;
        if (bgComp[0] == -1) bgComp[0] = pow(map_scaleXS_F[it->first],2); else bgComp[0] += pow(map_scaleXS_F[it->first],2);
        if (bgComp[1] == -1) bgComp[1] = pow(map_scaleXS_eb[it->first],2); else bgComp[1] += pow(map_scaleXS_eb[it->first],2);
        if (bgComp[2] == -1) bgComp[2] = pow(map_scaleXS_SF[it->first],2); else bgComp[2] += pow(map_scaleXS_SF[it->first],2);
        if (bgComp[3] == -1) bgComp[3] = pow(map_scaleXS_XS[it->first],2); else bgComp[3] += pow(map_scaleXS_XS[it->first],2);
    }
	if (bgComp[0] != -1) systematics_f["Background Composition"].push_back(sqrt(bgComp[0]));
	if (bgComp[1] != -1) systematics_eb["Background Composition"].push_back(sqrt(bgComp[1]));
	if (bgComp[2] != -1) systematics_SF["Background Composition"].push_back(sqrt(bgComp[2]));
	if (bgComp[3] != -1) systematics_xs["Background Composition"].push_back(sqrt(bgComp[3]));

	// Q2 scaling composition
    
    double scale[4]={-1,-1,-1,-1};
    for (std::map<std::string, double>::iterator it = map_scale_eb.begin(); it != map_scale_eb.end(); ++it) {
        cout << "Q2 - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["Factorisation Scale"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["Factorisation Scale"],2)) << endl;
        cout << "Q2 - SF - " << it->first << " " << map_scale_SF[it->first] << endl;
        cout << endl;
        //cout << it->first << " " << it->second << endl;
        if (scale[0] == -1) scale[0] = pow(map_scale_F[it->first],2); else scale[0] += pow(map_scale_F[it->first],2);
        if (scale[1] == -1) scale[1] = pow(map_scale_eb[it->first],2); else scale[1] += pow(map_scale_eb[it->first],2);
        if (scale[2] == -1) scale[2] = pow(map_scale_SF[it->first],2); else scale[2] += pow(map_scale_SF[it->first],2);
        if (scale[3] == -1) scale[3] = pow(map_scale_XS[it->first],2); else scale[3] += pow(map_scale_XS[it->first],2);
    }
	if (scale[0] != -1) systematics_f["Factorisation Scale"].push_back(sqrt(scale[0]));
	if (scale[1] != -1) systematics_eb["Factorisation Scale"].push_back(sqrt(scale[1]));
	if (scale[2] != -1) systematics_SF["Factorisation Scale"].push_back(sqrt(scale[2]));
	if (scale[3] != -1) systematics_xs["Factorisation Scale"].push_back(sqrt(scale[3]));

    // ME-PS matching composition
    
    double match[4]={-1,-1,-1,-1};
    for (std::map<std::string, double>::iterator it = map_match_eb.begin(); it != map_match_eb.end(); ++it) {
        cout << "MEPS - EB - " << it->first << " " << it->second << "^2 * " << bias_eb["ME-PS Matching Threshold"] << "^2 = " << sqrt(pow(it->second,2)+pow(bias_eb["ME-PS Matching Threshold"],2)) << endl;
        cout << "MEPS - SF - " << it->first << " " << map_match_SF[it->first] << endl;
        cout << endl;
        //cout << it->first << " " << it->second << endl;
        if (match[0] == -1) match[0] = pow(map_match_F[it->first],2); else match[0] += pow(map_match_F[it->first],2);
        if (match[1] == -1) match[1] = pow(map_match_eb[it->first],2); else match[1] += pow(map_match_eb[it->first],2);
        if (match[2] == -1) match[2] = pow(map_match_SF[it->first],2); else match[2] += pow(map_match_SF[it->first],2);
        if (match[3] == -1) match[3] = pow(map_match_XS[it->first],2); else match[3] += pow(map_match_XS[it->first],2);
    }
	if (match[0] != -1) systematics_f["ME-PS Matching Threshold"].push_back(sqrt(match[0]));
	if (match[1] != -1) systematics_eb["ME-PS Matching Threshold"].push_back(sqrt(match[1]));
	if (match[2] != -1) systematics_SF["ME-PS Matching Threshold"].push_back(sqrt(match[2]));
	if (match[3] != -1) systematics_xs["ME-PS Matching Threshold"].push_back(sqrt(match[3]));


 	/*cout << "Pdf uncertainty Eb: " << max(sqrt(sumplusE),sqrt(summinE)) << endl;
 	cout << "Pdf uncertainty XS: " << max(sqrt(sumplusXS),sqrt(summinXS)) << endl;
	*/
	if (nPdf > 0) systematics_f["PDF Uncertainties"].push_back(max(sqrt(sumplusF),sqrt(summinF)));
	if (nPdf > 0) systematics_eb["PDF Uncertainties"].push_back(max(sqrt(sumplusE),sqrt(summinE)));
	if (nPdf > 0) systematics_xs["PDF Uncertainties"].push_back(max(sqrt(sumplusXS),sqrt(summinXS)));
	if (nPdf > 0) systematics_SF["PDF Uncertainties"].push_back(max(sqrt(sumplusSF),sqrt(summinSF)));
    
    systematics_xs_ass["PDF Uncertainties"].push_back(sqrt(sumplusXS));
    systematics_xs_ass["PDF Uncertainties"].push_back(sqrt(summinXS));
    
    //systematics_xs_ass["Luminosity"].push_back(7.2);
    //systematics_xs_ass["Luminosity"].push_back(-7.2);

    //systematics_xs_ass["PDF"].push_back(2.3);
    //systematics_xs_ass["PDF"].push_back(2.7);

    //systematics_xs["Muon ID/Eff and Trigger eff SF"].push_back(2.3);
    //systematics_xs_ass["Muon ID/Eff and Trigger eff SF"].push_back(2.3);
    //systematics_xs_ass["Muon ID/Eff and Trigger eff SF"].push_back(-2.3);
    
    //systematics_xs["PDF Uncertainties"].push_back(2.1);
    //systematics_xs_ass["PDF Uncertainties"].push_back(1.9);
    //systematics_xs_ass["PDF Uncertainties"].push_back(2.1);
    
    // use 2011 values!!!
    
    /*if (chan == "Mu") {
        //systematics_xs["UE"].push_back(nomxs*(1.7/161.4));
        //systematics_xs_ass["UE"].push_back(nomxs*(1.7/161.4));
        //systematics_xs_ass["UE"].push_back(-nomxs*(1.7/161.4));
        
        systematics_xs_ass["PDF Uncertainties"].clear();
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.015));
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.019));
        
        systematics_xs_ass["ME-PS Matching Threshold"].clear();
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.022));
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.059));
        
        systematics_xs_ass["Factorisation Scale"].clear();
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.061));
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.038));
        
        systematics_xs_ass["Top Quark Mass"].clear();
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.000));
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.011));
        
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        
        //systematics_SF["Top Quark Mass"].clear();
        //systematics_SF["Top Quark Mass"].push_back(0.011);
        
        //systematics_SF["Method bias"].clear();
        //systematics_SF["Method bias"].push_back(4.6);
        
        //systematics_xs_ass["Background Composition"].clear();
        //systematics_xs_ass["Background Composition"].push_back(nomxs*(1.9/161.4));
        //systematics_xs_ass["Background Composition"].push_back(nomxs*(1.8/161.4));
        
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(1.4/225.2));
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(1.4/225.2));
        
        systematics_xs_ass["QCD Normalization"].push_back((1/225.2)*nomxs);
        systematics_xs_ass["QCD Normalization"].push_back((1/225.2)*nomxs);
        
    } 
    
    else if (chan == "El") {
        
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        
        systematics_xs_ass["Top Quark Mass"].clear();
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.007));
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.019));
        
        systematics_xs_ass["PDF Uncertainties"].clear();
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.018));
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.021));
        
        systematics_xs_ass["ME-PS Matching Threshold"].clear();
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.054));
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.041));
        
        systematics_xs_ass["Factorisation Scale"].clear();
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.087));
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.025));
        
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(0.9/223.2));
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(0.9/223.2));
        
        systematics_xs_ass["QCD Normalization"].push_back((2/223.2)*nomxs);
        systematics_xs_ass["QCD Normalization"].push_back((2/223.2)*nomxs);

        
    } 
    else if (chan == "Comb") {
        
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        systematics_xs_ass["btag"].push_back(nomxs*(0.075));
        
        systematics_xs_ass["Top Quark Mass"].clear();
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.003));
        systematics_xs_ass["Top Quark Mass"].push_back(nomxs*(0.014));
        
        systematics_xs_ass["PDF Uncertainties"].clear();
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.016));
        systematics_xs_ass["PDF Uncertainties"].push_back(nomxs*(0.020));
    
        systematics_xs_ass["ME-PS Matching Threshold"].clear();
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.046));
        systematics_xs_ass["ME-PS Matching Threshold"].push_back(nomxs*(0.031));
        
        systematics_xs_ass["Factorisation Scale"].clear();
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.062));
        systematics_xs_ass["Factorisation Scale"].push_back(nomxs*(0.021));
        
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(2.1/224.1));
        systematics_xs_ass["W+jets template from 7TeV"].push_back(nomxs*(2.1/224.1));
        
        systematics_xs_ass["QCD Normalization"].push_back((2/223.2)*nomxs);
        systematics_xs_ass["QCD Normalization"].push_back((2/223.2)*nomxs);

        
    }*/
    //systematics_xs["Luminosity"].push_back(159.5*0.04);
    //systematics_xs_ass["Luminosity"].push_back(159.5*0.04);
    //systematics_xs_ass["Luminosity"].push_back(-159.5*0.04);
    
    // merge theory systematics again
    
    std::map<std::string, std::vector<double> > tmp_systematics_xs_ass;
    std::map<std::string, std::vector<double> > systematics_xs_ass_perc;
    
    double xs_sum_p = 0;
    double xs_sum_m = 0;
    
    double q2_sum_p = 0;
    double q2_sum_m = 0;
    
    double match_sum_p = 0;
    double match_sum_m = 0;
    
    for (std::map<std::string, std::vector<double> >::const_iterator it=systematics_xs_ass.begin(); it != systematics_xs_ass.end(); ++it)
        tmp_systematics_xs_ass[it->first]=it->second;

    systematics_xs_ass.clear();

    for (std::map<std::string, std::vector<double> >::const_iterator it=tmp_systematics_xs_ass.begin(); it != tmp_systematics_xs_ass.end(); ++it) {
        
        if (it->first.find("Cross Section") != string::npos) {
            xs_sum_p = xs_sum_p + pow(it->second[0],2);
            xs_sum_m = xs_sum_m + pow(it->second[1],2);
        }
        else if (it->first.find("Factorisation") != string::npos) {
            q2_sum_p = q2_sum_p + pow(it->second[0],2);
            q2_sum_m = q2_sum_m + pow(it->second[1],2);
        }
        else if (it->first.find("ME-PS") != string::npos) {
            match_sum_p = match_sum_p + pow(it->second[0],2);
            match_sum_m = match_sum_m + pow(it->second[1],2);
        }
        
        else
            systematics_xs_ass[it->first]=it->second;

    
    }
    
    systematics_xs_ass["Background Composition"].push_back(sqrt(xs_sum_p));
    systematics_xs_ass["Background Composition"].push_back(sqrt(xs_sum_m));

    systematics_xs_ass["Factorisation Scale"].push_back(sqrt(q2_sum_p));
    systematics_xs_ass["Factorisation Scale"].push_back(sqrt(q2_sum_m));
    
    systematics_xs_ass["ME-PS Matching Threshold"].push_back(sqrt(match_sum_p));
    systematics_xs_ass["ME-PS Matching Threshold"].push_back(sqrt(match_sum_m));
    
    for (std::map<std::string, std::vector<double> >::const_iterator it=systematics_xs_ass.begin(); it != systematics_xs_ass.end(); ++it) {
        systematics_xs_ass_perc[it->first].push_back((it->second[0]/nomxs)*100);
        systematics_xs_ass_perc[it->first].push_back((it->second[1]/nomxs)*100);
    }
    
    /*systematics_xs_ass_perc["Lepton ID/eff"].push_back(1.0);
    systematics_xs_ass_perc["Lepton ID/eff"].push_back(-1.0);
    systematics_xs_ass_perc["Lepton Trigger eff"].push_back(1.0);
    systematics_xs_ass_perc["Lepton Trigger eff"].push_back(-1.0);

    
    systematics_xs_ass["Lepton ID/eff"].push_back(0.01*nomxs);
    systematics_xs_ass["Lepton ID/eff"].push_back(-0.01*nomxs);
    systematics_xs_ass["Lepton Trigger eff"].push_back(0.01*nomxs);
    systematics_xs_ass["Lepton Trigger eff"].push_back(-0.01*nomxs);*/

    cout << nomxs << endl;

	std::map<std::string,double> empty;
    
	printSystematics("F_CS",systematics_f,empty,100);	
	printSystematics("E_b",systematics_eb,empty,10);	
	printSystematics("SF",systematics_SF,empty,10);	
	printSystematics("XS",systematics_xs,empty,10);	
	printAssymSystematics("XS",systematics_xs_ass,empty,10);	
	printAssymSystematics("XS",systematics_xs_ass_perc,empty,10);	
	
    //cout << "PU syst " << systematics_xs_ass["PileUp"][0] << endl;
    //cout << "PU syst " << systematics_xs_ass["PileUp"][1] << endl;
	cout << endl;
    
    double sup = 0;
    double sdown = 0;
    for (int i=0;i<pdf_up.size(); i++) {
        
        sup+=pow(maxPDF(pdf_up[i],pdf_down[i]),2);

        sdown+=pow(maxPDF(-pdf_up[i],-pdf_down[i]),2);

    }
    
    //cout << sqrt(sup) << " " << sqrt(sdown) << endl;
    	
	return 0;
	
}
