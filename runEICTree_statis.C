#include <iostream>
#include <vector>
#include <sstream>
#include <string>


#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>

#include "Riostream.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

#include "TLorentzVector.h"
#include "TBranchElement.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <TH1I.h>

#include "TSpline.h"
//#include "acceptance.h"
//#include "RootUtils.C"
#define PI            3.1415926

#define MASS_MUON     0.1056
#define MASS_ELECTRON 0.000511
#define MASS_JPSI     3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8756129
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208
#define ENERGY_TARGET 490.0
#define NUMBER_EVNT   100000


//#define MASS_NUCLEON  0.93903412 //average (N*MASS_NEUTRON+Z*MASS_PROTON)/Z+N

using namespace std;
using namespace erhic;

//map<string, TH1I*> hist_distribution_tracks;

//const vector<string> GreyTracks = {"1","2","3","4","5","6","7","8"};

//double Q2_data[82]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

void runEICTree_statis(string filename="", const int nEvents =100000000){

    double meansss;
    double scale_charged_target        ;
    double scale_charged_central       ;
    double scale_charged_projectile    ;
    double scale_positive_target       ;
    double scale_positive_central      ;
    double scale_positive_projectile   ;
    double scale_negative_target       ;
    double scale_negative_central      ;
    double scale_negative_projectile   ;
    double scale_charged_forward       ;
    double scale_charged_backward      ;
    double scale_charged               ;
    double scale_positive              ;
    double scale_negative              ;
    double scale_charged_nong          ;
    double scale_positive_nong         ;
    double scale_negative_nong         ;

    std::ifstream fileIN,fileINmeans,fileINmeans_D; 
    
    string filename_D = "scale_g0_qhat_05_1E8_TEST"; //change every time - text file with the scale

    string filename_1 ;

    fileINmeans.open(("muD-means/"+filename_D+".txt").c_str());


        fileINmeans >>scale_charged_target     
                    >>scale_charged_central    
                    >>scale_charged_projectile 
                    >>scale_positive_target    
                    >>scale_positive_central   
                    >>scale_positive_projectile
                    >>scale_negative_target    
                    >>scale_negative_central   
                    >>scale_negative_projectile
                    >>scale_charged       
                    >>scale_positive
                    >>scale_negative
                    >>scale_charged_nong
                    >>scale_positive_nong
                    >>scale_negative_nong;      

                    cout << "scale_charged_target      = " <<scale_charged_target      << endl;
                    cout << "scale_charged_central     = " <<scale_charged_central     << endl;
                    cout << "scale_charged_projectile  = " <<scale_charged_projectile  << endl;
                    cout << "scale_positive_target     = " <<scale_positive_target     << endl;
                    cout << "scale_positive_central    = " <<scale_positive_central    << endl;
                    cout << "scale_positive_projectile = " <<scale_positive_projectile << endl;
                    cout << "scale_negative_target     = " <<scale_negative_target     << endl;
                    cout << "scale_negative_central    = " <<scale_negative_central    << endl;
                    cout << "scale_negative_projectile = " <<scale_negative_projectile << endl;
                    cout << "scale_charged             = " <<scale_charged             << endl;
                    cout << "scale_positive            = " <<scale_positive            << endl;
                    cout << "scale_negative            = " <<scale_negative            << endl;
                    cout << "scale_charged_nong        = " <<scale_charged_nong        << endl;
                    cout << "scale_positive_nong       = " <<scale_positive_nong       << endl;
                    cout << "scale_negative_nong       = " <<scale_negative_nong       << endl;
        fileINmeans_D.close();
		fileINmeans.close();
    
	TString basepath = "/gpfs/mnt/gpfs02/eic/crob/BeAGLE-dev-2021-04-27/BeAGLE/output/e665/qhat_05/sim-ene-22/muXe/g0/"; //path root files - change at each option (/g1 or /g2 or /g3)	
	TString basename = "Output_input_temp";	//basename of the rooot files
    std::cout << "Creating TChain for " << filename << std::endl;
    TChain *tree = new TChain("EICTree");
    //loop over the root files using TChain - mixing the trees in one.
	//idx = number of root files
     for (int idx=1; idx<1000+1; ++idx) {
         TString path = Form("%s/%s_%d.root",basepath.Data(),basename.Data(),idx);
		 cout << "Adding path = " << path << endl;
         tree->Add(path);
     }


    std::cout << "All TTrees added." << std::endl;
    EventBeagle* event(NULL);
    tree->SetBranchAddress("event", &event);
    
    TFile *out = new TFile ("analysis/Analysis_ratio_combined_Output_muXe_g0_qhat_05_1E8_updated.root", "RECREATE"); //root output file change at each option (/g1 or /g2 or /g3)
    //####################  Creating Histograms ####################################
    double ng_bmin[11] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
    double zng_bmax[11] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0};
    double bin_ng[11] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0};
    
    TH1F *h_nTracks     = new TH1F("h_nTracks", "h_nTracks", 100, 90, 300);
    TH1F *h_pdg         = new TH1F("h_pdg", "h_pdg", 500,-2500, 2500);
    TH1F *h_Energy      = new TH1F("h_Energy", "h_Energy", 100, 0,6);
    TH1F *h_dE_dN       = new TH1F("h_idE_dN", "h_dE_dN", 100, 0, 200);
    TH1F *h_Status      = new TH1F("h_Status", "h_Status", 50, 0, 25);
    TH1F *h_pt          = new TH1F("h_pt","h_pt", 100,0, 2);
    TH1F *h_pz          = new TH1F("h_pz","h_pz", 80,-40, 40);
    TH1F *h_eta         = new TH1F("h_eta", "h_eta", 100, -10, 10);
    
    TH1F *h_rapidity            = new TH1F("h_rapidity","h_rapidity", 100, -8, 5);
    TH1F *h_b                   = new TH1F("h_b", "h_b", 100, 0, 10);
    TH1F *h_trueQ2              = new TH1F("h_trueQ2", "h_trueQ2", 100, 0, 100);
    TH1F *h_trueQ2_weighted     = new TH1F("h_trueQ2_weighted", "h_trueQ2_weighted", 100, 0, 100);
    TH1F *h_Energy_proton       = new TH1F("h_Energy_proton", "h_Energy_proton", 100, 0, 3);
    TH1F *h_Energy_proton_gt    = new TH1F("h_Energy_proton_gt", "h_Energy_proton_gt", 100, 0, 3);
    TH1F *h_xf                  = new TH1F("h_xf", "h_xf", 100, -1.1, 1.1);
	TH1F *h_Energy_parton 		= new TH1F("h_Energy_parton","h_Energy_parton",100,0,6);    
    TH1F *h_w2          = new TH1F("h_W2","h_W2", 100, 0, 1000);
    TH1F *h_nu          = new TH1F("h_nu", "h_nu",100,50,400);
    TH1F *h_mom         = new TH1F("h_mom", "h_mom", 60, -30, 30);
    TH1F *h_p_cms       = new TH1F("h_p_cms", "h_p_cms", 60,0,30);
    TH1F *h_Energy_cms  = new TH1F("h_Energy_cms","h_Energy_cms", 100,0,6);
    TH1F *h_pz_cms      = new TH1F("h_pz_cms","h_pz_cms",80, -40,40);
    TH1F *h_pt_cms      = new TH1F("h_pt_cms","h_pt_cms",100, 0,2);
    TH1F *h_rap_cms     = new TH1F("h_rap_cms","h_rap_cms",48,-4,4);
    TH1F *h_xf_cms      = new TH1F("h_xf_cms","h_xf_cms",100,-1.1,1.1);
    TH1F *h_x_bj        = new TH1F("h_x_bj","h_x_bj",100,0,1);
    TH1F *h_theta       = new TH1F("h_theta","h_theta",100,0,360);
    TH1F *h_rap_lab             = new TH1F("h_rap_lab","h_rap_lab", 48,-4,4);
    TH1F *h_rap_lab_computed    = new TH1F("h_rap_lab_computed","h_rap_lab_computed", 48,-4,4);
    TH1F *h_mom_charged = new TH1F("h_mom_charged","h_mom_charged",100,-5,30);
    TH1F *h_1_over_x_bj     = new TH1F("h_1_over_x_bj","h_1_over_x_bj",100,-0.5,1);

	TH1F *h_pt_proton       = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
    TH1F *h_pt_proton_gt    = new TH1F("h_pt_proton_gt", "h_pt_proton_gt", 100, 0, 2);
    TH1F *h_nTracks_gt      = new TH1F("h_nTracks_gt","h_nTracks_gt",12, 0, 12);
    TH1F *h_pt_gt           = new TH1F("h_pt_gt", "h_pt_gt", 100, 0, 1);
    TH1F *h_Energy_gt       = new TH1F("h_Energy_gt", "h_Energy_gt", 100,0, 6);
    TH1F *h_rapidity_gt     = new TH1F("h_rapidity_gt","h_rapidity_gt",100,-10,10);
    TH1F *h_nTracks_p_gt    = new TH1F("h_nTracks_p_gt","h_nTracks_p_gt", 11, -0.5, 10.5); //****
    TH1F *h_nTracks_proton  = new TH1F("h_nTracks_proton", "h_nTracks_proton", 11, -0.5, 10.5); //****
    TH1F *h_pz_proton       = new TH1F("h_pz_proton","h_pz_proton", 40,-10,10);
    TH1F *h_pz_proton_cms   = new TH1F("h_pz__proton_cms","h_pz__proton_cms", 40,-10,10);

    TH1F *h_negative_Tracks = new TH1F("h_negative_Tracks","h_negative_Tracks",24,0,12);
    TH1F *h_positive_Tracks = new TH1F("h_positive_Tracks","h_positive_Tracks",24,0,12);
    TH1F *h_charged_Tracks  = new TH1F("h_charged_Tracks","h_charged_Tracks",24,0,12);
    TH1F *h_negative_nong_Tracks = new TH1F("h_negative_nong_Tracks","h_negative_nong_Tracks",24,0,12);
    TH1F *h_positive_nong_Tracks = new TH1F("h_positive_nong_Tracks","h_positive_nong_Tracks",24,0,12);
    TH1F *h_charged_nong_Tracks  = new TH1F("h_charged_nong_Tracks","h_charged_nong_Tracks",24,0,12); 



    TH1F *h_negative_nTracks_target     = new TH1F("h_negative_nTracks_target","Histogram title ; X-axis title ; Y-axis title", 24, 0, 12);
    TH1F *h_negative_nTracks_central    = new TH1F("h_negative_nTracks_central","h_negative_nTracks_central", 24, 0, 12);
    TH1F *h_negative_nTracks_projectile = new TH1F("h_negative_nTracks_projectile","h_negative_nTracks_projectile", 24, 0, 12);
    TH1F *h_positive_nTracks_target     = new TH1F("h_positive_nTracks_target","h_positive_nTracks_target", 24, 0, 12);
    TH1F *h_positive_nTracks_central    = new TH1F("h_positive_nTracks_central","h_positive_nTracks_central", 24, 0, 12);
    TH1F *h_positive_nTracks_projectile = new TH1F("h_positive_nTracks_projectile","h_positive_nTracks_projectile", 24, 0, 12);
    TH1F *h_charged_nTracks_target      = new TH1F("h_charged_nTracks_target","h_charged_nTracks_target", 24, 0, 12);
    TH1F *h_charged_nTracks_central     = new TH1F("h_charged_nTracks_central","h_charged_nTracks_central", 24, 0, 12);
    TH1F *h_charged_nTracks_projectile  = new TH1F("h_charged_nTracks_projectile","h_charged_nTracks_projectile", 24, 0, 12);

    


    double bins[14]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5};
   
    TH1F *h_nTracks_pos_pions                   = new TH1F("h_nTracks_pos_pions","h_nTracks_pos_pions",21,-0.5, 20.5);
    TH2F *h_nTracks_Grey_tracks                 = new TH2F("h_nTracks_Grey_tracks","h_nTracks_Grey_tracks", 11,0,11,60,0,15);
    TProfile *prof_nTracks_Grey_tracks          = new TProfile("prof_nTracks_Grey_tracks","prof_nTracks_Grey_tracks", 13,bins); 
    TProfile *prof_nTracks_target_positive      = new TProfile("prof_nTracks_target_positive","prof_nTracks_target_positive", 13,bins);
    TProfile *prof_nTracks_central_positive     = new TProfile("prof_nTracks_central_positive","prof_nTracks_central_positive", 13,bins);
    TProfile *prof_nTracks_projectile_positive  = new TProfile("prof_nTracks_projectile_positive","prof_nTracks_projectile_positive", 13,bins);
    TProfile *prof_nTracks_target_negative      = new TProfile("prof_nTracks_target_negative","prof_nTracks_target_negative", 13,bins);
    TProfile *prof_nTracks_central_negative     = new TProfile("prof_nTracks_central_negative","prof_nTracks_central_negative", 13,bins);
    TProfile *prof_nTracks_projectile_negative  = new TProfile("prof_nTracks_projectile_negative","prof_nTracks_projectile_negative", 13,bins);
    TProfile *prof_nTracks_target_charged       = new TProfile("prof_nTracks_target_charged","prof_nTracks_target_charged", 13,bins);
    TProfile *prof_nTracks_central_charged      = new TProfile("prof_nTracks_central_charged","prof_nTracks_central_charged", 13,bins);
    TProfile *prof_nTracks_projectile_charged   = new TProfile("prof_nTracks_projectile_charged","prof_nTracks_projectile_charged", 13,bins);
    TProfile *prof_tracks_backward_charged      = new TProfile("prof_tracks_backward_charged","prof_tracks_backward_charged",13,bins);
    TProfile *prof_tracks_forward_charged      = new TProfile("prof_tracks_forward_charged","prof_tracks_forward_charged",13,bins);

    TProfile *prof_nTracks_target_positive_contaminated      = new TProfile("prof_nTracks_target_positive_contaminated","prof_nTracks_target_positive_contaminated", 13,bins);
    TProfile *prof_nTracks_central_positive_contaminated     = new TProfile("prof_nTracks_central_positive_contaminated","prof_nTracks_central_positive_contaminated", 13,bins);
    TProfile *prof_nTracks_projectile_positive_contaminated  = new TProfile("prof_nTracks_projectile_positive_contaminated","prof_nTracks_projectile_positive_contaminated", 13,bins);
    TProfile *prof_nTracks_target_negative_contaminated      = new TProfile("prof_nTracks_target_negative_contaminated","prof_nTracks_target_negative_contaminated", 13,bins);
    TProfile *prof_nTracks_central_negative_contaminated     = new TProfile("prof_nTracks_central_negative_contaminated","prof_nTracks_central_negative_contaminated", 13,bins);
    TProfile *prof_nTracks_projectile_negative_contaminated  = new TProfile("prof_nTracks_projectile_negative_contaminated","prof_nTracks_projectile_negative_contaminated", 13,bins);
    TProfile *prof_nTracks_target_charged_contaminated       = new TProfile("prof_nTracks_target_charged_contaminated","prof_nTracks_target_charged_contaminated", 13,bins);
    TProfile *prof_nTracks_central_charged_contaminated      = new TProfile("prof_nTracks_central_charged_contaminated","prof_nTracks_central_charged_contaminated", 13,bins);
    TProfile *prof_nTracks_projectile_charged_contaminated   = new TProfile("prof_nTracks_projectile_charged_contaminated","prof_nTracks_projectile_charged_contaminated", 13,bins);

    TProfile *prof_nTracks_target_positive_weighted      = new TProfile("prof_nTracks_target_positive_weighted","prof_nTracks_target_positive_weighted", 13,bins);
    TProfile *prof_nTracks_central_positive_weighted     = new TProfile("prof_nTracks_central_positive_weighted","prof_nTracks_central_positive_weighted", 13,bins);
    TProfile *prof_nTracks_projectile_positive_weighted  = new TProfile("prof_nTracks_projectile_positive_weighted","prof_nTracks_projectile_positive_weighted", 13,bins);
    TProfile *prof_nTracks_target_negative_weighted      = new TProfile("prof_nTracks_target_negative_weighted","prof_nTracks_target_negative_weighted", 13,bins);
    TProfile *prof_nTracks_central_negative_weighted     = new TProfile("prof_nTracks_central_negative_weighted","prof_nTracks_central_negative_weighted", 13,bins);
    TProfile *prof_nTracks_projectile_negative_weighted  = new TProfile("prof_nTracks_projectile_negative_weighted","prof_nTracks_projectile_negative_weighted", 13,bins);
    TProfile *prof_nTracks_target_charged_weighted       = new TProfile("prof_nTracks_target_charged_weighted","prof_nTracks_target_charged_weighted", 13,bins);
    TProfile *prof_nTracks_central_charged_weighted      = new TProfile("prof_nTracks_central_charged_weighted","prof_nTracks_central_charged_weighted", 13,bins);
    TProfile *prof_nTracks_projectile_charged_weighted   = new TProfile("prof_nTracks_projectile_charged_weighted","prof_nTracks_projectile_charged_weighted", 13,bins);

    TProfile *prof_charged_tracks_rapidity      = new TProfile("prof_charged_tracks_rapidity","prof_charged_tracks_rapidity",20,-5,5);
    TProfile *prof_positive_tracks_rapidity     = new TProfile("prof_positive_tracks_rapidity","prof_positive_tracks_rapidity",20,-5,5);
    TProfile *prof_negative_tracks_rapidity     = new TProfile("prof_negative_tracks_rapidity","prof_negative_tracks_rapidity",20,-5,5);
	TProfile *prof_charged_nong_rapidity      = new TProfile("prof_charged_nong_rapidity","prof_charged_nong_rapidity",20,-5,5);
    TProfile *prof_positive_nong_rapidity      = new TProfile("prof_positive_nong_rapidity","prof_positive_nong_rapidity",20,-5,5);
	TProfile *prof_negative_nong_rapidity      = new TProfile("prof_negative_nong_rapidity","prof_negative_nong_rapidity",20,-5,5);
	
	TProfile *prof_Q_hadronic_charge_gt         = new TProfile("prof_Q_hadronic_charge_gt","prof_Q_hadronic_charge_gt", 10,bins);
    //TProfile *prof_Q_hadronic_charge_gt_backward =new TProfile("prof_Q_hadronic_charge_gt_backward","prof_Q_hadronic_charge_gt_backward", 10,bins);
    TProfile *prof_Q_hadronic_charge_backward_gt    = new TProfile("prof_Q_hadronic_charge_backward_gt","prof_Q_hadronic_charge_backward_gt",13,bins);
    TProfile *prof_Q_hadronic_charge_forward_gt     = new TProfile("prof_Q_hadronic_charge_forward_gt","prof_Q_hadronic_charge_forward_gt",13,bins);
    TProfile *prof_distance_grey_tracks             = new TProfile("prof_distance_grey_tracks","prof_distance_grey_tracks",13,bins);
    TProfile *prof_distance_grey_tracks_pions       = new TProfile("prof_distance_grey_tracks_pions","prof_distance_grey_tracks_pions",13,bins);
	TProfile *prof_Energygt_gt                      = new TProfile("prof_Energygt_gt","prof_Energygt_gt",13,bins);
	TH1F *h_average_distance                        = new TH1F("h_average_distance","h_average_distance",13,bins);
    TH1F *h_impact_parameter                        = new TH1F("h_impact_parameter","h_impact_parameter",13,0,13);
    TH1F *h_thickness                               = new TH1F("h_thickness","h_thickness",13,0,13);

    
      
    TProfile *prof  = new TProfile("prof", "Rapidity;<ng>", 21, -7, 3);
    TH2F *hist      = new TH2F("hist","Rapidity;<ng>", 21, -4, 4,11,-0.5,10.5);
    TH2F *h_dist_ng = new TH2F("h_dist_ng","h_dist_ng",13,-0.5,12.5,13,0,13);
    TH2F *h_dist_ng_pions = new TH2F("h_dist_ng_pions","h_dist_ng_pions",13,-0.5,12.5,13,0,13);
	TH2F *h_Energygt_gt = new TH2F("h_Energygt_gt", "h_Energygt_gt",13,-0.5,12.5, 30, 0, 3);

    int counter_total_particles = 0;
    int counter_event           = 0;
	int counter_total_ng		= 0;
	int counter_total_contamination = 0;
	int counter_total_pions = 0;
	int counter_total_kaons = 0;
    string title_collision = "e - Xe";

    TLorentzVector mu_scattered;
    TLorentzVector e_beam;
    TLorentzVector p_gamma;
    TLorentzVector p_final;
    TLorentzVector p_beam;
    TLorentzVector p_gamma_beagle;

    TLorentzVector p_final_rotated_x;
    TLorentzVector p_final_rotated_z;
    TLorentzVector p_final_boosted_z; 

    TVector3 ux;
    TVector3 uz;
    TVector3 uz_beagle;
    TVector3 ux_beagle;

    TVector3 p_final_rotated;
    TVector3 p_gamma_rotated;
    TVector3 p_gamma_rotated_2;

    TMatrixD rot_x(3,3);
    TMatrixD rot_z(3,3);
    TMatrixD boost_z(4,4);
    TVectorD p_gamma_boosted(4);
    TVectorD p_gamma_boosted_(4);
    
    TMatrixD P_1;
    double phi_rot1 = 0;
    double theta_rot2=0;
    
    double rap_lab              = 0;
    double rap_cms              = 0;
    double rap_lab_computed     = 0;
    double xf_cms               = 0;
    double W2_boosted           = 0;
    double rap_e_lab            = 0;
    double Beta                 = 0;
    int counter_final_e         = 0;
    int p_gamma_counter         = 0;
    int trys =0, trys1=0;
    int protons_lab = 0, neutrons_lab = 0;
    int struck_neutron          = 0;
    int struck_proton           = 0;

    double Q_hadronic_charge         = 0;
    double Q_hadronic_charge_backward= 0;
    double Q_hadronic_charge_forward = 0;

    double sum_N_prot_gt_0 = 0;
    double average_N_prot  = 0;
    
    
    int S_boosted               = 0;
    double S                    = 0;
    double gamma                = 0;
    int cms_protons_xf_cut      = 0;

    double acc,ww =0;

    //auto graph_n_vs_ng = TH1TOTGraph(prof_nTracks_Grey_tracks);

    double acc_x[82]={1.487,1.561,1.639,1.720,1.794,1.870,1.936,2.005,2.062,2.120,2.195,2.257,2.320,2.386,2.453,2.540,2.685,2.819,2.960,3.086,3.218,3.355,3.450,3.572,3.698,
                            3.829,3.937,4.048,4.191,4.310,4.462,4.620,4.784,4.885,5.022,5.164,5.347,5.575,5.812,6.060,6.362,6.634,6.916,7.262,7.624,8.060,8.462,8.885,9.328,9.793,
                            10.354,11.100,11.735,12.407,13.117,14.062,15.075,16.161,17.086,18.316,19.500,20.760,21.948,23.204,24.532,25.577,27.804,30.225,32.856,35.223,37.761,
                            39.370,42.798,47.836,45.247,50.574,53.098,55.747,58.938,63.184,67.735,74.146 }; 
    double acc_y[82]={0.386,0.391,0.396,0.400,0.405,0.410,0.415,0.420,0.425,0.432,0.439,0.447,0.454,0.461,0.468,0.478,0.490,0.505,0.519,0.534,0.546,0.561,0.575,0.587,0.604,0.619,
                            0.633,0.650,0.665,0.682,0.697,0.711,0.726,0.740,0.757,0.774,0.789,0.801,0.811,0.823,0.833,0.837,0.842,0.847,0.852,0.854,0.859,0.862,0.867,0.867,0.867,0.869,
                            0.871,0.871,0.871,0.871,0.876,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.881,0.876,0.874,0.874,
                            0.874,0.874,0.874,0.874};

    TGraph *tg_acc_2 = new TGraph(82,acc_x,acc_y);
    TF1* fcn = new TF1("fcn", "[0] / (1 + TMath::Exp(-[1] * (x - [2])))", 1, 100);
    TGraph *tg_acc_1 = new TGraph(82,acc_x,acc_y);
    TSpline5 *s = new TSpline5("tg3_s",tg_acc_1);
    tg_acc_2->Fit("fcn");

    double pzlep = 0;
    double pztarg = 0;
    int struck_nucleon = 0;
    double MASS_NUCLEON;
    double trueQ2=0;       
    double trueW2=0;       
    double trueX =0;       
    double trueY =0;        
    double trueNu =0;
    double s_hat =0;       
    double t_hat =0;   
    double u_hat =0;        
    double photon_flux =0;  
    int event_process =0 ;  
    int nParticles =0;
    double Tb = 0;
    double impact_parameter=0;
    double distance =0;         
    int N_nevap     =0;        
    int N_pevap     =0; 
    double x_bj     =0;

	int charged_total = 0, positive_total = 0, negative_total =0;
	int negative_hadrons_total = 0;
    //----------------------- particle loop initializations -------------------------------
    const erhic::ParticleMC* particle;

    int pdg         = 0;
    int status      = 0;
    int index       = 0;
    double pt       = 0;
    double eta      = 0;
    double phi      = 0;
    double rap      = 0;
    double mass     = 0;
    double theta    = 0;
    theta = theta*1000.0;
    double p        = 0;
    double Energy   = 0;
    double xf       = 0;
    double px       = 0;
    double py       = 0;
    double pz       = 0;
    int parent_index = 0;
    double parent_id = 0;
            
           
    int counter_pos_pions   = 0;
    int counter_neg_pions   = 0;
    double pt_cms           = 0;
    double pt_lab           = 0;
    double xf_2             = 0;
    double p_final_pz       = 0;  

    trys = 0;
    trys1= 0;    
    int counter_p_gt = 0;
	int count_process_99 = 0;
	int count_process_131_132 = 0;
	int count_process_135_136 = 0;
	int cuts_count_process_99 = 0;
	int cuts_count_process_131_132 = 0;
	int cuts_count_process_135_136 = 0;
    int count_process_91_94=0;
	int count_process_11_13=0;
	int count_process_95=0;	
 
	int count_process_99_shad = 0;
    int count_process_131_132_shad = 0;
    int count_process_135_136_shad = 0;
	int count_process_91_94_shad=0;
    int count_process_11_13_shad=0;
    int count_process_95_shad=0;
    int count_process_99_non_shad = 0;
    int count_process_131_132_non_shad = 0;
    int count_process_135_136_non_shad = 0;
    int count_process_91_94_non_shad=0;
    int count_process_11_13_non_shad=0;
    int count_process_95_non_shad=0;
	int count_process_99_between = 0;
    int count_process_131_132_between = 0;
    int count_process_135_136_between = 0;
	int count_process_91_94_between=0;
    int count_process_11_13_between=0;
    int count_process_95_between=0;
	int counter_event_cuts = 0;
	double Energy_proton_gt     = 0;
	double Energy_pion_gt       = 0;
	double Energy_kaon_gt       = 0;
	double Energy_contamination_gt = 0;

    //Note Initializations: If you want to count a total number of particles in all the events initialized outside of the event loop, if you want to count particles per event only, initiialized inside the event loop.

    cout << "------------------------ Starting new Event -------------------------------- " << endl;
    for(int i(0); i < nEvents; ++i ) {
        // Read the next entry from the tree.
        tree->GetEntry(i);
        trys=0;
        trys1=0;

        int counter_p_gt = 0; //grey tracks per event
        int counter_p = 0;//final protons
        int gt_counter          = 0; //grey tracks all particles per event
		int counter_pions = 0;//grey tracks pions
		int counter_kaons = 0;
		int counter_contamination = 0;
        int gt_protons_forw = 0,        gt_protons_zero = 0,          gt_protons_back = 0;
        int negative_nTracks_target= 0, negative_nTracks_central = 0, negative_nTracks_projectile = 0;
        int positive_nTracks_target= 0 ,positive_nTracks_central = 0, positive_nTracks_projectile = 0;
        int charged_nTracks_target = 0, charged_nTracks_central  = 0, charged_nTracks_projectile  = 0;
        int negative_tracks = 0,        positive_tracks = 0,          charged_tracks = 0;
        int tracks_forward_charged   = 0;
        int tracks_backward_charged  = 0;
     	int negative_nhadrons_target = 0, negative_nhadrons_central = 0, negative_nhadrons_projectile = 0;
        int negative_hadrons = 0; //tracks only negative hadrons

        pzlep        = event->pzlep;
        pztarg       = event->pztarg;
        struck_nucleon  = event->nucleon;
        MASS_NUCLEON = MASS_PROTON;
        if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;
        trueQ2       = event->GetTrueQ2();
        trueW2       = event->GetTrueW2();
        trueX        = event->GetTrueX();
        trueY        = event->GetTrueY();
        trueNu       = event->nu;
        s_hat        = event->GetHardS();
        t_hat        = event->t_hat;
        u_hat        = event->GetHardU();
        photon_flux  = event->GetPhotonFlux();
        event_process   = event->GetProcess();
        nParticles      = event->GetNTracks();       

        impact_parameter = event->b;
        Tb               = event->Thickness;
        distance         = event->d1st; //average distance in the medium of the struck quark 
        N_nevap             = event->Nnevap;
        N_pevap             = event->Npevap;

		counter_event++;
		x_bj = trueQ2/(2*MASS_NUCLEON* trueNu);
        //event cuts
		//DIS process for pythia
        if( event_process == 99 ){
			count_process_99++;
			if(x_bj<0.02){
				count_process_99_shad++;
			}
			if(x_bj>0.02){
				count_process_99_non_shad++;
			}
			if(x_bj>0.002 && x_bj<0.3){
				count_process_99_between++;
			}				
		}
		//QCD compton
		if( event_process ==131 || event_process ==132){
			count_process_131_132++;
			if(x_bj<0.02){
                count_process_131_132_shad++;
            }
            if(x_bj>0.02){
                count_process_131_132_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                count_process_131_132_between++;
            }
		}
		//photon-gluon fusion
		if( event_process ==135 || event_process ==136){
			count_process_135_136++;
			if(x_bj<0.02){
                count_process_135_136_shad++;
            }
            if(x_bj>0.02){
                count_process_135_136_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                count_process_135_136_between++;
            }
		}
		//vector meson production
		if(event_process ==91 || event_process ==92 || event_process==93 || event_process ==94){
			count_process_91_94++;
			if(x_bj<0.02){
                count_process_91_94_shad++;
            }
            if(x_bj>0.02){
                count_process_91_94_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                count_process_91_94_between++;
            }
		
		}
		//Resolved processes
		if(event_process ==11 || event_process ==12 || event_process==13 || event_process ==28 ||event_process ==53 || event_process ==68){
            count_process_11_13++;
            if(x_bj<0.02){
                count_process_11_13_shad++;
            }
            if(x_bj>0.02){
                count_process_11_13_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                count_process_11_13_between++;
            }

        }
		//Low Pt processes
        if(event_process ==95 ){
            count_process_95++;
            if(x_bj<0.02){
                count_process_95_shad++;
            }
            if(x_bj>0.02){
                count_process_95_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                count_process_95_between++;
            }

        }
        if( trueQ2 < 1 ) continue;
        if( trueNu < 50 || trueNu > 400 )continue;
        if( trueW2 < 16 || trueW2 > 900) continue;
        //x_bj = trueQ2/(2*MASS_NUCLEON* trueNu);
        if(x_bj < 0.0001) continue;
        counter_event_cuts++; //counting events after all cuts

		//DIS process for pythia
		if( event_process == 99 ){
			cuts_count_process_99++;
		}
		//QCD compton
		if( event_process ==131 | event_process ==132){
	    	cuts_count_process_131_132++;
		}
		//photon-gluon fusion
		if( event_process ==135 | event_process ==136){
			cuts_count_process_135_136++;
		}

        //#####ACEPTANCE & WEIGHT CALCULATION

        acc = fcn->Eval(trueQ2); 
        ww = 1 / acc; //weight
 
      
        for(int j(0); j < nParticles; ++j ) {


            particle = event->GetTrack(j);
            pdg      = particle->GetPdgCode();
            status   = particle->GetStatus();
            index    = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
            pt       = particle->GetPt();
            eta      = particle->GetEta();
            phi      = particle->GetPhi();
            rap      = particle->GetRapidity();
            mass     = particle->GetM();
            theta    = particle->GetTheta(); 
            theta = theta*1000.0; //change to mrad;
            p        = particle->GetP();
            Energy   = particle->GetE();
            xf       = particle->GetXFeynman();
            px       = particle->GetPx();
            py       = particle->GetPy();
            pz       = particle->GetPz();
            parent_index = particle->GetParentIndex();
            parent_id =particle->GetParentId();
            
            
            //double event_particle = particle->GetEvent();

            p_beam.SetPxPyPzE(0.,0.,0.,sqrt(MASS_NUCLEON*MASS_NUCLEON)); //fixed target
            e_beam.SetPxPyPzE(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON)); //initial electron
            S = (e_beam+p_beam)*(e_beam+p_beam);
            
                 
            if (pdg==-13 && status==1 &&trys==0 ){

                mu_scattered.SetPxPyPzE(px,py,-pz,Energy);
                trys++;
                counter_final_e++; 

            }
            
            if (pdg == 22 && status==21 && trys1==0){

                p_gamma.SetPxPyPzE(px,py,-pz,Energy); //virtual photon
                trys1++;
                //Beta = sqrt((p_gamma.E()*p_gamma.E())+trueQ2)/(p_gamma.E()+MASS_NUCLEON);
                p_gamma_counter++;           
            }

            phi_rot1 = atan(p_gamma.Py()/p_gamma.Pz());//is in radians

            if (pdg==2112 && status==21){struck_neutron++;}

            if (pdg==2212 && status==21){struck_proton++;}          

            //particle cuts
            //if (status != 1) continue; //final state particles 

            //Initializing Particle 

            p_final.SetPxPyPzE(px,py,-pz,Energy);
            pt_lab = sqrt(p_final.Px()*p_final.Px()+p_final.Py()*p_final.Py());            
            Beta = sqrt((p_gamma.E()*p_gamma.E())+trueQ2)/(p_gamma.E()+MASS_NUCLEON);
            gamma = 1/sqrt(1-(Beta*Beta));            
            
            //LAB FRAME 
           
            if (pdg==2212 && status==1)
            {
        
                rap_lab = -rap;//putting the minus just to keep direction in positive z
                rap_lab_computed = 0.5*log((p_final.E()+p_final.Z())/(p_final.E()-p_final.Z()));
                p_final_pz = (p_final.Pz());
                h_pz_proton->Fill(p_final_pz);             
                protons_lab++;
                h_xf->Fill(xf);
                h_rap_lab->Fill(rap_lab);
                h_rap_lab_computed->Fill(rap_lab_computed);
                
            }
			//--------------partons from PyQM---------
			if (status ==1 | status ==2 | status ==3){
				if(pdg==1 | pdg==2 | pdg==3 | pdg ==4| pdg==-1| pdg==-2| pdg==-3| pdg==-4 |pdg==21){
					h_Energy_parton->Fill(p_final.E());
				}
			}

            if (pdg==2112) {neutrons_lab++; } 
            //cout << "virtual photon lab:(Px,Py,Pz,E) =" << p_gamma.Px() << " , "<< p_gamma.Py()<< " , " << p_gamma.Pz() << p_gamma.E() << endl;                     

            h_Energy->Fill(p_final.E());
            h_pt->Fill(pt_lab);
            h_pz->Fill(p_final.Pz());
            h_mom->Fill(p_final.P());

            rot_x[0][0] = 1;
            rot_x[1][1] = cos(phi_rot1);
            rot_x[1][2] = -sin(phi_rot1);
            rot_x[2][1] = sin(phi_rot1);
            rot_x[2][2] = cos(phi_rot1);

            p_final_rotated = rot_x*(p_final.Vect());
            p_gamma_rotated = rot_x*(p_gamma.Vect());
            theta_rot2 = atan(-p_gamma_rotated.Px()/p_gamma_rotated.Pz());//is in radians

            rot_z[0][0] = cos(theta_rot2);//is a rotation in y
            rot_z[0][2] = sin(theta_rot2);
            rot_z[2][0] = -sin(theta_rot2);
            rot_z[2][2] = cos(theta_rot2);
            rot_z[1][1] = 1;
           
            p_final_rotated = rot_z*(p_final_rotated);
            p_gamma_rotated_2 = rot_z*(p_gamma_rotated);

            p_final.SetPxPyPzE(p_final_rotated.X(),p_final_rotated.Y(),p_final_rotated.Z(),Energy);


            TVectorD p_final_boost(4);
            TVectorD p_final_boost_after(4);

            p_final_boost(0) = p_final.Px();
            p_final_boost(1) = p_final.Py();
            p_final_boost(2) = p_final.Pz();
            p_final_boost(3) = p_final.E();
            p_gamma_boosted(0) = p_gamma_rotated_2(0);
            p_gamma_boosted(1) = p_gamma_rotated_2(1);
            p_gamma_boosted(2) = p_gamma_rotated_2(2);
            p_gamma_boosted(3) = p_gamma.E();           

            boost_z[0][0] = 1;
            boost_z[1][1] = 1;
            boost_z[2][2] = gamma;
            boost_z[2][3] = -gamma*Beta;
            boost_z[3][2] = -gamma*Beta;
            boost_z[3][3] = gamma;

            p_final_boost_after = boost_z*p_final_boost;
            p_gamma_boosted_ = boost_z*p_gamma_boosted;
            p_final.SetPxPyPzE(p_final_boost_after(0),p_final_boost_after(1),p_final_boost_after(2),p_final_boost_after(3));
            //cout << "virtual photon cms:(Px,Py,Pz,E) =" << p_gamma_boosted_(0) << " , "<< p_gamma_boosted_(1)<< " , " << p_gamma_boosted_(2) << p_gamma_boosted_(3) << endl;
            if (pdg==-13 && status==1 ){
            //cout << "Scattered e :(Px,Py,Pz,E) =" << mu_scattered.Px() << " , "<< mu_scattered.Py() << " , " << mu_scattered.Pz()  << mu_scattered.E() << endl;
            }

            //----------- CMS FRAME -----------------------------------------
			//----------- Partons ---------------------------
	

            if (status != 1) continue; //final state particles 

            pt_cms = sqrt(p_final.Px()*p_final.Px()+p_final.Py()*p_final.Py());//***            
            h_p_cms->Fill(p_final.P());
            h_Energy_cms->Fill(p_final.E());
            h_pz_cms->Fill(p_final.Pz());
            h_pt_cms->Fill(pt_cms);
            h_theta->Fill(theta); 
                        

            //Rapidity in the new frame
            
            rap_cms = 0.5*log((p_final.E()+p_final.Pz())/(p_final.E()-p_final.Pz()));//signos cambiados en z

            h_rap_cms->Fill(rap_cms);

            if (rap_cms < -1){gt_protons_back++;}                         
            
            if (rap_cms > -0.5 && rap_cms < 0.5 ){gt_protons_zero++;}
                            
            if (rap_cms > 2){gt_protons_forw++;}

            //##################################### Other Rapidity cuts ###########################################

            if (rap_cms > 0){
                if(pdg ==2212 || pdg==211 || pdg==321 || pdg==11 || pdg==-211 || pdg==-321 || pdg==-13 || pdg==13){

                    tracks_forward_charged++;  
                }   
            }
            if (rap_cms < 0){
                if(pdg ==2212 || pdg==211 || pdg==321 || pdg==11 || pdg==-211 || pdg==-321 || pdg==-13 || pdg==13){

                    tracks_backward_charged++;  
                }
            }
            // average mom charged particles
            if (pdg ==2212 || pdg==211 || pdg==321 || pdg==11 || pdg==-211 || pdg==-321 || pdg==-13 || pdg==13)
            {
                h_mom_charged->Fill(p);
            }

            //############ Rapidity cuts by charged, positive and negative particles ###############



           if (pdg == -13 || pdg == -211 || pdg==-321 || pdg==11 || pdg==-2212)
           {
                negative_tracks++;
				negative_total++;

                if (rap_cms < -1){
                    negative_nTracks_target++;
                }                       
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){
                    negative_nTracks_central++;
                }
                            
                if (rap_cms > 2){
                    negative_nTracks_projectile++;
                }

           }
            
           if (pdg == 2212 || pdg == 211 || pdg==321)
           //if(pdg==2212)
           {
                positive_tracks++;
				positive_total++;
                if (rap_cms < -1){
                    positive_nTracks_target++;
                }                         
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){
                    positive_nTracks_central++;
                }
                            
                if (rap_cms > 2){
                    positive_nTracks_projectile++;
                }
           }

           if (pdg == 2212 || pdg == 211 || pdg ==11 || pdg ==-211 || pdg==321 || pdg==-321 || pdg==-13 || pdg==13 || pdg==-2212)
           {
                charged_tracks++;
				charged_total++;
                if (rap_cms < -1){
                    charged_nTracks_target++;
                }                         
            
                if (rap_cms > -0.5 && rap_cms < 0.5 ){
                    charged_nTracks_central++;
                }
                            
                if (rap_cms > 2){
                    charged_nTracks_projectile++;
                }

                h_charged_Tracks->Fill(charged_tracks);
                
           }
		   //##################################################################################

		   if (pdg==-211 || pdg ==-321 || pdg==-2212)
		   {
		    negative_hadrons++;
			negative_hadrons_total++;
		   		if (rap_cms < -1){
		    		negative_nhadrons_target++;
		   		}   

		   		if (rap_cms > -0.5 && rap_cms < 0.5 ){
		   			negative_nhadrons_central++;
		   		}
		   		if (rap_cms > 2){
		   			negative_nhadrons_projectile++;
		   		}
		   }
           //###################################################################################

            if (pdg==2212)
            {
                counter_p++;
                xf_cms = 2*p_final.Pz()/sqrt(trueW2);             
                h_nTracks_proton->Fill(counter_p);
                cms_protons_xf_cut++;


                    if (p>0.2 && p < 0.6 ){

                        counter_p_gt++;    //counts the grey tracks per event
						counter_total_ng++; //counts the total of grey tracks in all events
                        Energy_proton_gt = Energy;
						h_Energy_proton_gt->Fill(Energy_proton_gt);
                        h_pt_proton_gt->Fill(pt);
                        h_rapidity_gt->Fill(rap);
                        h_nTracks_p_gt->Fill(counter_p_gt); 

                    }

            h_Energy_proton->Fill(Energy);
            h_pt_proton->Fill(pt);
            h_pz_proton_cms->Fill(p_final.Pz());         
            h_xf_cms->Fill(xf_cms);
            }

																						            
            if (pdg == 211){ counter_pos_pions++;} 
            if (pdg == -211) {counter_neg_pions++;}

            
            if (p > 0.2 && p < 0.6) //Grey Tracks All Particles
            {
                gt_counter++;
				if (pdg==211){
					Energy_pion_gt = Energy;
					counter_pions++;
					counter_total_pions++;
				} 
				if (pdg==321){
					Energy_kaon_gt = Energy;
					counter_kaons++;
					counter_total_kaons++;
				}
				if (pdg==2212 || pdg==211 || pdg==321){
					Energy_contamination_gt= Energy;
					counter_contamination++; //counts grey tracks with contamination per event
					counter_total_contamination++; //counts grey tracks with contamination for the total of events
				}
                h_pt_gt->Fill(pt);
                h_Energy_gt->Fill(Energy);
                        
            }

            //############### Particles without the Pt Cut ###########

            S_boosted = (e_beam+p_beam)*(e_beam+p_beam);

               
            if (pdg==11 && status==1  )
            {

                W2_boosted = (p_beam+p_gamma)*(p_beam+p_gamma);
           
            }

            //############### Fill histograms particle loop  ###########
                    
            counter_total_particles++;
            
              
            h_pdg->Fill(pdg);
            h_Status->Fill(status);         
            h_eta->Fill(eta);
            
            h_rapidity->Fill(rap);        

        //    prof_charged_tracks_rapidity->Fill(rap_cms,charged_tracks);
        //    prof_positive_tracks_rapidity->Fill(rap_cms,positive_tracks);
        //    prof_negative_tracks_rapidity->Fill(rap_cms,negative_tracks);
		
               

        }  
        //########### end of particle loop ###########
       

        // Q_hadronic_charge = positive_tracks - negative_tracks;
        // Q_hadronic_charge_backward=tracks_backward_positive - tracks_backward_negative;
        // Q_hadronic_charge_forward=tracks_forward_positive - tracks_forward_negative;

        //########### Fill histograms and Prof ###########

        // prof_Q_hadronic_charge_gt->Fill(counter_p_gt,Q_hadronic_charge);
        // prof_Q_hadronic_charge_backward_gt->Fill(counter_p_gt,Q_hadronic_charge_backward);
        // prof_Q_hadronic_charge_forward_gt->Fill(counter_p_gt,Q_hadronic_charge_forward);
        prof_charged_tracks_rapidity->Fill(rap_cms,charged_tracks);
        prof_positive_tracks_rapidity->Fill(rap_cms,positive_tracks);
        prof_negative_tracks_rapidity->Fill(rap_cms,negative_tracks);            
        prof_charged_nong_rapidity->Fill(rap_cms,charged_tracks-counter_p_gt);
		prof_positive_nong_rapidity->Fill(rap_cms,positive_tracks-counter_p_gt);
		prof_negative_nong_rapidity->Fill(rap_cms,negative_tracks-counter_p_gt);
        
        h_nTracks->Fill(nParticles);
        h_b->Fill(impact_parameter);
        h_trueQ2->Fill(trueQ2);
        h_trueQ2_weighted->Fill(trueQ2,ww);
         
        h_nTracks_pos_pions->Fill(counter_pos_pions);

        h_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);
        prof_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);

        prof_tracks_forward_charged->Fill(counter_p_gt,tracks_forward_charged);
        prof_tracks_backward_charged->Fill(counter_p_gt,tracks_backward_charged);

        prof_nTracks_target_positive->Fill(counter_p_gt,positive_nTracks_target);
        prof_nTracks_central_positive->Fill(counter_p_gt,positive_nTracks_central);
        prof_nTracks_projectile_positive->Fill(counter_p_gt,positive_nTracks_projectile);
        prof_nTracks_target_negative->Fill(counter_p_gt,negative_nTracks_target);
        prof_nTracks_central_negative->Fill(counter_p_gt,negative_nTracks_central);
        prof_nTracks_projectile_negative->Fill(counter_p_gt,negative_nTracks_projectile);
        prof_nTracks_target_charged->Fill(counter_p_gt,charged_nTracks_target);
        prof_nTracks_central_charged->Fill(counter_p_gt,charged_nTracks_central);
        prof_nTracks_projectile_charged->Fill(counter_p_gt,charged_nTracks_projectile);
        
		prof_nTracks_target_positive_contaminated->Fill(counter_contamination,positive_nTracks_target);
        prof_nTracks_central_positive_contaminated->Fill(counter_contamination,positive_nTracks_central);
        prof_nTracks_projectile_positive_contaminated->Fill(counter_contamination,positive_nTracks_projectile);
        prof_nTracks_target_negative_contaminated->Fill(counter_contamination,negative_nTracks_target);
        prof_nTracks_central_negative_contaminated->Fill(counter_contamination,negative_nTracks_central);
        prof_nTracks_projectile_negative_contaminated->Fill(counter_contamination,negative_nTracks_projectile);
        prof_nTracks_target_charged_contaminated->Fill(counter_contamination,charged_nTracks_target);
        prof_nTracks_central_charged_contaminated->Fill(counter_contamination,charged_nTracks_central);
        prof_nTracks_projectile_charged_contaminated->Fill(counter_contamination,charged_nTracks_projectile);

        prof_nTracks_target_positive_weighted->Fill(counter_p_gt,positive_nTracks_target,ww);
        prof_nTracks_central_positive_weighted->Fill(counter_p_gt,positive_nTracks_central,ww);
        prof_nTracks_projectile_positive_weighted->Fill(counter_p_gt,positive_nTracks_projectile,ww);
        prof_nTracks_target_negative_weighted->Fill(counter_p_gt,negative_nTracks_target,ww);
        prof_nTracks_central_negative_weighted->Fill(counter_p_gt,negative_nTracks_central,ww);
        prof_nTracks_projectile_negative_weighted->Fill(counter_p_gt,negative_nTracks_projectile,ww);
        prof_nTracks_target_charged_weighted->Fill(counter_p_gt,charged_nTracks_target,ww);
        prof_nTracks_central_charged_weighted->Fill(counter_p_gt,charged_nTracks_central,ww);
        prof_nTracks_projectile_charged_weighted->Fill(counter_p_gt,charged_nTracks_projectile,ww);

        h_nTracks_gt->Fill(counter_p_gt); 
        h_negative_nTracks_target->Fill(negative_nTracks_target);
        h_negative_nTracks_central->Fill(negative_nTracks_central);
        h_negative_nTracks_projectile->Fill(negative_nTracks_projectile);
        h_positive_nTracks_projectile->Fill(positive_nTracks_projectile);
        h_positive_nTracks_central->Fill(positive_nTracks_central);
        h_positive_nTracks_target->Fill(positive_nTracks_target);
        h_charged_nTracks_target->Fill(charged_nTracks_target);
        h_charged_nTracks_central->Fill(charged_nTracks_central);
        h_charged_nTracks_projectile->Fill(charged_nTracks_projectile);

        h_positive_Tracks->Fill(positive_tracks);
        h_negative_Tracks->Fill(negative_tracks);
        h_charged_Tracks->Fill(charged_tracks);
        h_charged_nong_Tracks->Fill(charged_tracks-counter_p_gt);
        h_positive_nong_Tracks->Fill(positive_tracks-counter_p_gt);
        h_negative_nong_Tracks->Fill(negative_tracks-counter_p_gt);

        prof_distance_grey_tracks->Fill(counter_p_gt,distance);
        prof_distance_grey_tracks_pions->Fill(counter_pions,distance);
		h_average_distance->Fill(distance);
		h_dist_ng_pions->Fill(counter_pions,distance);
        h_dist_ng->Fill(counter_p_gt,distance);
        h_impact_parameter->Fill(impact_parameter);
        h_thickness->Fill(Tb);
	
    	h_Energygt_gt->Fill(counter_p_gt,Energy_proton_gt);    
		prof_Energygt_gt->Fill(counter_p_gt,Energy_proton_gt);

        h_w2->Fill(trueW2);
        h_nu->Fill(trueNu);

        h_x_bj->Fill(x_bj);
		h_1_over_x_bj->Fill(1/x_bj);

    }// END EVENT LOOP
	cout << "---------------------------------------------------------------------------------" << endl;
	cout << "PROCESS MODELS INFO PYTHIA                     " << endl;
	cout << "Total Events							: " << counter_event << endl;
	cout << "Total Events Process 99				: " << count_process_99 << endl;
	cout << "Total Events Processes 131 - 132		: " << count_process_131_132 << endl;
	cout << "Total Events Processes 135 - 136		: " << count_process_135_136 << endl;
	cout << "Total Events Process 91-94				: " << count_process_91_94 << endl;
    cout << "Total Events Processes Resolved		: " << count_process_11_13 << endl;
    cout << "Total Events Processes 95				: " << count_process_95 << endl;
	cout << "--------------------------------Region x < 0.02---------------------------------------" << endl;
    
    cout << "Events Process 99 shadowing reg			: " << count_process_99_shad << endl;
    cout << "Events Processes 131 - 132 shadowing reg	: " << count_process_131_132_shad << endl;
    cout << "Events Processes 135 - 136 shadowing reg	: " << count_process_135_136_shad << endl;
    cout << "Events Process 91-94 shadowing reg			: " << count_process_91_94_shad << endl;
    cout << "Events Processes 11 - 13 shadowing reg		: " << count_process_11_13_shad << endl;
    cout << "Events Processes 95  shadowing reg			: " << count_process_95_shad << endl;
    cout << "-------------------------------Region x > 0.02 ---------------------------------------" << endl;
    cout << "Events Process 99 No shadowing reg			: " << count_process_99_non_shad << endl;
    cout << "Events Processes 131 - 132 No shadowing reg: " << count_process_131_132_non_shad << endl;
    cout << "Events Processes 135 - 136 No shadowing reg: " << count_process_135_136_non_shad << endl;
	cout << "Events Process 91-94 No shadowing reg		: " << count_process_91_94_non_shad << endl;
    cout << "Events Processes 11-13 No shadowing reg	: " << count_process_11_13_non_shad << endl;
    cout << "Events Processes 95 No shadowing reg		: " << count_process_95_non_shad << endl;
	cout << "-------------------------------Region  x>0.002 and x<0.3 ---------------------------------" << endl;
    cout << "Events Process 99 between				: " << count_process_99_between << endl;
    cout << "Events Processes 131 - 132 between		: " << count_process_131_132_between << endl;
    cout << "Events Processes 135 - 136 between		: " << count_process_135_136_between << endl;
    cout << "Events Process 91-94 between			: " << count_process_91_94_between << endl;
    cout << "Events Processes 11-13 between			: " << count_process_11_13_between << endl;
    cout << "Events Processes 95 between			: " << count_process_95_between << endl;

	cout << "----------------------------------------------------------------------------------------" << endl;
	cout << "AFTER E665 EVENTS CUTS                       : " << endl;
	cout << "After E665 cuts - Events Process 99          : " << cuts_count_process_99 << endl;
    cout << "After E665 cuts - Events Processes 131 - 132 : " << cuts_count_process_131_132 << endl;
    cout << "After E665 cuts - Events Processes 135 - 136 : " << cuts_count_process_135_136 << endl;
	cout << "Ratios Events Process 99                     : " << (cuts_count_process_99*100/counter_event) << endl;
    cout << "Ratios Events Processes 131 - 132 			  : " << (cuts_count_process_131_132*100/counter_event) << endl;
    cout << "Ratios Events Processes 135 - 136 			  : " << (cuts_count_process_135_136*100/counter_event) << endl;

	cout << "------------------------------------------------------------------------------------------" << endl;

    cout << "Number of particles                          : " << counter_total_particles << endl;
    cout << "Events after event cuts                      : " << counter_event_cuts << endl;
    cout << "mean value of All protons                    : " << h_nTracks_proton->GetMean() << endl;
    cout << "number scattered electrons (1 per event)     : " << counter_final_e << endl;
    cout << "number virtual photons (1 per event)         : " << p_gamma_counter << endl;
    cout << "number struck nucleons(1 per event)          : " << struck_proton+struck_neutron << endl;
    cout << "number struck neutron                        : " << struck_neutron <<" percentage: " << (100*struck_neutron)/(struck_proton+struck_neutron)<< endl;
    cout << "number struck proton                         : " << struck_proton << " percentage: "<<(100*struck_proton)/(struck_proton+struck_neutron) << endl;
    cout <<"--------------------------- GREY TRACKS 0.2 < p < 0.6 GeV ------------------------------------" << endl;
    
	cout << "Total grey tracks for all the events (protons) 			: " << counter_total_ng << endl;
	cout << "Total grey tracks with contamination for all the events 	: " << counter_total_contamination << endl;
	cout << "Total grey tracks pions									: " << counter_total_pions << endl;
	cout << "Total grey tracks kaons									: " << counter_total_kaons << endl;	
	cout << "--------------------------------------------------------------------------------------------- " << endl;
	cout << "number of total protons lab frame & final state: " << protons_lab << endl;
    cout << "number protons xf cut                          : " << cms_protons_xf_cut <<endl;
    // cout << "number of protons rap < -1: " << gt_protons_back << endl;
    // cout << "number of protons -0.5 < rap < 0.5 : " << gt_protons_zero << endl;
    // cout << "number of protons rap > 2: " << gt_protons_forw << endl;
    cout << "--------------------------------------------------------------------------------------------- " << endl;
	cout << "Total charged particles		 				 :" << charged_total << endl;
	cout << "Total positive particles                        :" << positive_total << endl;
    cout << "Total negative particles                        :" << negative_total << endl;
	cout << "Total negative hadrons particles                :" << negative_hadrons_total << endl;
	
	cout << "--------------------------------------------------------------------------------------------- " << endl;
	cout << "W2 after:  " << W2_boosted << endl;

    cout << "mean value of CHARGED PARTICLES TARGET         : "  <<h_charged_nTracks_target->GetMean() << endl;
    cout << "mean value of CHARGED PARTICLES CENTRAL        : "  <<h_charged_nTracks_central->GetMean() << endl;
    cout << "mean value of CHARGED PARTICLES PROJECTILE     : "  <<h_charged_nTracks_projectile->GetMean() << endl;

    cout << "mean value of POSITIVE PARTICLES TARGET        : "  <<h_positive_nTracks_target->GetMean() << endl;
    cout << "mean value of  POSITIVE PARTICLES CENTRAL      : "  <<h_positive_nTracks_central->GetMean() << endl;
    cout << "mean value of  POSITIVE PARTICLES PROJECTILE   : "  <<h_positive_nTracks_projectile->GetMean() << endl;

    cout << "mean value of NEGATIVE PARTICLES TARGET        : "  <<h_negative_nTracks_target->GetMean() << endl;  
    cout << "mean value of NEGATIVE PARTICLES CENTRAL       : "  <<h_negative_nTracks_central->GetMean() << endl;
    cout << "mean value of NEGATIVE PARTICLES PROJECTILE    : "  <<h_negative_nTracks_projectile->GetMean() << endl;
    

    cout << "mean value of CHARGED PARTICLES ALL REGIONS    : "  << h_charged_Tracks->GetMean() << endl;
    cout << "mean value of NEGATIVE PARTICLES ALL REGIONS   : "  << h_negative_Tracks->GetMean() << endl;
    cout << "mean value of POSITIVE PARTICLES ALL REGIONS   : "  << h_positive_Tracks->GetMean() << endl;
    cout << "mean value of CHARGED  NO NG PARTICLES ALL REGIONS   : "  << h_charged_nong_Tracks->GetMean() << endl;
    cout << "mean value of POSITIVE NO NG PARTICLES ALL REGIONS   : "  << h_positive_nong_Tracks->GetMean() << endl;
    cout << "mean value of NEGATIVE NO NG PARTICLES ALL REGIONS   : "  << h_negative_nong_Tracks->GetMean() << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Average Momentum of CHARGED Particles          : "  << h_mom_charged->GetMean() <<endl; 
    //------------------------------- Scaling profiles--------------------------------------------------
     TAxis *xaxis = prof_nTracks_target_positive->GetXaxis();
     int nbins = xaxis->GetNbins();
     float xmin = xaxis->GetXmin();
     float xmax =xaxis->GetXmax();
     TH1F *hp = new TH1F("hp","hp",nbins,xmin,xmax);
     for (int bin = 0; bin <= nbins; bin++)
     {
         hp->SetBinContent(bin,prof_nTracks_target_positive->GetBinContent(bin));
         hp->SetBinError(bin,prof_nTracks_target_positive->GetBinError(bin));
     }
     hp->SetEntries(prof_nTracks_target_positive->GetEntries());
     hp->Scale(2);
     hp->Draw();

    //--------------------------------------------------------------------------------------------------


    //################################### ACCEPTANCE TREATEMENT###########################
    TCanvas *ac = new TCanvas("ac","Data Acceptance", 600,600);
    ac->Divide(2,1);
    ac->cd(1);
    //TGraph *tg_acc_1 = new TGraph(82,acc_x,acc_y);
    gPad->SetLogx();
    tg_acc_1->SetTitle("TSpline; Q^{2} ; Acceptance");
    tg_acc_1->GetXaxis()->SetLimits(1,100);
    tg_acc_1->SetMaximum(2.0);//for y axis
    tg_acc_1->SetMinimum(0.0);
    tg_acc_1->Draw("AP*");
    //TSpline5 *s = new TSpline5("tg3_s",tg_acc_1);
    s->SetLineColor(kRed);
    //s->Eval(3);
    //cout << "findX: " << s->Eval(3)<<endl;
    s->Draw("same");
    s->Print();
    TLegend *l_acc_1 = new TLegend(0.35, 0.75, 0.4, 0.92);
    l_acc_1->SetBorderSize(0); // no border
    l_acc_1->SetFillStyle(0);
    l_acc_1->SetFillColor(0); // Legend background should be white
    l_acc_1->SetTextFont(32);
    l_acc_1->SetTextSize(0.08); // Increase entry font size!
    l_acc_1->AddEntry(tg_acc_1,"TSpline");
    l_acc_1->Draw();

    ac->cd(2);
    //TGraph *tg_acc_2 = new TGraph(82,acc_x,acc_y);
    //TF1* fcn = new TF1("fcn", "[0] / (1 + TMath::Exp(-[1] * (x - [2])))", 1, 100);
    gPad->SetLogx();
    tg_acc_2->SetTitle("Exponential; Q^{2} ; Acceptance");
    tg_acc_2->GetXaxis()->SetLimits(1,100);
    tg_acc_2->SetMaximum(2.0);//for y axis
    tg_acc_2->SetMinimum(0.0);
    tg_acc_2->Draw("AP*");
    //tg_acc_2->Fit("fcn");
    // auto acc = fcn->Eval(7.2);
    // auto ww = 1 / acc;

    // cout << "Acceptance = " << acc << endl;
    // cout << "weight = " << ww << endl;
    // for (int i = 0; i < 82; ++i)
    // {
    //     auto acc = fcn->Eval(acc_x[i]); //Q^2 from paper
    //     auto ww = 1 / acc; //weight

    //     //tg_acc_2->Fit(acc_x[i],ww);
    // }

    TLegend *l_acc_2 = new TLegend(0.35, 0.75, 0.4, 0.92);
    l_acc_2->SetBorderSize(0); // no border
    l_acc_2->SetFillStyle(0);
    l_acc_2->SetFillColor(0); // Legend background should be white
    l_acc_2->SetTextFont(32);
    l_acc_2->SetTextSize(0.08); // Increase entry font size!
    l_acc_2->AddEntry(tg_acc_2,"exp^{(-1(x - 2))}");
    l_acc_2->Draw();

    //################################### Drawing for e - Xe #############################
    //---------------------------------------- Multiplicity Ratios by Region ----------------------------------------------------
   // TCanvas *p0 = new TCanvas("p0","Multiplicty Ratio by region");
    //p0->Divide(2,1);
    //p0->cd(1);
    //prof_tracks_forward_charged->SetTitle(";Forward charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_tracks_forward_charged->Scale(1.0/scale_charged_forward);
    //prof_tracks_forward_charged->Draw();
    //p0->cd(2);
    //prof_tracks_backward_charged->SetTitle(";Backward charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_tracks_backward_charged->Scale(1.0/scale_charged_backward);
    //prof_tracks_backward_charged->Draw();



    //TCANvas p2,p3,p4 nTracks to be analized
    //---------------------------------------- Multiplicity Ratios by charged and Region ----------------------------------------------
    TCanvas *p4 = new TCanvas("p4","Multipl Ratio by Region and Charged part");
    p4->Divide(3,1);
    p4->cd(1);
    prof_nTracks_target_charged->SetTitle("Target charged;Target charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_target_charged->Scale(1.0/2.12425); //1.94
    prof_nTracks_target_charged->Scale(1.0/scale_charged_target); 
    prof_nTracks_target_charged->Draw();
    //TFitResultPtr r = prof_nTracks_target_charged->Fit("pol1","S");
    //r->Print("V");
    /*prof_nTracks_target_charged->Fit("pol1");
    TF1 *fit1 = prof_nTracks_target_charged->GetFunction("pol1");
    fit1->SetLineWidth(1);
    fit1->SetLineColor(kBlue);
    */
    //gStyle->SetLineStyleString(3,"400 200");
     

    p4->cd(2);
    prof_nTracks_central_charged->SetTitle("Central charged;Central Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_charged->Scale(1.0/1.16563); //1.69
    prof_nTracks_central_charged->Scale(1.0/scale_charged_central); 
    prof_nTracks_central_charged->Draw();
    /*
    prof_nTracks_central_charged->Fit("pol1");
    TF1 *fit2 = prof_nTracks_central_charged->GetFunction("pol1");
    fit2->SetLineWidth(1);
    fit2->SetLineColor(kBlue);
    */

    p4->cd(3);
    prof_nTracks_projectile_charged->SetTitle("Projectile charged;Projectile Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_charged->Scale(1.0/1.34466); //1.55
    prof_nTracks_projectile_charged->Scale(1.0/scale_charged_projectile); //1.55
    prof_nTracks_projectile_charged->Draw();
    /*
    prof_nTracks_projectile_charged->Fit("pol1");
    TF1 *fit3 = prof_nTracks_projectile_charged->GetFunction("pol1");
    fit3->SetLineWidth(1);
    fit3->SetLineColor(kBlue);
    */

    TCanvas *p3 = new TCanvas("p3","Multipl Ratio by Region and Negative part");
    p3->Divide(3,1);
    
    p3->cd(1);
    prof_nTracks_target_negative->SetTitle("Target negative;Target - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_target_negative->Scale(1.0/0.617756); 
    prof_nTracks_target_negative->Scale(1.0/scale_negative_target); 
    prof_nTracks_target_negative->Draw();
    /*
    prof_nTracks_target_negative->Fit("pol1");
    TF1 *fit4 = prof_nTracks_target_negative->GetFunction("pol1");
    fit4->SetLineWidth(1);
    fit4->SetLineColor(kBlue);
    */

    p3->cd(2);
    prof_nTracks_central_negative->SetTitle("Central negative;Central Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_negative->Scale(1.0/0.559972);
    prof_nTracks_central_negative->Scale(1.0/scale_negative_central);
    prof_nTracks_central_negative->Draw();
    /*
    prof_nTracks_central_negative->Fit("pol1");
    TF1 *fit5 = prof_nTracks_central_negative->GetFunction("pol1");
    fit5->SetLineWidth(1);
    fit5->SetLineColor(kBlue);
    */
    p3->cd(3);
    prof_nTracks_projectile_negative->SetTitle("Projectile negative;Projectile Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_negative->Scale(1.0/1.013);
    prof_nTracks_projectile_negative->Scale(1.0/scale_negative_projectile);
    prof_nTracks_projectile_negative->Draw();
    /*
    prof_nTracks_projectile_negative->Fit("pol1");
    TF1 *fit6 = prof_nTracks_projectile_negative->GetFunction("pol1");
    fit6->SetLineWidth(1);
    fit6->SetLineColor(kBlue);
    */

    TCanvas *p2 = new TCanvas("p2","Multipl Ratio by Region and Positive part");
    p2->Divide(3,1);
    
    p2->cd(1);
    prof_nTracks_target_positive->SetTitle("Target positive;Target Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}"); 
    //prof_nTracks_target_positive->Scale(1.0/1.5065);
    prof_nTracks_target_positive->Scale(1.0/scale_positive_target);
    prof_nTracks_target_positive->Draw();
    /*
    prof_nTracks_target_positive->Fit("pol1");
    TF1 *fit7 = prof_nTracks_target_positive->GetFunction("pol1");
    fit7->SetLineWidth(1);
    fit7->SetLineColor(kBlue);
    */

    p2->cd(2);
    prof_nTracks_central_positive->SetTitle("Central positive;Central Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_positive->Scale(1.0/0.605654);
    prof_nTracks_central_positive->Scale(1.0/scale_positive_central);
    prof_nTracks_central_positive->Draw();
    /*
    prof_nTracks_central_positive->Fit("pol1");
    TF1 *fit8 = prof_nTracks_central_positive->GetFunction("pol1");
    fit8->SetLineWidth(1);
    fit8->SetLineColor(kBlue);
    */

    p2->cd(3);
    prof_nTracks_projectile_positive->SetTitle("Projectile positive;Projectile + Region n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_positive->Scale(1.0/0.331667);
    prof_nTracks_projectile_positive->Scale(1.0/scale_positive_projectile);
    prof_nTracks_projectile_positive->Draw();
    /*
    prof_nTracks_projectile_positive->Fit("pol1");
    TF1 *fit9 = prof_nTracks_projectile_positive->GetFunction("pol1");
    fit9->SetLineWidth(1);
    fit9->SetLineColor(kBlue);
    */
//----------------------------------------------------------------------------------------------
//--------------------------Multiplicity ratios vs rapidity
//-----------------------------------------------------------------------------------------------

    TCanvas *m = new TCanvas("m","Ratios vs Rapidity");
    m->Divide(3,1);
    m->cd(1);
    prof_charged_tracks_rapidity->Scale(1.0/scale_charged);
    prof_charged_nong_rapidity->Scale(1.0/scale_charged);
    prof_charged_tracks_rapidity->Draw();
    prof_charged_nong_rapidity->Draw("same");
    m->cd(2);
    prof_positive_tracks_rapidity->Scale(1.0/scale_positive);
    prof_positive_nong_rapidity->Scale(1.0/scale_positive);
    prof_positive_tracks_rapidity->Draw();
    prof_positive_nong_rapidity->Draw("same");
    m->cd(3);
    prof_negative_tracks_rapidity->Scale(1.0/scale_negative);
    prof_negative_nong_rapidity->Scale(1.0/scale_negative);
    prof_negative_tracks_rapidity->Draw();
    prof_negative_nong_rapidity->Draw("same");

 
//---------------------------------contamination -------------------------------------------------
    TCanvas *pp4 = new TCanvas("pp4","Multipl Ratio by Region and Charged contaminated");
    pp4->Divide(3,1);
    pp4->cd(1);
    prof_nTracks_target_charged_contaminated->SetTitle("Target charged;Target charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_charged_contaminated->Scale(1.0/scale_charged_target);
    prof_nTracks_target_charged_contaminated->Draw();
   
	pp4->cd(2);
    prof_nTracks_central_charged_contaminated->SetTitle("Central charged;Central Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_charged_contaminated->Scale(1.0/scale_charged_central);
    prof_nTracks_central_charged_contaminated->Draw();

    pp4->cd(3);
    prof_nTracks_projectile_charged_contaminated->SetTitle("Projectile charged;Projectile Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_charged_contaminated->Scale(1.0/scale_charged_projectile); //1.55
    prof_nTracks_projectile_charged_contaminated->Draw();
	
	TCanvas *pp3 = new TCanvas("pp3","Multipl Ratio by Region and Negative contaminated");
    pp3->Divide(3,1);
   
    pp3->cd(1);
    prof_nTracks_target_negative_contaminated->SetTitle("Target negative;Target - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_negative_contaminated->Scale(1.0/scale_negative_target);
    prof_nTracks_target_negative_contaminated->Draw();

    pp3->cd(2);
    prof_nTracks_central_negative_contaminated->SetTitle("Central negative;Central Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_negative_contaminated->Scale(1.0/scale_negative_central);
    prof_nTracks_central_negative_contaminated->Draw();
    pp3->cd(3);
    prof_nTracks_projectile_negative_contaminated->SetTitle("Projectile negative;Projectile Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_negative_contaminated->Scale(1.0/scale_negative_projectile);
    prof_nTracks_projectile_negative_contaminated->Draw();
	
	TCanvas *pp2 = new TCanvas("pp2","Multipl Ratio by Region and Positive contaminated");
    pp2->Divide(3,1);

    pp2->cd(1);
    prof_nTracks_target_positive_contaminated->SetTitle("Target positive;Target Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_positive_contaminated->Scale(1.0/scale_positive_target);
    prof_nTracks_target_positive_contaminated->Draw();

    p2->cd(2);
    prof_nTracks_central_positive_contaminated->SetTitle("Central positive;Central Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_positive_contaminated->Scale(1.0/scale_positive_central);
    prof_nTracks_central_positive_contaminated->Draw();

    p2->cd(3);
    prof_nTracks_projectile_positive_contaminated->SetTitle("Projectile positive;Projectile + Region n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_positive_contaminated->Scale(1.0/scale_positive_projectile);
    prof_nTracks_projectile_positive_contaminated->Draw();
//--------------------------------- Track analysis ------------------------------------------------
    TCanvas *a1 = new TCanvas("a1","track distribution target ",750,450);
    a1->Divide(3,1);
    a1->cd(1);
    h_positive_nTracks_target->SetTitle("track distribution target;Positive tracks - Target Region;nEvents");
    h_positive_nTracks_target->Draw("P");
    a1->cd(2);
    h_negative_nTracks_target->SetTitle("track distribution target;Negative tracks - Target Region;nEvents");
    h_negative_nTracks_target->Draw("P");
    a1->cd(3);
    h_charged_nTracks_target->SetTitle("track distribution target;Charged tracks - Target Region;nEvents ");
    h_charged_nTracks_target->Draw("P");

    //a1->SaveAs(("/Volumes/Tesla/projects/output/plots/qhat_05/track ditribution target_FILE_"+filename+".pdf").c_str()); 
    
    TCanvas *a2 = new TCanvas("a2","track distribution central",750,450);
    a2->Divide(3,1);
    a2->cd(1);
    h_positive_nTracks_central->SetTitle("track distribution central; Positive tracks - Central Region;nEvents");
    h_positive_nTracks_central->Draw("P");
    a2->cd(2);
    h_negative_nTracks_central->SetTitle(";Negative tracks - Central Region;nEvents");
    h_negative_nTracks_central->Draw("P");
    a2->cd(3);
    h_charged_nTracks_central->SetTitle(";Charged tracks - Central Region;nEvents");
    h_charged_nTracks_central->Draw("P");

    // a2->SaveAs(("/Volumes/Tesla/projects/output/plots/qhat_05/track ditribution central_FILE_"+filename+".pdf").c_str()); 


    TCanvas *a3 = new TCanvas("a3","track distribution projectile",750,450);
    a3->Divide(3,1);
    a3->cd(1);
    h_positive_nTracks_projectile->SetTitle("track distribution projectile;Positive tracks - Projectile Region;nEvents");
    h_positive_nTracks_projectile->Draw("P");
    a3->cd(2);
    h_negative_nTracks_projectile->SetTitle(";Negative tracks - Projectile Region;nEvents");
    h_negative_nTracks_projectile->Draw("P");
    a3->cd(3);
    h_charged_nTracks_projectile->SetTitle(";Charged tracks - Projectile Region;nEvents");
    h_charged_nTracks_projectile->Draw("P");

    //a3->SaveAs(("/Volumes/Tesla/projects/output/plots/qhat_05/track ditribution projectile_FILE_"+filename+".pdf").c_str()); 


    std::cout << "Attempting to write output file" << std::endl;
    out->Write();
    std::cout << "Attempting to write output file: Ok!" << std::endl;

}

