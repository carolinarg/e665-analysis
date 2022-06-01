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

map<string, TH1I*> hist_distribution_tracks;

const vector<string> GreyTracks = {"1","2","3","4","5","6","7","8"};


void runEICTree_D_statis(string option="scale_g1_qhat_05_1E8_TEST", const int nEvents = 100000){


    std::ifstream fileIN,fileINmeans;
   	TString basepath = "/gpfs/mnt/gpfs02/eic/crob/BeAGLE-dev-2021-04-27/BeAGLE/output/e665/qhat_05/sim-ene-22/muD/g1/";
	TString basename = "Output_input_temp";

        std::ofstream flux(("muD-means/"+option+".txt").c_str()); 

        TChain *tree = new TChain("EICTree");

		for(int idx=1; idx < 1+1; ++idx ){
			TString path = Form("%s/%s_%d.root",basepath.Data(),basename.Data(),idx);
			cout <<"Adding path = " << path << endl;
			tree->Add(path);
		}
      

        EventBeagle* event(NULL);
        tree->SetBranchAddress("event", &event);
    
        TFile *out = new TFile (std::string("analysis/Analysis_ratio_muD_"+option+".root").c_str(), "RECREATE");
        //####################  Creating Histograms ####################################

        //TTree *t = new TTree("t","t filter");

        for (const string& gt : GreyTracks) {
        hist_distribution_tracks[gt] = new TH1I(Form("hist_distribution_tracks%s", gt.c_str()), "title;Multiplicity Distributions;Events", 100, -0.5, 10.5);
        }
    
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
        TH1F *h_trueQ2              = new TH1F("h_trueQ2", "h_trueQ2", 100, 0, 20);
        TH1F *h_Energy_proton       = new TH1F("h_Energy_proton", "h_Energy_proton", 100, 0, 3);
        TH1F *h_Energy_proton_gt    = new TH1F("h_Energy_proton_gt", "h_Energy_proton_gt", 100, 0, 3);
        TH1F *h_xf                  = new TH1F("h_xf", "h_xf", 100, -1.1, 1.1);
        
        TH1F *h_pt_proton       = new TH1F("h_pt_proton", "h_pt_proton", 100, 0, 2);
        TH1F *h_pt_proton_gt    = new TH1F("h_pt_proton_gt", "h_pt_proton_gt", 100, 0, 2);
        TH1F *h_nTracks_gt      = new TH1F("h_nTracks_gt","h_nTracks_gt",100, 0, 50);
        TH1F *h_pt_gt           = new TH1F("h_pt_gt", "h_pt_gt", 100, 0, 1);
        TH1F *h_Energy_gt       = new TH1F("h_Energy_gt", "h_Energy_gt", 100,0, 6);
        TH1F *h_rapidity_gt     = new TH1F("h_rapidity_gt","h_rapidity_gt",100,-10,10);
        TH1F *h_nTracks_p_gt    = new TH1F("h_nTracks_p_gt","h_nTracks_p_gt", 11, -0.5, 10.5); //****
        TH1F *h_nTracks_proton  = new TH1F("h_nTracks_proton", "h_nTracks_proton", 11, -0.5, 10.5); //****
        TH1F *h_pz_proton       = new TH1F("h_pz_proton","h_pz_proton", 40,-10,10);
        TH1F *h_pz_proton_cms   = new TH1F("h_pz__proton_cms","h_pz__proton_cms", 40,-10,10);

        TH1F *h_negative_nTracks = new TH1F("h_negative_nTracks","h_negative_nTracks",24,0,12);
        TH1F *h_positive_nTracks = new TH1F("h_positive_nTracks","h_positive_nTracks",24,0,12);
        TH1F *h_charged_nTracks  = new TH1F("h_charged_nTracks","h_charged_nTracks",24,0,12);
        TH1F *h_negative_nong_nTracks = new TH1F("h_negative_nong_nTracks","h_negative_nong_nTracks",24,0,12);
        TH1F *h_positive_nong_nTracks = new TH1F("h_positive_nong_nTracks","h_positive_nong_nTracks",24,0,12);
        TH1F *h_charged_nong_nTracks  = new TH1F("h_charged_nong_nTracks","h_charged_nong_nTracks",24,0,12);

        TH1F *h_negative_nTracks_target      = new TH1F("h_negative_nTracks_target","h_negative_nTracks_target", 12, 0, 12.5);
        TH1F *h_negative_nTracks_central     = new TH1F("h_negative_nTracks_central","h_negative_nTracks_central", 12, 0, 12.5);
        TH1F *h_negative_nTracks_projectile = new TH1F("h_negative_nTracks_projectile","h_negative_nTracks_projectile", 12, 0, 12.5);
        TH1F *h_positive_nTracks_target     = new TH1F("h_positive_nTracks_target","h_positive_nTracks_target", 12, 0, 12.5);
        TH1F *h_positive_nTracks_central    = new TH1F("h_positive_nTracks_central","h_positive_nTracks_central", 12, 0, 12.5);
        TH1F *h_positive_nTracks_projectile = new TH1F("h_positive_nTracks_projectile","h_positive_nTracks_projectile", 12, 0, 12.5);
        TH1F *h_charged_nTracks_target      = new TH1F("h_charged_nTracks_target","h_charged_nTracks_target", 12, 0, 12.5);
        TH1F *h_charged_nTracks_central     = new TH1F("h_charged_nTracks_central","h_charged_nTracks_central", 12, 0, 12.5);
        TH1F *h_charged_nTracks_projectile  = new TH1F("h_charged_nTracks_projectile","h_charged_nTracks_projectile", 12, 0, 12.5);
        TH1F *h_charged_nTracks_forward     = new TH1F("h_charged_nTracks_forward","h_charged_nTracks_forward", 12, 0, 12.5);
        TH1F *h_charged_nTracks_backward  = new TH1F("h_charged_nTracks_backward","h_charged_nTracks_backward", 12, 0, 12.5);



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

        TProfile *prof_charged_tracks_rapidity      = new TProfile("prof_charged_tracks_rapidity","prof_charged_tracks_rapidity",10,-5,5);
        TProfile *prof_positive_tracks_rapidity     = new TProfile("prof_positive_tracks_rapidity","prof_positive_tracks_rapidity",10,-5,5);
        TProfile *prof_negative_tracks_rapidity     = new TProfile("prof_negative_tracks_rapidity","prof_negative_tracks_rapidity",10,-5,5);
        TProfile *prof_Q_hadronic_charge_gt         = new TProfile("prof_Q_hadronic_charge_gt","prof_Q_hadronic_charge_gt", 10,bins);
        //TProfile *prof_Q_hadronic_charge_gt_backward =new TProfile("prof_Q_hadronic_charge_gt_backward","prof_Q_hadronic_charge_gt_backward", 10,bins);
        TProfile *prof_Q_hadronic_charge_backward_gt    = new TProfile("prof_Q_hadronic_charge_backward_gt","prof_Q_hadronic_charge_backward_gt",13,bins);
        TProfile *prof_Q_hadronic_charge_forward_gt     = new TProfile("prof_Q_hadronic_charge_forward_gt","prof_Q_hadronic_charge_forward_gt",13,bins);
		TProfile *prof_Energygt_gt                      = new TProfile("prof_Energygt_gt","prof_Energygt_gt",13,bins);
		        
        TH1F *h_w2          = new TH1F("h_W2","h_W2", 100, 0, 1000);
        TH1F *h_nu          = new TH1F("h_nu", "h_nu",100,50,400);
        TH1F *h_mom         = new TH1F("h_mom", "h_mom", 60, -30, 30);
        TH1F *h_p_cms       = new TH1F("h_p_cms", "h_p_cms", 60,0,30);
        TH1F *h_Energy_cms  = new TH1F("h_Energy_cms","h_Energy_cms", 100,0,6);
        TH1F *h_pz_cms      = new TH1F("h_pz_cms","h_pz_cms",80, -40,40);
        TH1F *h_pt_cms      = new TH1F("h_pt_cms","h_pt_cms",100, 0,2);
        TH1F *h_rap_cms     = new TH1F("h_rap_cms","h_rap_cms",48,-4,4);
        TH1F *h_xf_cms      = new TH1F("h_xf_cms","h_xf_cms",100,-1.1,1.1);
        TH1F *h_x_bj        = new TH1F("h_x_bj","h_x_bj",100,0,0.15);
        TH1F *h_theta       = new TH1F("h_theta","h_theta",100,0,360);
        TH1F *h_rap_lab             = new TH1F("h_rap_lab","h_rap_lab", 48,-4,4);
        TH1F *h_rap_lab_computed    = new TH1F("h_rap_lab_computed","h_rap_lab_computed", 48,-4,4);

        
        TProfile *prof  = new TProfile("prof", "Rapidity;<ng>", 21, -7, 3);
        TH2F *hist      = new TH2F("hist","Rapidity;<ng>", 21, -4, 4,11,-0.5,10.5);
    	TH2F *h_dist_ng_pions = new TH2F("h_dist_ng_pions","h_dist_ng_pions",13,-0.5,12.5, 30, 0, 3);
    	TH2F *h_Energygt_gt = new TH2F("h_Energygt_gt", "h_Energygt_gt",13,-0.5,12.5, 30, 0, 3);
    
        int counter_total_particles = 0;
        int counter_event           = 0; 
		int counter_ALL_event       = 0;
        int counter_total_contamination = 0;
        string title_collision = "e - Xe";

        TLorentzVector e_scattered;
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
        double Q_hadronic_charge    = 0;
        double Q_hadronic_charge_backward=0;
        double Q_hadronic_charge_forward =0;

        double sum_N_prot_gt_0=0;
        double average_N_prot = 0;
        
        
        int S_boosted               = 0;
        double S                    = 0;
        double gamma                = 0;
        int cms_protons_xf_cut      = 0;
        

        //auto graph_n_vs_ng = TH1TOTGraph(prof_nTracks_Grey_tracks);


        double pzlep = 0;
        double pztarg = 0;
        int struck_nucleon = 0;
        double MASS_NUCLEON;
        double trueQ2   = 0;
        double trueW2   = 0;
        double trueX    = 0;
        double trueY    = 0;
        double trueNu   = 0;
        double s_hat    = 0;
        double t_hat    = 0;
        double u_hat    = 0;
        double photon_flux  = 0;
        int event_process   = 0;
        int nParticles      = 0;
        double impact_parameter = 0;
        double Tb       = 0;
        double distance = 0;
        int N_nevap     = 0;
        int N_pevap     = 0;
        double x_bj     = 0; 

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
        double p        = 0;
        double Energy   = 0;
        double xf       = 0;
        double px       = 0;
        double py       = 0;
        double pz       = 0;
        int parent_index = 0;
        double parent_id = 0;

        int counter           = 0;              
        //wint gt_counter        = 0;       
        int counter_pos_pions = 0;
        int counter_neg_pions = 0;                              
        double pt_cms         = 0;
        double pt_lab         = 0;
        double xf_cms_s       = 0;
        double p_final_pz     = 0;  
        trys = 0;
        trys1= 0;
    	double Energy_proton_gt     = 0;
    	double Energy_pion_gt       = 0;
    	double Energy_kaon_gt       = 0;
    	double Energy_contamination = 0;    
 		int count_process_99 = 0;
    	int count_process_131_132 = 0;
    	int count_process_135_136 = 0;
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

        int cuts_count_process_99 = 0;
        int cuts_count_process_131_132 = 0;
        int cuts_count_process_135_136 = 0;
        int cuts_cuts_count_process_99 = 0;
        int cuts_cuts_count_process_131_132 = 0;
        int cuts_cuts_count_process_135_136 = 0;
        int cuts_count_process_91_94=0;
        int cuts_count_process_11_13=0;
        int cuts_count_process_95=0;

        int cuts_count_process_99_shad = 0;
        int cuts_count_process_131_132_shad = 0;
        int cuts_count_process_135_136_shad = 0;
        int cuts_count_process_91_94_shad=0;
        int cuts_count_process_11_13_shad=0;
        int cuts_count_process_95_shad=0;
        int cuts_count_process_99_non_shad = 0;
        int cuts_count_process_131_132_non_shad = 0;
        int cuts_count_process_135_136_non_shad = 0;
        int cuts_count_process_91_94_non_shad=0;
        int cuts_count_process_11_13_non_shad=0;
        int cuts_count_process_95_non_shad=0;
        int cuts_count_process_99_between = 0;
        int cuts_count_process_131_132_between = 0;
        int cuts_count_process_135_136_between = 0;
        int cuts_count_process_91_94_between=0;
        int cuts_count_process_11_13_between=0;
        int cuts_count_process_95_between=0;
        
        cout << "------------------------ Starting new Event -------------------------------- " << endl;
        for(int i(0); i < nEvents; ++i ) {
            
            // Read the next entry from the tree.
            tree->GetEntry(i);
            trys=0;
            trys1=0;
            int counter_p_gt = 0;
            int counter_p = 0;
	        int gt_counter          = 0; //grey tracks all particles per event
        	int counter_pions = 0;//grey tracks pions
        	int counter_kaons = 0;
       	 	int counter_contamination = 0;
            int gt_protons_forw=0, gt_protons_zero=0,gt_protons_back=0;
            int negative_nTracks_target= 0, negative_nTracks_central = 0, negative_nTracks_projectile = 0;
            int positive_nTracks_target=0,positive_nTracks_central=0,positive_nTracks_projectile=0;
            int charged_nTracks_target=0,charged_nTracks_central=0,charged_nTracks_projectile=0;
            int negative_tracks=0, positive_tracks=0,charged_tracks=0;
            int tracks_forward_charged=0, tracks_backward_charged=0;

            pzlep = event->pzlep;
            pztarg = event->pztarg;
            struck_nucleon = event->nucleon;
            MASS_NUCLEON = MASS_PROTON;
            if( struck_nucleon==2112 ) MASS_NUCLEON = MASS_NEUTRON;
            trueQ2 = event->GetTrueQ2();
            trueW2 = event->GetTrueW2();
            trueX = event->GetTrueX();
            trueY = event->GetTrueY();
            trueNu = event->nu;
	        s_hat = event->GetHardS();
            t_hat = event->t_hat;
            u_hat = event->GetHardU();
            photon_flux = event->GetPhotonFlux();
            event_process = event->GetProcess();
            nParticles = event->GetNTracks();      
            impact_parameter = event->b;
            Tb = event->Thickness;
            distance = event->d1st;
            N_nevap = event->Nnevap;
            N_pevap = event->Npevap;
            counter_ALL_event++;
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

            //event cuts
	
            if( trueQ2 < 1 ) continue;
            if( trueNu < 50 || trueNu > 400 )continue;
            if( trueW2 < 16 || trueW2 > 900) continue;
            x_bj = trueQ2/(2*MASS_NUCLEON*trueNu);
            if(x_bj < 0.0001) continue;
            counter_event++; //counting events after all cuts
             //DIS process for pythia
        if( event_process == 99 ){
            cuts_count_process_99++;
            if(x_bj<0.02){
                cuts_count_process_99_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_99_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_99_between++;
            }
        }
        //QCD compton
        if( event_process ==131 || event_process ==132){
            cuts_count_process_131_132++;
            if(x_bj<0.02){
                cuts_count_process_131_132_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_131_132_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_131_132_between++;
            }
        }
        //photon-gluon fusion
        if( event_process ==135 || event_process ==136){
            cuts_count_process_135_136++;
            if(x_bj<0.02){
                cuts_count_process_135_136_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_135_136_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_135_136_between++;
            }
        }
		 //vector meson production
        if(event_process ==91 || event_process ==92 || event_process==93 || event_process ==94){
            cuts_count_process_91_94++;
            if(x_bj<0.02){
                cuts_count_process_91_94_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_91_94_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_91_94_between++;
            }

        }
        //Resolved processes
        if(event_process ==11 || event_process ==12 || event_process==13 || event_process ==28 ||event_process ==53 || event_process ==68){
            cuts_count_process_11_13++;
            if(x_bj<0.02){
                cuts_count_process_11_13_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_11_13_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_11_13_between++;
            }

        }

		//Low Pt processes
        if(event_process ==95 ){
            cuts_count_process_95++;
            if(x_bj<0.02){
                cuts_count_process_95_shad++;
            }
            if(x_bj>0.02){
                cuts_count_process_95_non_shad++;
            }
            if(x_bj>0.002 && x_bj<0.3){
                cuts_count_process_95_between++;
            }

        }
            
          
            for(int j(0); j < nParticles; ++j ) {


                particle    = event->GetTrack(j);
                pdg         = particle->GetPdgCode();
                status      = particle->GetStatus();
                index       = particle->GetIndex();//index 1 and 2 are incoming particle electron and proton.
                pt          = particle->GetPt();
                eta         = particle->GetEta();
                phi         = particle->GetPhi();
                rap         = particle->GetRapidity();
                mass        = particle->GetM();
                theta       = particle->GetTheta(); 
                theta       = theta*1000.0; //change to mrad;
                p           = particle->GetP();
                Energy      = particle->GetE();
                xf          = particle->GetXFeynman();
                px          = particle->GetPx();
                py          = particle->GetPy();
                pz          = particle->GetPz();
                parent_index = particle->GetParentIndex();
                parent_id =particle->GetParentId();
                //event_particle = particle->GetEvent();

                p_beam.SetPxPyPzE(0.,0.,0.,sqrt(MASS_NUCLEON*MASS_NUCLEON)); //fixed target
                e_beam.SetPxPyPzE(0.,0.,pzlep,sqrt(pzlep*pzlep+MASS_ELECTRON*MASS_ELECTRON)); //initial electron
                S = (e_beam+p_beam)*(e_beam+p_beam);
            
                 
                if (pdg==-13 && status==1 && trys==0 ){

                    e_scattered.SetPxPyPzE(px,py,-pz,Energy);
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

                if (pdg==2112) {neutrons_lab++; }                      


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

                rot_z[0][0] = cos(theta_rot2);
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
              
                //----------- CMS FRAME -----------------------------------------
                if (status != 1) continue; //final state particles 

                pt_cms = sqrt(p_final.Px()*p_final.Px()+p_final.Py()*p_final.Py());//***            
                h_p_cms->Fill(p_final.P());
                h_Energy_cms->Fill(p_final.E());
                h_pz_cms->Fill(p_final.Pz());
                h_pt_cms->Fill(pt_cms);
                        

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


            //############ Rapidity cuts by charged, positive and negative particles ###############


               if (pdg == -13 || pdg == -211 || pdg==-321 || pdg==11 || pdg==-2212)
               {
                    negative_tracks++;

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
               {
                    positive_tracks++;
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
                    if (rap_cms < -1){
                        charged_nTracks_target++;
                    }                         
                
                    if (rap_cms > -0.5 && rap_cms < 0.5 ){
                        charged_nTracks_central++;
                    }
                                
                    if (rap_cms > 2){
                        charged_nTracks_projectile++;
                    }
                    
               }
               //###################################################################################

                if (pdg==2212)
                {
                    counter_p++;
                    xf_cms = 2*p_final.Pz()/sqrt(trueW2);              
                    h_nTracks_proton->Fill(counter_p);


                    //if(xf_cms > -0.2)continue;

                    cms_protons_xf_cut++;

                        if (p>0.2 && p < 0.6 ){

					counter_p_gt++;
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
                //h_rap_cms->Fill(rap_cms);
                }

                if (pdg == 211){ counter_pos_pions++;} 
                if (pdg == -211) {counter_neg_pions++;}

                
                if (p > 0.2 && p < 0.6) //Grey Tracks All Particles
                {
                    gt_counter++;
                if (pdg==211){
                    Energy_pion_gt = Energy;
                    counter_pions++;
                }
                if (pdg==321){
                    Energy_kaon_gt = Energy;
                    counter_kaons++;
                }
                if (pdg==2212 | pdg==211 | pdg==321){
                    Energy_contamination= Energy;
                    counter_contamination++; //counts grey tracks with contamination per event
                    counter_total_contamination++; //counts grey tracks with contamination for the total of events
                }
                    h_pt_gt->Fill(pt);
                    h_Energy_gt->Fill(Energy);
                            
                }

                //############### Particles without the Pt Cut ###########
 

                S_boosted = (e_beam+p_beam)*(e_beam+p_beam);
                   
                if (pdg==11 && status==1  ){

                W2_boosted = (p_beam+p_gamma)*(p_beam+p_gamma);
               
               }
       

            //Fill histograms particle loop
 
                counter++;
                counter_total_particles++;
                
                h_theta->Fill(theta);   
                h_pdg->Fill(pdg);
                h_Status->Fill(status);         
                h_eta->Fill(eta);
                
                h_rapidity->Fill(rap);       

                prof_charged_tracks_rapidity->Fill(rap,charged_tracks);
                prof_positive_tracks_rapidity->Fill(rap,positive_tracks);
                prof_negative_tracks_rapidity->Fill(rap,negative_tracks);


                
                        

            } // end of particle loop

            // if(counter_p_gt==0)
            // {
            //     sum_N_prot_gt_0 = counter_p+sum_N_prot_gt_0;

            //     // average_N_prot=sum_N_prot/counter_p;
            //     cout << "counter_p_gt: " << counter_p_gt << endl;
            //     cout << "counter_p: " << counter_p << endl;
            //     cout << "sum_N_prot : "<<sum_N_prot_gt_0 << endl;
            // }

            // cout << "sum_N_prot OUTSIDE LOOP IF: "<<sum_N_prot_gt_0 << endl;

            // cout << "average: " << average_N_prot << endl;

            // Q_hadronic_charge = positive_tracks - negative_tracks;
            // Q_hadronic_charge_backward=tracks_backward_positive - tracks_backward_negative;
            // Q_hadronic_charge_forward=tracks_forward_positive - tracks_forward_negative;

            //cout << "hadronic charge: " << Q_hadronic_charge << endl;

            // prof_Q_hadronic_charge_gt->Fill(counter_p_gt,Q_hadronic_charge);
            // prof_Q_hadronic_charge_backward_gt->Fill(counter_p_gt,Q_hadronic_charge_backward);
            // prof_Q_hadronic_charge_forward_gt->Fill(counter_p_gt,Q_hadronic_charge_forward);
                
            
            //fill histograms
            h_nTracks->Fill(nParticles);
            h_b->Fill(impact_parameter);
            h_trueQ2->Fill(trueQ2);
            
            h_nTracks_gt->Fill(gt_counter);
            h_nTracks_pos_pions->Fill(counter_pos_pions);

            h_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);
            prof_nTracks_Grey_tracks->Fill(counter_p_gt,counter_p);

            prof_nTracks_target_positive->Fill(counter_p_gt,positive_nTracks_target);
            prof_nTracks_central_positive->Fill(counter_p_gt,positive_nTracks_central);
            prof_nTracks_projectile_positive->Fill(counter_p_gt,positive_nTracks_projectile);
            prof_nTracks_target_negative->Fill(counter_p_gt,negative_nTracks_target);
            prof_nTracks_central_negative->Fill(counter_p_gt,negative_nTracks_central);
            prof_nTracks_projectile_negative->Fill(counter_p_gt,negative_nTracks_projectile);
            prof_nTracks_target_charged->Fill(counter_p_gt,charged_nTracks_target);
            prof_nTracks_central_charged->Fill(counter_p_gt,charged_nTracks_central);
            prof_nTracks_projectile_charged->Fill(counter_p_gt,charged_nTracks_projectile);
            
            h_negative_nTracks_target->Fill(negative_nTracks_target);
            h_negative_nTracks_central->Fill(negative_nTracks_central);
            h_negative_nTracks_projectile->Fill(negative_nTracks_projectile);
            h_positive_nTracks_projectile->Fill(positive_nTracks_projectile);
            h_positive_nTracks_central->Fill(positive_nTracks_central);
            h_positive_nTracks_target->Fill(positive_nTracks_target);
            h_charged_nTracks_target->Fill(charged_nTracks_target);
            h_charged_nTracks_central->Fill(charged_nTracks_central);
            h_charged_nTracks_projectile->Fill(charged_nTracks_projectile);
            h_charged_nTracks_forward->Fill(tracks_forward_charged);
            h_charged_nTracks_backward->Fill(tracks_backward_charged);

            h_positive_nTracks->Fill(positive_tracks);
            h_negative_nTracks->Fill(negative_tracks);
            h_charged_nTracks->Fill(charged_tracks);
            h_charged_nong_nTracks->Fill(charged_tracks-counter_p_gt);
            h_positive_nong_nTracks->Fill(positive_tracks-counter_p_gt);
            h_negative_nong_nTracks->Fill(negative_tracks-counter_p_gt);

            h_Energygt_gt->Fill(counter_p_gt,Energy_proton_gt);
        	prof_Energygt_gt->Fill(counter_p_gt,Energy_proton_gt);    


            h_w2->Fill(trueW2);
            h_nu->Fill(trueNu);
            h_x_bj->Fill(x_bj);

        }// end event loop

        cout << "Number of particles 										: " << counter_total_particles << endl;
        cout << "Events after event cuts 									: " << counter_event << endl;
        cout << "mean value of All protons									: " << h_nTracks_proton->GetMean() << endl;
        cout << "number scattered electrons (1 per event)					: " << counter_final_e << endl;
        cout << "number virtual photons (1 per event)						: " << p_gamma_counter << endl;
        
        cout << "number of total neutrons lab frame & final state			: " << neutrons_lab << endl;
        cout << "number struck nucleons(1 per event) 						: " << struck_proton+struck_neutron << endl;
        cout << "number struck neutron : " << struck_neutron <<	" percentage: " << (100*struck_neutron)/(struck_proton+struck_neutron)<< endl;
        cout << "number struck proton : " << struck_proton << 	" percentage: " << (100*struck_proton)/(struck_proton+struck_neutron) << endl;
        
        //cout << "number of total protons CMS frame & final state: " << counter_p << endl;
        cout << "number of total protons lab frame & final state			: " << protons_lab << endl;
        cout << "number protons xf cut										: " << cms_protons_xf_cut <<endl;
        // cout << "number of protons rap < -1: " << gt_protons_back << endl;
        // cout << "number of protons -0.5 < rap < 0.5 : " << gt_protons_zero << endl;
        // cout << "number of protons rap > 2: " << gt_protons_forw << endl;

        cout << "W2 after													: " << W2_boosted << endl;

        cout << "mean value of CHARGED PARTICLES TARGET         			: " <<h_charged_nTracks_target->GetMean() << endl;
        cout << "mean value of CHARGED PARTICLES CENTRAL        			: " <<h_charged_nTracks_central->GetMean() << endl;
        cout << "mean value of CHARGED PARTICLES PROJECTILE     			: " <<h_charged_nTracks_projectile->GetMean() << endl;

        cout << "mean value of POSITIVE PARTICLES TARGET        			: " <<h_positive_nTracks_target->GetMean() << endl;
        cout << "mean value of  POSITIVE PARTICLES CENTRAL      			: " <<h_positive_nTracks_central->GetMean() << endl;
        cout << "mean value of  POSITIVE PARTICLES PROJECTILE   			: " <<h_positive_nTracks_projectile->GetMean() << endl;

        cout << "mean value of NEGATIVE PARTICLES TARGET        			: " <<h_negative_nTracks_target->GetMean() << endl;  
        cout << "mean value of NEGATIVE PARTICLES CENTRAL       			: " <<h_negative_nTracks_central->GetMean() << endl;
        cout << "mean value of NEGATIVE PARTICLES PROJECTILE    			: " <<h_negative_nTracks_projectile->GetMean() << endl;
        
        cout << "mean value of CHARGED PARTICLES FORWARD REGION 			: " << h_charged_nTracks_forward->GetMean() << endl;
        cout << "mean value of CHARGED PARTICLES BACKWARD REGION			: " << h_charged_nTracks_backward->GetMean() << endl;

        cout << "mean value of CHARGED PARTICLES ALL REGIONS    			: " << h_charged_nTracks->GetMean() << endl;
        cout << "mean value of NEGATIVE PARTICLES ALL REGIONS   			: " << h_negative_nTracks->GetMean() << endl;
        cout << "mean value of POSITIVE PARTICLES ALL REGIONS   			: " << h_positive_nTracks->GetMean() << endl;
        cout << "mean value of CHARGED  NO NG PARTICLES ALL REGIONS   : "  << h_charged_nong_nTracks->GetMean() << endl;
        cout << "mean value of POSITIVE NO NG PARTICLES ALL REGIONS   : "  << h_positive_nong_nTracks->GetMean() << endl;
        cout << "mean value of NEGATIVE NO NG PARTICLES ALL REGIONS   : "  << h_negative_nong_nTracks->GetMean() << endl; 

        cout << "---------------------------------------------------------------------------------" << endl;
    	cout << "PROCESS MODELS INFO PYTHIA                     " << endl;
    	cout << "Total Events                           : " << counter_ALL_event << endl;
    	cout << "Total Events Process 99                : " << count_process_99 << endl;
    	cout << "Total Events Processes 131 - 132       : " << count_process_131_132 << endl;
    	cout << "Total Events Processes 135 - 136       : " << count_process_135_136 << endl;
    	cout << "Total Events Process 91-94             : " << count_process_91_94 << endl;
    	cout << "Total Events Processes Resolved        : " << count_process_11_13 << endl;
    	cout << "Total Events Processes 95              : " << count_process_95 << endl;
    	cout << "--------------------------------Region x < 0.02---------------------------------------" << endl;

    	cout << "Events Process 99 shadowing reg            : " << count_process_99_shad << endl;
    	cout << "Events Processes 131 - 132 shadowing reg   : " << count_process_131_132_shad << endl;
    	cout << "Events Processes 135 - 136 shadowing reg   : " << count_process_135_136_shad << endl;
    	cout << "Events Process 91-94 shadowing reg         : " << count_process_91_94_shad << endl;
    	cout << "Events Processes 11 - 13 shadowing reg     : " << count_process_11_13_shad << endl;
    	cout << "Events Processes 95  shadowing reg         : " << count_process_95_shad << endl;
    	cout << "-------------------------------Region x > 0.02 ---------------------------------------" << endl;
    	cout << "Events Process 99 No shadowing reg         : " << count_process_99_non_shad << endl;
    	cout << "Events Processes 131 - 132 No shadowing reg: " << count_process_131_132_non_shad << endl;
    	cout << "Events Processes 135 - 136 No shadowing reg: " << count_process_135_136_non_shad << endl;
    	cout << "Events Process 91-94 No shadowing reg      : " << count_process_91_94_non_shad << endl;
    	cout << "Events Processes 11-13 No shadowing reg    : " << count_process_11_13_non_shad << endl;
    	cout << "Events Processes 95 No shadowing reg       : " << count_process_95_non_shad << endl;
    	cout << "-------------------------------Region  x>0.002 and x<0.3 ---------------------------------" << endl;
    	cout << "Events Process 99 between              : " << count_process_99_between << endl;
    	cout << "Events Processes 131 - 132 between     : " << count_process_131_132_between << endl;
    	cout << "Events Processes 135 - 136 between     : " << count_process_135_136_between << endl;
    	cout << "Events Process 91-94 between           : " << count_process_91_94_between << endl;
    	cout << "Events Processes 11-13 between         : " << count_process_11_13_between << endl;
    	cout << "Events Processes 95 between            : " << count_process_95_between << endl;

        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "PROCESS MODELS INFO PYTHIA  after cuts	   " << endl;
        cout << "Total Events                           : " << counter_event << endl;
        cout << "Total Events Process 99                : " << cuts_count_process_99 << endl;
        cout << "Total Events Processes 131 - 132       : " << cuts_count_process_131_132 << endl;
        cout << "Total Events Processes 135 - 136       : " << cuts_count_process_135_136 << endl;
        cout << "Total Events Process 91-94             : " << cuts_count_process_91_94 << endl;
        cout << "Total Events Processes Resolved        : " << cuts_count_process_11_13 << endl;
        cout << "Total Events Processes 95              : " << cuts_count_process_95 << endl;
        cout << "--------------------------------Region x < 0.02---------------------------------------" << endl;

        cout << "Events Process 99 shadowing reg            : " << cuts_count_process_99_shad << endl;
        cout << "Events Processes 131 - 132 shadowing reg   : " << cuts_count_process_131_132_shad << endl;
        cout << "Events Processes 135 - 136 shadowing reg   : " << cuts_count_process_135_136_shad << endl;
        cout << "Events Process 91-94 shadowing reg         : " << cuts_count_process_91_94_shad << endl;
        cout << "Events Processes 11 - 13 shadowing reg     : " << cuts_count_process_11_13_shad << endl;
        cout << "Events Processes 95  shadowing reg         : " << cuts_count_process_95_shad << endl;
        cout << "-------------------------------Region x > 0.02 ---------------------------------------" << endl;
        cout << "Events Process 99 No shadowing reg         : " << cuts_count_process_99_non_shad << endl;
        cout << "Events Processes 131 - 132 No shadowing reg: " << cuts_count_process_131_132_non_shad << endl;
        cout << "Events Processes 135 - 136 No shadowing reg: " << cuts_count_process_135_136_non_shad << endl;
        cout << "Events Process 91-94 No shadowing reg      : " << cuts_count_process_91_94_non_shad << endl;
        cout << "Events Processes 11-13 No shadowing reg    : " << cuts_count_process_11_13_non_shad << endl;
        cout << "Events Processes 95 No shadowing reg       : " << cuts_count_process_95_non_shad << endl;
        cout << "-------------------------------Region  x>0.002 and x<0.3 ---------------------------------" << endl;
        cout << "Events Process 99 between              : " << cuts_count_process_99_between << endl;
        cout << "Events Processes 131 - 132 between     : " << cuts_count_process_131_132_between << endl;
        cout << "Events Processes 135 - 136 between     : " << cuts_count_process_135_136_between << endl;
        cout << "Events Process 91-94 between           : " << cuts_count_process_91_94_between << endl;
        cout << "Events Processes 11-13 between         : " << cuts_count_process_11_13_between << endl;
        cout << "Events Processes 95 between            : " << cuts_count_process_95_between << endl;


        //...............
        flux <<h_charged_nTracks_target->GetMean() << endl;
        flux <<h_charged_nTracks_central->GetMean() << endl;
        flux <<h_charged_nTracks_projectile->GetMean() << endl;

        flux <<h_positive_nTracks_target->GetMean() << endl;
        flux <<h_positive_nTracks_central->GetMean() << endl;
        flux <<h_positive_nTracks_projectile->GetMean() << endl;

        flux <<h_negative_nTracks_target->GetMean() << endl;  
        flux <<h_negative_nTracks_central->GetMean() << endl;
        flux <<h_negative_nTracks_projectile->GetMean() << endl;

        flux <<h_charged_nTracks->GetMean() << endl;
        flux <<h_positive_nTracks->GetMean() << endl;
        flux <<h_negative_nTracks->GetMean() << endl;
        flux <<h_charged_nong_nTracks->GetMean() << endl;
        flux <<h_positive_nong_nTracks->GetMean() << endl;
        flux <<h_negative_nong_nTracks->GetMean() << endl;

        //flux << h_charged_tracks->GetMean() << endl;
        //flux << h_negative_tracks->GetMean() << endl;
        //flux << h_positive_tracks->GetMean() << endl;
        //...............


        flux.close();


/*
    TCanvas *p4 = new TCanvas("p4","p4");
    p4->Divide(3,1);
    
    p4->cd(1);
    prof_nTracks_target_charged->SetTitle("Target charged;Target charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_target_charged->Scale(1.0/2.12425); //1.94
    prof_nTracks_target_charged->Scale(1.0/scale_charged_target); 
    prof_nTracks_target_charged->Draw();
    p4->cd(2);
    prof_nTracks_central_charged->SetTitle("Central charged;Central Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_charged->Scale(1.0/1.16563); //1.69
    prof_nTracks_central_charged->Scale(1.0/scale_charged_central); 
    prof_nTracks_central_charged->Draw();
    p4->cd(3);
    prof_nTracks_projectile_charged->SetTitle("Projectile charged;Projectile Region charged n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_charged->Scale(1.0/1.34466); //1.55
    prof_nTracks_projectile_charged->Scale(1.0/scale_charged_projectile); //1.55
    prof_nTracks_projectile_charged->Draw();

    TCanvas *p3 = new TCanvas("p3","p3");
    p3->Divide(3,1);
    
    p3->cd(1);
    prof_nTracks_target_negative->SetTitle("Target negative;Target - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_target_negative->Scale(1.0/0.617756); 
    prof_nTracks_target_negative->Scale(1.0/scale_negative_target); 
    prof_nTracks_target_negative->Draw();
    p3->cd(2);
    prof_nTracks_central_negative->SetTitle("Central negative;Central Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_negative->Scale(1.0/0.559972);
    prof_nTracks_central_negative->Scale(1.0/scale_negative_central);
    prof_nTracks_central_negative->Draw();
    p3->cd(3);
    prof_nTracks_projectile_negative->SetTitle("Projectile negative;Projectile Region - n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_negative->Scale(1.0/1.013);
    prof_nTracks_projectile_negative->Scale(1.0/scale_negative_projectile);
    prof_nTracks_projectile_negative->Draw();

    TCanvas *p2 = new TCanvas("p2","p2");
    p2->Divide(3,1);
    
    p2->cd(1);
    prof_nTracks_target_positive->SetTitle("Target positive;Target Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}"); 
    //prof_nTracks_target_positive->Scale(1.0/1.5065);
    prof_nTracks_target_positive->Scale(1.0/scale_positive_target);
    prof_nTracks_target_positive->Draw();
    p2->cd(2);
    prof_nTracks_central_positive->SetTitle("Central positive;Central Region + n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_central_positive->Scale(1.0/0.605654);
    prof_nTracks_central_positive->Scale(1.0/scale_positive_central);
    prof_nTracks_central_positive->Draw();
    p2->cd(3);
    prof_nTracks_projectile_positive->SetTitle("Projectile positive;Projectile + Region n_{g};R = #LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    //prof_nTracks_projectile_positive->Scale(1.0/0.331667);
    prof_nTracks_projectile_positive->Scale(1.0/scale_positive_projectile);
    prof_nTracks_projectile_positive->Draw();



*/


    
	out->Write();
	cout << "Creating output file " << endl;
        //}
       // fileIN.close();

         
  
    //out->Write();

}

