// .L rootlogon.C

// .L rootutils.C
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "TFile.h"
#include <string>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "Riostream.h"
//#include "RootUtils.C"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
const int n = 7;

void draw_final_ratios(string filename = ""){

    std::ifstream fileIN;
    // #### open txt file that include all the root files ####

    //fileIN.open("/Volumes/Tesla/projects/output/gather-plots/qhat_05/analysis_ratio.txt");
    fileIN.open("/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/analysis_ratio.txt");
    //fileIN.open("/Volumes/Tesla/projects/acceptance/data/acceptance_data.txt");
    
    while(!fileIN.eof()){
    // #### giving the path of the txt filethat include all the root files to the word filename####
        fileIN >> filename;
        cout << filename << endl;

    // #### openning root files for a specific qhat - change everytime the directory of qhat

        auto f = TFile::Open(std::string("/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/"+filename).c_str(), "READ"); //reading from runEICTree_final_copy.C

        TFile *out = new TFile (std::string("/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_"+filename).c_str(), "RECREATE");//will be read by overlaping_plots.C


    double ng_x[9][n] = {

        {0.0,1.014,1.986,3.0,4.014,5.028,6.0},
        {0.0,1.007,2.014,3.021,3.986,4.993,5.958},
        {0.0164,0.997,2.02,3.04,4.03,5.04,6.03},
        {0.0187,1.01,2.01,3.0,4.01,5.02,6.0},
        {0.0,1.01,2.0,3.01,4.02,5.01,6.0},
        {0.0,1.02,2.01,3.01,4.01,5.0,6.0},
        {0.0,1.01,1.97,2.98,4.03,4.99,6.0},
        {0.0272,0.991,2.0,3.01,4.0,5.0,6.01},
        {0.0281,1.0,2.01,3.02,4.0,5.0,6.01}

    };

    double ng_x_shifted[n] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0};

    double R_y[9][n]= {

        {1.12,2.48,4.08,5.05,6.32,6.47,7.73},
        {1.25,3.18,5.4,7.18,8.57,9.52,11.2},
        {1.01,1.52,1.84,2.02,2.68,1.65,2.41},
        {1.144,1.238,1.331,1.278,1.616,1.026,1.852},
        {1.135,1.252,1.417,1.437,1.506,0.9456,1.522},
        {1.15,1.27,1.24,1.17,1.78,1.18,2.31},
        {0.972,0.972,0.9,0.852,0.792,0.396,0.456},
        {0.971,0.933,1.05,0.965,0.447,0.373,0.419},
        {0.994,1.01,0.734,0.73,1.15,0.41,0.478}
        

    };


    double eR_y_1[9][n]={
        {0.2232, 0.0893,0.3572,0.1785,0.0,0.0446,0.2232},
        {0.0,0.0,0.0,0.0,0.1793,0.269,0.224},
        {0.0,0.0897,0.0,0.0,0.0,0.0,0.0448},
        {0.0,0.0,0.0,0.130,0.113,0.0904,0.4746},
        {0.0,0.0,0.0,0.0451,0.1241,0.6576,0.2481},
        {0.0,0.013,0.0,0.068,0.1472,0.1585,0.8377},
        {0.0,0.0,0.04469,0.05587,0.15644,0.17319,0.21788},
        {0.0,0.00559,0.0727,0.1509,0.17877,0.20112,0.40224},
        {0.0,0.0335,0.07823,0.25703,0.3463,0.37989,0.22418}

    };

    double eR_y_2[9][n]={
        {0.0,0.0,0.0,0.0,0.1786,0.0893,0.0893},
        {0.49328,0.0448,0.0,0.0,0.2691,0.2242,0.0359},
        {0.0,0.0,0.0,0.0448,0.0448,0.0,0.1794},
        {0.0,0.0,0.0,0.1930,0.0904,0.09039,0.4859},
        {0.0,0.3423,0.0,0.0564,0.1241,0.10151,0.2255},
        {0.0,0.0,0.0,0.09062,0.1471,0.15846,0.8378},
        {0.0,0.0,0.01676,0.05578,0.14525,0.15643,0.18995},
        {0.0,0.02233,0.06146,0.13407,0.21229,0.17877,0.2905},
        {0.0,0.08938,0.07821,0.134408,0.33712,0.368718,0.21229}

    };

    double eR_x_1[9][n]={
        {0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.}
    };
    double eR_x_2[9][n]={
        {0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.},{0.,0.,0.,0.,0.,0.,0.}
    };

    
  

    string namxaxis[9] = {"ng_charged_target","ng_positive_target","ng_negative_target","ng_charged_central","ng_positive_central","ng_negative_central","ng_charged_projectile","ng_positive_projectile","ng_negative_projectile"};
    string namyaxis[9] = {"R_charged_target","R_positive_target","R_negative_target","R_charged_central","R_positive_central","R_negative_central","R_charged_projectile","R_positive_projectile","R_negative_projectile"};


    auto prof_nTracks_target_positive       = (TGraph*)f->Get("prof_nTracks_target_positive");
    auto prof_nTracks_central_positive      = (TGraph*)f->Get("prof_nTracks_central_positive");
    auto prof_nTracks_projectile_positive   = (TGraph*)f->Get("prof_nTracks_projectile_positive");
    auto prof_nTracks_target_negative       = (TGraph*)f->Get("prof_nTracks_target_negative");
    auto prof_nTracks_central_negative      = (TGraph*)f->Get("prof_nTracks_central_negative");
    auto prof_nTracks_projectile_negative   = (TGraph*)f->Get("prof_nTracks_projectile_negative");
    auto prof_nTracks_target_charged        = (TGraph*)f->Get("prof_nTracks_target_charged");
    auto prof_nTracks_central_charged       = (TGraph*)f->Get("prof_nTracks_central_charged");
    auto prof_nTracks_projectile_charged    = (TGraph*)f->Get("prof_nTracks_projectile_charged");
    auto prof_charged_tracks_rapidity       = (TGraph*)f->Get("prof_charged_tracks_rapidity");
    auto prof_positive_tracks_rapidity      = (TGraph*)f->Get("prof_positive_tracks_rapidity");
    auto prof_negative_tracks_rapidity      = (TGraph*)f->Get("prof_negative_tracks_rapidity");

    auto prof_nTracks_target_positive_weighted       = (TGraph*)f->Get("prof_nTracks_target_positive_weighted");
    auto prof_nTracks_central_positive_weighted      = (TGraph*)f->Get("prof_nTracks_central_positive_weighted");
    auto prof_nTracks_projectile_positive_weighted   = (TGraph*)f->Get("prof_nTracks_projectile_positive_weighted");
    auto prof_nTracks_target_negative_weighted       = (TGraph*)f->Get("prof_nTracks_target_negative_weighted");
    auto prof_nTracks_central_negative_weighted      = (TGraph*)f->Get("prof_nTracks_central_negative_weighted");
    auto prof_nTracks_projectile_negative_weighted   = (TGraph*)f->Get("prof_nTracks_projectile_negative_weighted");
    auto prof_nTracks_target_charged_weighted        = (TGraph*)f->Get("prof_nTracks_target_charged_weighted");
    auto prof_nTracks_central_charged_weighted       = (TGraph*)f->Get("prof_nTracks_central_charged_weighted");
    auto prof_nTracks_projectile_charged_weighted    = (TGraph*)f->Get("prof_nTracks_projectile_charged_weighted");


    prof_nTracks_target_positive    ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_positive   ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_positive->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_negative    ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_negative   ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_negative->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_charged     ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_charged    ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_charged ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");

    prof_nTracks_target_positive_weighted     ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_positive_weighted    ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_positive_weighted ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_negative_weighted     ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_negative_weighted    ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_negative_weighted ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_target_charged_weighted      ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_central_charged_weighted     ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");
    prof_nTracks_projectile_charged_weighted  ->SetTitle(";n_{g};#LTn#GT_{e-Xe}/#LTn#GT_{e-D}");

    // prof_charged_tracks_rapidity    ->SetTitle(";y;#LTn#GT_{e-Xe}/#LTn#GT_{e-D} - charged");
    // prof_positive_tracks_rapidity   ->SetTitle(";y;#LTn#GT_{e-Xe}/#LTn#GT_{e-D} - positive");
    // prof_negative_tracks_rapidity   ->SetTitle(";y;#LTn#GT_{e-Xe}/#LTn#GT_{e-D} - negative");


    std::vector<TGraph*> profs = {
        // target
        prof_nTracks_target_charged,
        prof_nTracks_target_positive,
        prof_nTracks_target_negative
        ,
        // central
        prof_nTracks_central_charged,
        prof_nTracks_central_positive,
        prof_nTracks_central_negative
        ,
        // projectile
        prof_nTracks_projectile_charged,
        prof_nTracks_projectile_positive,
        prof_nTracks_projectile_negative
        ,

        //all regions
        prof_charged_tracks_rapidity,
        prof_positive_tracks_rapidity,
        prof_negative_tracks_rapidity,

        // //hadronic charge
        // prof_Q_hadronic_charge_gt,
        // prof_Q_hadronic_charge_backward_gt,
        // prof_Q_hadronic_charge_forward_gt
    
    
    };

    std::vector<TGraph*> profs_ww = {
        // target
        prof_nTracks_target_charged_weighted ,
        prof_nTracks_target_positive_weighted ,
        prof_nTracks_target_negative_weighted 
        ,
        // central
        prof_nTracks_central_charged_weighted ,
        prof_nTracks_central_positive_weighted ,
        prof_nTracks_central_negative_weighted 
        ,
        // projectile
        prof_nTracks_projectile_charged_weighted ,
        prof_nTracks_projectile_positive_weighted ,
        prof_nTracks_projectile_negative_weighted 
        
    
    };
    for (int i = 0; i < 9; ++i)
    {
        if (profs[i] == nullptr) 
        {
            cout<< "Problems with path directory" << endl;
            cout << "CHECK where are you oppenning the root files" << endl;
            return -1; // better exit and fix this problem.
        }
    }
    
  
    TGraph *tg[9];
    TGraphAsymmErrors *tg2[9];

    TCanvas *c3= new TCanvas("c3","Adams Ratio",800,800);
    c3->DivideSquare(9,0,0);
    for(int i=0;i<9;i++){
        c3->cd(i+1);
        //auto gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
        tg2[i] = new TGraphAsymmErrors(n,ng_x_shifted,R_y[i],eR_x_1[i],eR_x_2[i],eR_y_1[i],eR_y_2[i]);
        tg2[i]->GetXaxis()->SetLimits(-0.5,9.2);
        tg2[i]->SetMaximum(12.0);//for y axis
        tg2[i]->SetMinimum(-1);
        //tg2[i]->SetMarkerStyle(kOpenTriangleUp);
        
        if(i>2){
        tg2[i]->SetMaximum(2.5);
        tg2[i]->SetMinimum(0.5);
        }
        
        if(i>5)
        {
        tg2[i]->SetMaximum(1.4);
        tg2[i]->SetMinimum(0.0);
        }

        tg2[i]->Draw("AP*");
    }    
        
    //cout << "PASSED FIRST TGraph" << endl;

    TCanvas *c2= new TCanvas("c2","data Ratios",550,550);
    c2->DivideSquare(9,0,0);
    for (int i=0;i<9;i++){        
        
        c2->cd(i+1);
        tg2[i] = new TGraphAsymmErrors(n,ng_x_shifted,R_y[i],eR_x_1[i],eR_x_2[i],eR_y_1[i],eR_y_2[i]);
        //tg[i] = new TGraph(n,ng_x_shifted,R_y[i]);
            profs[i]->SetMarkerColor(1);
            profs[i]->SetLineColor(1);  //setting the color to the error bars
            profs[i]->SetMarkerStyle(22);
            tg2[i]->SetMarkerColor(1);
            tg2[i]->SetMarkerStyle(kOpenCircle);
            tg2[i]->GetXaxis()->SetLimits(-0.5,12.5);
            tg2[i]->SetMaximum(12.0);//for y axis
            tg2[i]->SetMinimum(-1);
        if(i>2){
            tg2[i]->SetMaximum(2.5);
            tg2[i]->SetMinimum(0.5);
        }
        if(i>5)
        {
            tg2[i]->SetMaximum(1.4);
            tg2[i]->SetMinimum(0.0);
        }
        //profs_ww[i]->SetMarkerStyle(kOpenCircle);
        tg2[i]->Draw("APE");
        profs[i]->Draw("PSAME");
        //profs_ww[i]->Draw("PSAME");
        tg2[i]->Write();
        profs[i]->Write();
        //profs_ww[i]->Write();

    }
    //cout << "PASSED C2" << endl;

    TLegend *lg = new TLegend(0.45, 0.85, 0.85, 0.95);
    lg->SetBorderSize(0); // no border
    lg->SetFillStyle(0);
    lg->SetFillColor(0); // Legend background should be white
    lg->SetTextFont(32);
    lg->SetTextSize(0.04); // Increase entry font size!
    lg->AddEntry(profs[3], "  No Gluons - BeAGLE", "lep");
    //lg->AddEntry(profs_ww[3], "  BeAGLE with acceptance", "lep");
    lg->AddEntry(tg2[3],"  Data E665","lep");
    
    lg->Draw();

  

    //c2->SaveAs(("/Volumes/Tesla/projects/output/plots/qhat_05/Ratio_"+filename+".pdf").c_str()); 

    

  
    out->Write();
    }
    fileIN.close();
  
}


