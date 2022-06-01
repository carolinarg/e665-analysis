
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "vector"
#include "TCanvas.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"


using namespace std;

int numFiles        = 5;
int numPlotsFile    = 9;
int numCanvas       = 4;
int total_plots     = 36;

void myText(Double_t x, Double_t y, const char *text, Color_t color = kBlack, Double_t tsize = 0.03, Int_t tfont = 42, Double_t langle = 0)
{

    TLatex l;
    l.SetTextAlign(12);
    l.SetTextFont(tfont);
    l.SetTextSize(tsize);
    l.SetTextAngle(langle);
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x, y, text);
}

void overlaping(){

    const vector<string> Options ={

    //"/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g0_qhat_05_1E8.root",
    //"/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g1_qhat_05_1E8.root",
    //"/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g2_qhat_05_1E8.root",
    //"/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g3_qhat_05_1E8.root",
    //"/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g0_qhat_00_1E8.root"};

    "/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g0_qhat_05_1E8.root",
    "/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g1_qhat_05_1E8.root",
    "/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g2_qhat_05_1E8.root",
    "/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g3_qhat_05_1E8.root",
    "/Volumes/Tesla/projects/output/muXe/qhat_05/sim-mars-2022/Post_Analysis_ratio_combined_Output_muXe_g0_qhat_00_1E8.root"};


    cout << "OPTIONS :" << Options[0] << endl;
    cout << "OPTIONS :" << Options[1] << endl;
    cout << "OPTIONS :" << Options[2] << endl;
    cout << "OPTIONS :" << Options[3] << endl;
    cout << "OPTIONS :" << Options[4] << endl;

    const vector<string> prof_name = {
        // target
        "prof_nTracks_target_charged",
        "prof_nTracks_target_positive",
        "prof_nTracks_target_negative"
        ,
        // central
        "prof_nTracks_central_charged",
        "prof_nTracks_central_positive",
        "prof_nTracks_central_negative"
        ,
        // projectile
        "prof_nTracks_projectile_charged",
        "prof_nTracks_projectile_positive",
        "prof_nTracks_projectile_negative",
    };



    double region_number[9] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};

    TFile *f[numFiles];

    TProfile *new_profs[45];
    TGraphAsymmErrors *new_tg[45];
    TProfile *av_distance[4];
    TH1F *distance[4];

    int k = 0;
    cout << "STARTING:" << endl;

    for (int i=0; i< numFiles;++i){
        f[i]= TFile::Open(Options[i].c_str());
        f[i]->Print();
        cout << "File Number:" << i <<endl;
        cout << "OPTIONS: " << Options[i] << endl;
           
    }

    for (int i = 0; i < numFiles; ++i) //FILES LOOP
    { 
        cout << "INSIDE FILES LOOP" << endl;
        if (i==0){
            k=0;
        };
        if (i==1){
            k=9;
        };
        if (i==2){
            k=18;
        };
        if (i==3){
            k=27;
        }
        if (i==4){
            k=36;    
        };
 

        for(int j=0;j<numPlotsFile;j++){ //plots loop

            new_profs[j + k] = (TProfile*)f[i]->Get(std::string(prof_name[j]).c_str());
            new_tg[j + k] = (TGraphAsymmErrors*)f[i]->Get(std::string("Graph;" + to_string(j + 1)).c_str());
            // if (new_tg[j + k] == nullptr) {
            //     cout<< "Problems with: " << std::string("Graph;" + to_string(j + 1)).c_str() << endl;
            //     cout << "check path directory" << endl;
            //     return -1; // better exit and fix this problem.
            // }

        }//closing plots loop
    }//closing files loop
    
    //TCanvas *c = new TCanvas("Overlaping Ratios","c",845,845);
    TCanvas *c = new TCanvas("Overlaping Ratios","c",550,550);
            
    c->DivideSquare(9,0,0);

    for (int n = 0; n < 9; ++n){
        c->cd(n+1);
        if (n>=0 && n<=2){
        new_profs[n]->SetTitle(";n_{g};R(y* < 1)");

        }
        if (n>=3 && n<=5){
        new_profs[n]->SetTitle(";n_{g};R(-0.5 < y* < 0.5)");
        }
        if (n>=6 && n<=9){
        new_profs[n]->SetTitle(";n_{g};R(y* > 2)");
        }

        gStyle->SetErrorX(0.0001);
        new_tg[n]->GetYaxis()->SetTitleOffset(3.5);
        new_tg[n]->GetXaxis()->SetTitleOffset(3.0);
        new_profs[n]->GetYaxis()->CenterTitle(true);
        new_profs[n]->GetYaxis()->SetTitleOffset(2.5);
        new_profs[n]->GetXaxis()->SetTitleOffset(2.5);
        new_tg[n]->SetMarkerColor(1);
        new_tg[n]->SetMarkerStyle(2);
        new_profs[n]->SetMarkerColor(2);
        new_profs[n]->SetLineColor(1);  //setting the color to the error bars
        new_profs[n]->SetMarkerStyle(43);
        new_profs[n+9]->SetMarkerColor(6);
        new_profs[n+9]->SetLineColor(1);
        new_profs[n+9]->SetMarkerStyle(43);
        new_profs[n+18]->SetMarkerColor(8);
        new_profs[n+18]->SetLineColor(1);
        new_profs[n+18]->SetMarkerStyle(43);
        new_profs[n+27]->SetMarkerStyle(43);
        new_profs[n+27]->SetMarkerColor(4);
        new_profs[n+27]->SetLineColor(1);
        new_profs[n+36]->SetMarkerStyle(43);
        new_profs[n+36]->SetMarkerColor(1);
        new_profs[n+36]->SetLineColor(1);
        
        if (n>0 && n<2){
        new_profs[n]->SetTitle(";n_{g};R 1");
        }
        new_profs[n]->GetXaxis()->SetLimits(-0.5,13.5);
        new_profs[n]->SetMaximum(13.5);//for y axis
        new_profs[n]->SetMinimum(-1);
            
        if(n>2 & n<6){
            new_profs[n]->SetMaximum(2.7999);
            new_profs[n]->SetMinimum(0.5);
        }
            
        if(n>5){
            new_profs[n]->SetMaximum(1.3999);
            new_profs[n]->SetMinimum(0.0);
        }
        new_profs[n]->Draw();
        new_tg[n]->Draw("PSAME");

        new_profs[n+9]->Draw("PSAME"); 
        new_profs[n+18]->Draw("PSAME");
        new_profs[n+27]->Draw("PSAME");
        new_profs[n+36]->Draw("PSAME");
         
 
        c->Update();       
                
    }
    c->cd(7);

    // For the legend on the top
    // gPad->SetTopMargin(0.1);
    // gPad->SetRightMargin(0.1);
    myText(0.25, 0.23, "Charged", kBlack, 14, 43);
    myText(0.25, 0.43, "Target", kBlack, 14, 43);
    c->cd(8);
    // gPad->SetTopMargin(0.1);
    // gPad->SetRightMargin(0.1);
    myText(0.12, 0.23, "Positive", kBlack, 14, 43);
    myText(0.25, 0.43, "Central", kBlack, 14, 43);
    c->cd(9);
    // gPad->SetTopMargin(0.1);
    // gPad->SetRightMargin(0.1);
    myText(0.12, 0.23, "Negative", kBlack, 14, 43);
    myText(0.25, 0.43, "Projectile", kBlack, 14, 43);
    
    c->cd(3);
    TLegend *lg = new TLegend(0.1, 0.350, 0.4, 0.85);
    lg->SetBorderSize(0); // no border
    lg->SetFillStyle(0);
    lg->SetFillColor(0); // Legend background should be white
    lg->SetTextFont(42);
    lg->SetTextSize(0.05); // Increase entry font size!
    lg->AddEntry(new_profs[8], "   No Gluons", "p");
    lg->AddEntry(new_profs[17], "   1 Gluon", "p");
    lg->AddEntry(new_profs[26], "   1 hard + softs Gluons", "p");
    lg->AddEntry(new_profs[35], "   Softs Gluons", "p");
    lg->AddEntry(new_profs[44], "   #hat{q} = 0.0 GeV^{2}/fm", "p");
    lg->AddEntry(new_tg[8],"   Data E665","lep");
    lg->Draw();

    myText(0.1, 0.90, "PyQM BeAGLE, #hat{q} = 0.5 GeV^{2}/fm ", kBlack, 10, 43);


    c->SaveAs("/Volumes/Tesla/projects/output/plots/sim-mars-2022/TOTAL_ratios_qhat_05_qhat_05_1E8_updated_mars_24_smaller.pdf"); 

    cout << "END LOOP" << endl;
}

