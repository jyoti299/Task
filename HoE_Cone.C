#define HoE_Cone_cxx
#include "HoE_Cone.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TProfile.h>
int main(){
    HoE_Cone a;
    a.Loop();
    return 0;
}
float DeltaR(float ecalEta, float ecalPhi, float hcalEta, float hcalPhi)
{
float dR=999.9;
  const double pi = 3.14159;

  float deltaphi = fabs(ecalPhi - hcalPhi);
  if (deltaphi > pi) deltaphi = 2. * pi - deltaphi;

  float deltaeta = fabs(ecalEta - hcalEta);
  dR = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);

  return dR;
}
void HoE_Cone::Loop()
{
//   In a ROOT session, you can do:
//      root> .L HoE_Cone.C
//      root> HoE_Cone t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

    outputFile = new TFile("test_check.root","RECREATE");

   //////////////////////////////
   //  /// cone size deltaR=0.15 ////
   //    ////////////////////////////// 
   
   TProfile*   prof_H_tot_m_max_2x2_vs_genPt_allDepths = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_allDepths","prof_H_tot_m_max_2x2_vs_genPt_allDepths",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Noise_allDepths = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Noise_allDepths","prof_H_tot_m_max_2x2_vs_genPt_Noise_allDepths",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth1 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth1","prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth1",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Depth1 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Depth1","prof_H_tot_m_max_2x2_vs_genPt_Depth1",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth2 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth2","prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth2",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Depth2 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Depth2","prof_H_tot_m_max_2x2_vs_genPt_Depth2",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth3 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth3","prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth3",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Depth3 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Depth3","prof_H_tot_m_max_2x2_vs_genPt_Depth3",40,0,200,0,2);
   TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth4 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth4","prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth4",40,0,200,0,2);
  TProfile*   prof_H_tot_m_max_2x2_vs_genPt_Depth4 = new TProfile("prof_H_tot_m_max_2x2_vs_genPt_Depth4","prof_H_tot_m_max_2x2_vs_genPt_Depth4",40,0,200,0,2);

   TH1D* h1_myHoE_conedR0p15 = new TH1D("h1_myHoE_conedR0p15", "myHoE_conedR0p15", 200, 0, 20);   
   TH1D* h1_myHoE_conedR0p15_Depth1 = new TH1D("h1_myHoE_conedR0p15_d1", "myHoE_conedR0p15_d1", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth2 = new TH1D("h1_myHoE_conedR0p15_d2", "myHoE_conedR0p15_d2", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth3 = new TH1D("h1_myHoE_conedR0p15_d3", "myHoE_conedR0p15_d3", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth4 = new TH1D("h1_myHoE_conedR0p15_d4", "myHoE_conedR0p15_d4", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_noise = new TH1D("h1_myHoE_conedR0p15_noise", "myHoE_conedR0p15_noise", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth1_noise = new TH1D("h1_myHoE_conedR0p15_d1_noise", "myHoE_conedR0p15_d1_noise", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth2_noise = new TH1D("h1_myHoE_conedR0p15_d2_noise", "myHoE_conedR0p15_d2_noise", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth3_noise = new TH1D("h1_myHoE_conedR0p15_d3_noise", "myHoE_conedR0p15_d3_noise", 200, 0, 20);
   TH1D* h1_myHoE_conedR0p15_Depth4_noise = new TH1D("h1_myHoE_conedR0p15_d4_noise", "myHoE_conedR0p15_d4_noise", 200, 0, 20);
   TH2* h2 = new TH2D("h2", "Supercluster Eta and DeltaR", 100, -2.5, 2.5, 50, 0, 5.0);
   TH2* h3 = new TH2D("h3", "Supercluster Eta and DPhi", 100, -2.5, 2.5, 50, 0, 5.0);
   TH1D* h4 = new TH1D("h4"," DPhi ",30, 0, 5);
   //TH2* h2 = new TH2D("h2", "Supercluster Eta and DeltaR",100,0,10, 50, -2.5, 2.5);
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
     cout<<" nentries "<<nentries<<endl;

      for(int iele=0; iele < pho_genmatch->size(); iele++) {
      //cout<<" photo pt size "<<pho_pt->size()<< " genmatch size "<<pho_genmatch->size()<<endl;
// declaring variables  
     
       int nhrh2=((*hcalRechitEta)[iele].size());
        float this_DeltaR = 999.9;
        float this_DeltaR_noise = 999.9;

         float hcalE_0p15=0; 
	  float hcalE_0p15_d1=0;
	  float hcalE_0p15_d2=0;
	  float hcalE_0p15_d3=0;
	  float hcalE_0p15_d4=0;
	  float hcalE_0p15_d5=0;
	  float hcalE_0p15_d6=0;
	  float hcalE_0p15_d7=0;
	  float hcalE_0p15_d8=0;
  
          float hcalE_0p15_noise=0; 
          float hcalE_0p15_d1_noise=0;
          float hcalE_0p15_d2_noise=0;
          float hcalE_0p15_d3_noise=0;
          float hcalE_0p15_d4_noise=0;
          float hcalE_0p15_d5_noise=0;
          float hcalE_0p15_d6_noise=0;
          float hcalE_0p15_d7_noise=0;
          float dPhi =0;
       for (int ihrh=0; ihrh<nhrh2; ihrh++) {
        this_DeltaR = DeltaR((*pho_sc_eta)[iele], (*pho_sc_phi)[iele], (*hcalAllRechitEta)[iele][ihrh], (*hcalAllRechitPhi)[iele][ihrh]);
        this_DeltaR_noise = DeltaR((*hcalRechitEta_noise)[iele][ihrh], (*hcalRechitPhi_noise)[iele][ihrh],(*hcalAllRechitEta)[iele][ihrh], (*hcalAllRechitPhi)[iele][ihrh]); 
           h3->Fill((*pho_sc_eta)[iele],(*hcal_diffFromPhi)[iele][ihrh]);
          // h4->Fill((*hcal_diffFromPhi)[iele][ihrh]);
            dPhi = fabs((*pho_sc_phi)[iele] - (*hcalRechitPhi_noise)[iele][ihrh]);
            h4->Fill(dPhi);
         
        int hcal_depth=(*hcalRechitDepth)[iele][ihrh];
        if ((*pho_gen_ecal)[iele] > 20) {
          if (this_DeltaR <= 0.15) {
          hcalE_0p15 = hcalE_0p15+(*hcalRechitEnergy)[iele][ihrh];
            if (hcal_depth==1) {
	      hcalE_0p15_d1=hcalE_0p15_d1+(*hcalRechitEnergy)[iele][ihrh];
	      } 	 else if (hcal_depth==2) {
              hcalE_0p15_d2=hcalE_0p15_d2+(*hcalRechitEnergy)[iele][ihrh];
              }  	 else if (hcal_depth==3) {
              hcalE_0p15_d3=hcalE_0p15_d3+(*hcalRechitEnergy)[iele][ihrh];
              }          else if (hcal_depth==4) {
              hcalE_0p15_d4=hcalE_0p15_d4+(*hcalRechitEnergy)[iele][ihrh];
              }          else if (hcal_depth==5) {
              hcalE_0p15_d5=hcalE_0p15_d5+(*hcalRechitEnergy)[iele][ihrh];
              }          else if (hcal_depth==6) {
              hcalE_0p15_d6=hcalE_0p15_d6+(*hcalRechitEnergy)[iele][ihrh];
              }          else if (hcal_depth==7) {
              hcalE_0p15_d7=hcalE_0p15_d7+(*hcalRechitEnergy)[iele][ihrh];
              }
      } //deltaR cone of 0.15
        if (this_DeltaR_noise < 0.15) {
           // h2->Fill((*pho_sc_eta)[iele],this_DeltaR_noise);
            hcalE_0p15_noise = hcalE_0p15_noise+(*hcalRechitNoise)[iele][ihrh];
         if (hcal_depth==1) {
              hcalE_0p15_d1_noise=hcalE_0p15_d1_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==2) {
              hcalE_0p15_d2_noise=hcalE_0p15_d2_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==3) {
              hcalE_0p15_d3_noise=hcalE_0p15_d3_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==4) {
              hcalE_0p15_d4_noise=hcalE_0p15_d4_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==5) {
              hcalE_0p15_d5_noise=hcalE_0p15_d5_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==6) {
              hcalE_0p15_d6_noise=hcalE_0p15_d6_noise+(*hcalRechitNoise)[iele][ihrh];
              }          else if (hcal_depth==7) {
              hcalE_0p15_d7_noise=hcalE_0p15_d7_noise+(*hcalRechitNoise)[iele][ihrh];
              }
          } //deltaR cone of 0.15 for noise       

       }// 20 GeV cut
     } //nhrh2 loop
       cout<<" ends here "<<endl; 
          float HoE_0p15 = hcalE_0p15/((*pho_sc_raw_energy)[iele]);
	  float HoE_0p15_d1 = hcalE_0p15_d1/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d2 = hcalE_0p15_d2/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d3 = hcalE_0p15_d3/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d4 = hcalE_0p15_d4/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d5 = hcalE_0p15_d5/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d6 = hcalE_0p15_d6/((*pho_sc_raw_energy)[iele]);
          float HoE_0p15_d7 = hcalE_0p15_d7/((*pho_sc_raw_energy)[iele]);
//cout<<" jentry "<<jentry <<" All "<<hcalE_0p15_noise<<" Depth1 "<<hcalE_0p15_d1_noise<<" Depth2 "<<hcalE_0p15_d2_noise<<" Depth 3 "<<hcalE_0p15_d3_noise<<" Depth 4 "<<hcalE_0p15_d4_noise<<endl;
          h1_myHoE_conedR0p15->Fill(hcalE_0p15);
          h1_myHoE_conedR0p15_Depth1->Fill(hcalE_0p15_d1);
          h1_myHoE_conedR0p15_Depth2->Fill(hcalE_0p15_d2);
          h1_myHoE_conedR0p15_Depth3->Fill(hcalE_0p15_d3);
          h1_myHoE_conedR0p15_Depth4->Fill(hcalE_0p15_d4);   
          h1_myHoE_conedR0p15_noise->Fill(hcalE_0p15_noise);
          h1_myHoE_conedR0p15_Depth1_noise->Fill(hcalE_0p15_d1_noise);
          h1_myHoE_conedR0p15_Depth2_noise->Fill(hcalE_0p15_d2_noise);
          h1_myHoE_conedR0p15_Depth3_noise->Fill(hcalE_0p15_d3_noise);
          h1_myHoE_conedR0p15_Depth4_noise->Fill(hcalE_0p15_d4_noise);

          float ECal_energy = (*pho_sc_raw_energy)[iele];//(*pho_gen_ecal)[iele];
          //h2->Fill((*pho_sc_eta)[iele],this_DeltaR_noise);
          //h2->Fill(this_DeltaR_noise,(*pho_sc_eta)[iele]);
          prof_H_tot_m_max_2x2_vs_genPt_Noise_allDepths->Fill(ECal_energy,hcalE_0p15_noise);
          prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth1->Fill(ECal_energy,hcalE_0p15_d1_noise);
          prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth2->Fill(ECal_energy,hcalE_0p15_d2_noise);
          prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth3->Fill(ECal_energy,hcalE_0p15_d3_noise);
          prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth4->Fill(ECal_energy,hcalE_0p15_d4_noise);

          prof_H_tot_m_max_2x2_vs_genPt_allDepths->Fill(ECal_energy,hcalE_0p15);
          prof_H_tot_m_max_2x2_vs_genPt_Depth1->Fill(ECal_energy,hcalE_0p15_d1);
          prof_H_tot_m_max_2x2_vs_genPt_Depth2->Fill(ECal_energy,hcalE_0p15_d2);
          prof_H_tot_m_max_2x2_vs_genPt_Depth3->Fill(ECal_energy,hcalE_0p15_d3);
          prof_H_tot_m_max_2x2_vs_genPt_Depth4->Fill(ECal_energy,hcalE_0p15_d4);

   }  //iele loop
  }

     prof_H_tot_m_max_2x2_vs_genPt_Noise_Depth1->Write();
}
