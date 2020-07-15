//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  7 21:57:05 2020 by ROOT version 6.14/09
// from TTree EventTree/EventData
// found on file: output_file_AllRechit_Cone.root
//////////////////////////////////////////////////////////

#ifndef HoE_Cone_h
#define HoE_Cone_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
using namespace std;
using namespace ROOT;
class HoE_Cone {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TFile *outputFile;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi_block;
   Int_t           event;
   Int_t           bunch_crossing;
   Int_t           orbit_number;
   Int_t           store_number;
   Int_t           n_pho;
   Int_t           n_hcalhit;
   Float_t         pu_true;
   Int_t           pu_obs;
   Float_t         rho;
   vector<int>     *pho_eb;
   vector<int>     *pho_ee;
   vector<int>     *pho_gap_eb_ee;
   vector<int>     *pho_gap_eb_eta;
   vector<int>     *pho_gap_eb_phi;
   vector<int>     *pho_gap_ee_dee;
   vector<int>     *pho_gap_ee_ring;
   vector<int>     *pho_golden;
   vector<int>     *pho_unknown;
   vector<int>     *pho_bigbrem;
   vector<int>     *pho_gap;
   vector<int>     *pho_badtrack;
   vector<int>     *pho_showering;
   vector<float>   *pho_track_fbrem;
   vector<float>   *pho_sc_fbrem;
   vector<float>   *pho_gen_ecal;
   vector<int>     *pho_nbrem;
   vector<int>     *pho_genmatch;
   vector<float>   *pho_dR_reco_gen;
   vector<float>   *pho_pt_ratio_reco_gen;
   vector<float>   *pho_sc_energy;
   vector<float>   *pho_sc_raw_energy;
   vector<float>   *pho_ecal_energy;
   vector<float>   *pho_seed_energy;
   vector<float>   *pho_seed_corr_energy;
   vector<float>   *pho_cmssw_hoe;
   vector<float>   *pho_cmssw_hoe_tower;
   vector<float>   *pho_cmssw_hoe_5x5;
   vector<float>   *pho_sc_eta;
   vector<float>   *pho_sc_phi;
   vector<float>   *pho_pt;
   vector<float>   *pho_eta;
   vector<float>   *pho_phi;
   vector<float>   *pho_sieie_5x5;
   vector<float>   *pho_r9_5x5;
   vector<float>   *pho_pfiso_pho;
   vector<float>   *pho_pfiso_neu;
   vector<float>   *pho_pfiso_cha;
   vector<float>   *pho_pfiso_pu;
   vector<float>   *pho_pfiso_hcal;
   vector<float>   *pho_pfiso_ecal;
   vector<float>   *pho_detiso03_ecalhit;
   vector<float>   *pho_detiso03_hcaltower1;
   vector<float>   *pho_detiso03_hcaltower2;
   vector<float>   *pho_detiso03_trk;
   vector<float>   *pho_detiso03_trk_heep;
   vector<int>     *pho_seed_detid;
   vector<int>     *pho_seed_subdetid;
   vector<int>     *pho_seed_ieta;
   vector<int>     *pho_seed_iphi;
   vector<float>   *pho_seed_eta;
   vector<float>   *pho_seed_phi;
   vector<int>     *pho_seed_raw_id;
   vector<int>     *pho_seed_hcal_ieta;
   vector<int>     *pho_seed_hcal_iphi;
   vector<int>     *hcalhit_ieta;
   vector<int>     *hcalhit_iphi;
   vector<float>   *hcalhit_energy;
   vector<int>     *hcalhit_seed_dieta;
   vector<int>     *hcalhit_seed_diphi;
   vector<int>     *hcalhit_raw_id;
   vector<int>     *hcalhit_depth;
   vector<int>     *hcalhit_pho_index;
   vector<float>   *hcalhit_eta;
   vector<float>   *hcalhit_phi;
   vector<vector<int> > *hcalRechitIeta;
   vector<vector<int> > *hcalRechitIphi;
   vector<vector<float> > *hcalRechitEnergy;
   vector<vector<int> > *hcalRechitAbsDIetaFromPhoSeed;
   vector<vector<int> > *hcalRechitAbsDIphiFromPhoSeed;
   vector<vector<int> > *hcalRechitRawID;
   vector<vector<int> > *hcalRechitDepth;
   vector<vector<float> > *hcalRechitNoise;
   vector<vector<float> > *hcalRechitEta;
   vector<float>   *diffFromPhi;
   vector<float>   *rechit_Noise;
   vector<vector<float> > *hcalRechitPhi;
   vector<vector<float> > *hcalRechitEta_noise;
   vector<vector<float> > *hcalRechitPhi_noise;
   vector<vector<float> > *hcal_diffFromPhi;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi_block;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunch_crossing;   //!
   TBranch        *b_orbit_number;   //!
   TBranch        *b_store_number;   //!
   TBranch        *b_n_pho;   //!
   TBranch        *b_n_hcalhit;   //!
   TBranch        *b_pu_true;   //!
   TBranch        *b_pu_obs;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pho_eb;   //!
   TBranch        *b_pho_ee;   //!
   TBranch        *b_pho_gap_eb_ee;   //!
   TBranch        *b_pho_gap_eb_eta;   //!
   TBranch        *b_pho_gap_eb_phi;   //!
   TBranch        *b_pho_gap_ee_dee;   //!
   TBranch        *b_pho_gap_ee_ring;   //!
   TBranch        *b_pho_golden;   //!
   TBranch        *b_pho_unknown;   //!
   TBranch        *b_pho_bigbrem;   //!
   TBranch        *b_pho_gap;   //!
   TBranch        *b_pho_badtrack;   //!
   TBranch        *b_pho_showering;   //!
   TBranch        *b_pho_track_fbrem;   //!
   TBranch        *b_pho_sc_fbrem;   //!
   TBranch        *b_pho_gen_ecal;   //!
   TBranch        *b_pho_nbrem;   //!
   TBranch        *b_pho_genmatch;   //!
   TBranch        *b_pho_dR_reco_gen;   //!
   TBranch        *b_pho_pt_ratio_reco_gen;   //!
   TBranch        *b_pho_sc_energy;   //!
   TBranch        *b_pho_sc_raw_energy;   //!
   TBranch        *b_pho_ecal_energy;   //!
   TBranch        *b_pho_seed_energy;   //!
   TBranch        *b_pho_seed_corr_energy;   //!
   TBranch        *b_pho_cmssw_hoe;   //!
   TBranch        *b_pho_cmssw_hoe_tower;   //!
   TBranch        *b_pho_cmssw_hoe_5x5;   //!
   TBranch        *b_pho_sc_eta;   //!
   TBranch        *b_pho_sc_phi;   //!
   TBranch        *b_pho_pt;   //!
   TBranch        *b_pho_eta;   //!
   TBranch        *b_pho_phi;   //!
   TBranch        *b_pho_sieie_5x5;   //!
   TBranch        *b_pho_r9_5x5;   //!
   TBranch        *b_pho_pfiso_pho;   //!
   TBranch        *b_pho_pfiso_neu;   //!
   TBranch        *b_pho_pfiso_cha;   //!
   TBranch        *b_pho_pfiso_pu;   //!
   TBranch        *b_pho_pfiso_hcal;   //!
   TBranch        *b_pho_pfiso_ecal;   //!
   TBranch        *b_pho_detiso03_ecalhit;   //!
   TBranch        *b_pho_detiso03_hcaltower1;   //!
   TBranch        *b_pho_detiso03_hcaltower2;   //!
   TBranch        *b_pho_detiso03_trk;   //!
   TBranch        *b_pho_detiso03_trk_heep;   //!
   TBranch        *b_pho_seed_detid;   //!
   TBranch        *b_pho_seed_subdetid;   //!
   TBranch        *b_pho_seed_ieta;   //!
   TBranch        *b_pho_seed_iphi;   //!
   TBranch        *b_pho_seed_eta;   //!
   TBranch        *b_pho_seed_phi;   //!
   TBranch        *b_pho_seed_raw_id;   //!
   TBranch        *b_pho_seed_hcal_ieta;   //!
   TBranch        *b_pho_seed_hcal_iphi;   //!
   TBranch        *b_hcalhit_ieta;   //!
   TBranch        *b_hcalhit_iphi;   //!
   TBranch        *b_hcalhit_energy;   //!
   TBranch        *b_hcalhit_seed_dieta;   //!
   TBranch        *b_hcalhit_seed_diphi;   //!
   TBranch        *b_hcalhit_raw_id;   //!
   TBranch        *b_hcalhit_depth;   //!
   TBranch        *b_hcalhit_pho_index;   //!
   TBranch        *b_hcalhit_eta;   //!
   TBranch        *b_hcalhit_phi;   //!
   TBranch        *b_hcalRechitIeta;   //!
   TBranch        *b_hcalRechitIphi;   //!
   TBranch        *b_hcalRechitEnergy;   //!
   TBranch        *b_hcalRechitAbsDIetaFromPhoSeed;   //!
   TBranch        *b_hcalRechitAbsDIphiFromPhoSeed;   //!
   TBranch        *b_hcalRechitRawID;   //!
   TBranch        *b_hcalRechitDepth;   //!
   TBranch        *b_hcalRechitNoise;   //!
   TBranch        *b_hcalRechitEta;   //!
   TBranch        *b_diffFromPhi;   //!
   TBranch        *b_rechit_Noise;   //!
   TBranch        *b_hcalRechitPhi;   //!
   TBranch        *b_hcalRechitEta_noise;   //!
   TBranch        *b_hcalRechitPhi_noise;   //!
   TBranch        *b_hcal_diffFromPhi;   //!

   HoE_Cone(TTree *tree=0);
   virtual ~HoE_Cone();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HoE_Cone_cxx
HoE_Cone::HoE_Cone(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_file_AllRechit_Cone.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_file_AllRechit_Cone.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("output_file_AllRechit_Cone.root:/demo");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

HoE_Cone::~HoE_Cone()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   outputFile->cd();
   outputFile->Write();
   outputFile->Close();
}

Int_t HoE_Cone::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HoE_Cone::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HoE_Cone::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pho_eb = 0;
   pho_ee = 0;
   pho_gap_eb_ee = 0;
   pho_gap_eb_eta = 0;
   pho_gap_eb_phi = 0;
   pho_gap_ee_dee = 0;
   pho_gap_ee_ring = 0;
   pho_golden = 0;
   pho_unknown = 0;
   pho_bigbrem = 0;
   pho_gap = 0;
   pho_badtrack = 0;
   pho_showering = 0;
   pho_track_fbrem = 0;
   pho_sc_fbrem = 0;
   pho_gen_ecal = 0;
   pho_nbrem = 0;
   pho_genmatch = 0;
   pho_dR_reco_gen = 0;
   pho_pt_ratio_reco_gen = 0;
   pho_sc_energy = 0;
   pho_sc_raw_energy = 0;
   pho_ecal_energy = 0;
   pho_seed_energy = 0;
   pho_seed_corr_energy = 0;
   pho_cmssw_hoe = 0;
   pho_cmssw_hoe_tower = 0;
   pho_cmssw_hoe_5x5 = 0;
   pho_sc_eta = 0;
   pho_sc_phi = 0;
   pho_pt = 0;
   pho_eta = 0;
   pho_phi = 0;
   pho_sieie_5x5 = 0;
   pho_r9_5x5 = 0;
   pho_pfiso_pho = 0;
   pho_pfiso_neu = 0;
   pho_pfiso_cha = 0;
   pho_pfiso_pu = 0;
   pho_pfiso_hcal = 0;
   pho_pfiso_ecal = 0;
   pho_detiso03_ecalhit = 0;
   pho_detiso03_hcaltower1 = 0;
   pho_detiso03_hcaltower2 = 0;
   pho_detiso03_trk = 0;
   pho_detiso03_trk_heep = 0;
   pho_seed_detid = 0;
   pho_seed_subdetid = 0;
   pho_seed_ieta = 0;
   pho_seed_iphi = 0;
   pho_seed_eta = 0;
   pho_seed_phi = 0;
   pho_seed_raw_id = 0;
   pho_seed_hcal_ieta = 0;
   pho_seed_hcal_iphi = 0;
   hcalhit_ieta = 0;
   hcalhit_iphi = 0;
   hcalhit_energy = 0;
   hcalhit_seed_dieta = 0;
   hcalhit_seed_diphi = 0;
   hcalhit_raw_id = 0;
   hcalhit_depth = 0;
   hcalhit_pho_index = 0;
   hcalhit_eta = 0;
   hcalhit_phi = 0;
   hcalRechitIeta = 0;
   hcalRechitIphi = 0;
   hcalRechitEnergy = 0;
   hcalRechitAbsDIetaFromPhoSeed = 0;
   hcalRechitAbsDIphiFromPhoSeed = 0;
   hcalRechitRawID = 0;
   hcalRechitDepth = 0;
   hcalRechitNoise = 0;
   hcalRechitEta = 0;
   diffFromPhi = 0;
   rechit_Noise = 0;
   hcalRechitPhi = 0;
   hcalRechitEta_noise = 0;
   hcalRechitPhi_noise = 0;
   hcal_diffFromPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi_block", &lumi_block, &b_lumi_block);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunch_crossing", &bunch_crossing, &b_bunch_crossing);
   fChain->SetBranchAddress("orbit_number", &orbit_number, &b_orbit_number);
   fChain->SetBranchAddress("store_number", &store_number, &b_store_number);
   fChain->SetBranchAddress("n_pho", &n_pho, &b_n_pho);
   fChain->SetBranchAddress("n_hcalhit", &n_hcalhit, &b_n_hcalhit);
   fChain->SetBranchAddress("pu_true", &pu_true, &b_pu_true);
   fChain->SetBranchAddress("pu_obs", &pu_obs, &b_pu_obs);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pho_eb", &pho_eb, &b_pho_eb);
   fChain->SetBranchAddress("pho_ee", &pho_ee, &b_pho_ee);
   fChain->SetBranchAddress("pho_gap_eb_ee", &pho_gap_eb_ee, &b_pho_gap_eb_ee);
   fChain->SetBranchAddress("pho_gap_eb_eta", &pho_gap_eb_eta, &b_pho_gap_eb_eta);
   fChain->SetBranchAddress("pho_gap_eb_phi", &pho_gap_eb_phi, &b_pho_gap_eb_phi);
   fChain->SetBranchAddress("pho_gap_ee_dee", &pho_gap_ee_dee, &b_pho_gap_ee_dee);
   fChain->SetBranchAddress("pho_gap_ee_ring", &pho_gap_ee_ring, &b_pho_gap_ee_ring);
   fChain->SetBranchAddress("pho_golden", &pho_golden, &b_pho_golden);
   fChain->SetBranchAddress("pho_unknown", &pho_unknown, &b_pho_unknown);
   fChain->SetBranchAddress("pho_bigbrem", &pho_bigbrem, &b_pho_bigbrem);
   fChain->SetBranchAddress("pho_gap", &pho_gap, &b_pho_gap);
   fChain->SetBranchAddress("pho_badtrack", &pho_badtrack, &b_pho_badtrack);
   fChain->SetBranchAddress("pho_showering", &pho_showering, &b_pho_showering);
   fChain->SetBranchAddress("pho_track_fbrem", &pho_track_fbrem, &b_pho_track_fbrem);
   fChain->SetBranchAddress("pho_sc_fbrem", &pho_sc_fbrem, &b_pho_sc_fbrem);
   fChain->SetBranchAddress("pho_gen_ecal", &pho_gen_ecal, &b_pho_gen_ecal);
   fChain->SetBranchAddress("pho_nbrem", &pho_nbrem, &b_pho_nbrem);
   fChain->SetBranchAddress("pho_genmatch", &pho_genmatch, &b_pho_genmatch);
   fChain->SetBranchAddress("pho_dR_reco_gen", &pho_dR_reco_gen, &b_pho_dR_reco_gen);
   fChain->SetBranchAddress("pho_pt_ratio_reco_gen", &pho_pt_ratio_reco_gen, &b_pho_pt_ratio_reco_gen);
   fChain->SetBranchAddress("pho_sc_energy", &pho_sc_energy, &b_pho_sc_energy);
   fChain->SetBranchAddress("pho_sc_raw_energy", &pho_sc_raw_energy, &b_pho_sc_raw_energy);
   fChain->SetBranchAddress("pho_ecal_energy", &pho_ecal_energy, &b_pho_ecal_energy);
   fChain->SetBranchAddress("pho_seed_energy", &pho_seed_energy, &b_pho_seed_energy);
   fChain->SetBranchAddress("pho_seed_corr_energy", &pho_seed_corr_energy, &b_pho_seed_corr_energy);
   fChain->SetBranchAddress("pho_cmssw_hoe", &pho_cmssw_hoe, &b_pho_cmssw_hoe);
   fChain->SetBranchAddress("pho_cmssw_hoe_tower", &pho_cmssw_hoe_tower, &b_pho_cmssw_hoe_tower);
   fChain->SetBranchAddress("pho_cmssw_hoe_5x5", &pho_cmssw_hoe_5x5, &b_pho_cmssw_hoe_5x5);
   fChain->SetBranchAddress("pho_sc_eta", &pho_sc_eta, &b_pho_sc_eta);
   fChain->SetBranchAddress("pho_sc_phi", &pho_sc_phi, &b_pho_sc_phi);
   fChain->SetBranchAddress("pho_pt", &pho_pt, &b_pho_pt);
   fChain->SetBranchAddress("pho_eta", &pho_eta, &b_pho_eta);
   fChain->SetBranchAddress("pho_phi", &pho_phi, &b_pho_phi);
   fChain->SetBranchAddress("pho_sieie_5x5", &pho_sieie_5x5, &b_pho_sieie_5x5);
   fChain->SetBranchAddress("pho_r9_5x5", &pho_r9_5x5, &b_pho_r9_5x5);
   fChain->SetBranchAddress("pho_pfiso_pho", &pho_pfiso_pho, &b_pho_pfiso_pho);
   fChain->SetBranchAddress("pho_pfiso_neu", &pho_pfiso_neu, &b_pho_pfiso_neu);
   fChain->SetBranchAddress("pho_pfiso_cha", &pho_pfiso_cha, &b_pho_pfiso_cha);
   fChain->SetBranchAddress("pho_pfiso_pu", &pho_pfiso_pu, &b_pho_pfiso_pu);
   fChain->SetBranchAddress("pho_pfiso_hcal", &pho_pfiso_hcal, &b_pho_pfiso_hcal);
   fChain->SetBranchAddress("pho_pfiso_ecal", &pho_pfiso_ecal, &b_pho_pfiso_ecal);
   fChain->SetBranchAddress("pho_detiso03_ecalhit", &pho_detiso03_ecalhit, &b_pho_detiso03_ecalhit);
   fChain->SetBranchAddress("pho_detiso03_hcaltower1", &pho_detiso03_hcaltower1, &b_pho_detiso03_hcaltower1);
   fChain->SetBranchAddress("pho_detiso03_hcaltower2", &pho_detiso03_hcaltower2, &b_pho_detiso03_hcaltower2);
   fChain->SetBranchAddress("pho_detiso03_trk", &pho_detiso03_trk, &b_pho_detiso03_trk);
   fChain->SetBranchAddress("pho_detiso03_trk_heep", &pho_detiso03_trk_heep, &b_pho_detiso03_trk_heep);
   fChain->SetBranchAddress("pho_seed_detid", &pho_seed_detid, &b_pho_seed_detid);
   fChain->SetBranchAddress("pho_seed_subdetid", &pho_seed_subdetid, &b_pho_seed_subdetid);
   fChain->SetBranchAddress("pho_seed_ieta", &pho_seed_ieta, &b_pho_seed_ieta);
   fChain->SetBranchAddress("pho_seed_iphi", &pho_seed_iphi, &b_pho_seed_iphi);
   fChain->SetBranchAddress("pho_seed_eta", &pho_seed_eta, &b_pho_seed_eta);
   fChain->SetBranchAddress("pho_seed_phi", &pho_seed_phi, &b_pho_seed_phi);
   fChain->SetBranchAddress("pho_seed_raw_id", &pho_seed_raw_id, &b_pho_seed_raw_id);
   fChain->SetBranchAddress("pho_seed_hcal_ieta", &pho_seed_hcal_ieta, &b_pho_seed_hcal_ieta);
   fChain->SetBranchAddress("pho_seed_hcal_iphi", &pho_seed_hcal_iphi, &b_pho_seed_hcal_iphi);
   fChain->SetBranchAddress("hcalhit_ieta", &hcalhit_ieta, &b_hcalhit_ieta);
   fChain->SetBranchAddress("hcalhit_iphi", &hcalhit_iphi, &b_hcalhit_iphi);
   fChain->SetBranchAddress("hcalhit_energy", &hcalhit_energy, &b_hcalhit_energy);
   fChain->SetBranchAddress("hcalhit_seed_dieta", &hcalhit_seed_dieta, &b_hcalhit_seed_dieta);
   fChain->SetBranchAddress("hcalhit_seed_diphi", &hcalhit_seed_diphi, &b_hcalhit_seed_diphi);
   fChain->SetBranchAddress("hcalhit_raw_id", &hcalhit_raw_id, &b_hcalhit_raw_id);
   fChain->SetBranchAddress("hcalhit_depth", &hcalhit_depth, &b_hcalhit_depth);
   fChain->SetBranchAddress("hcalhit_pho_index", &hcalhit_pho_index, &b_hcalhit_pho_index);
   fChain->SetBranchAddress("hcalhit_eta", &hcalhit_eta, &b_hcalhit_eta);
   fChain->SetBranchAddress("hcalhit_phi", &hcalhit_phi, &b_hcalhit_phi);
   fChain->SetBranchAddress("hcalRechitIeta", &hcalRechitIeta, &b_hcalRechitIeta);
   fChain->SetBranchAddress("hcalRechitIphi", &hcalRechitIphi, &b_hcalRechitIphi);
   fChain->SetBranchAddress("hcalRechitEnergy", &hcalRechitEnergy, &b_hcalRechitEnergy);
   fChain->SetBranchAddress("hcalRechitAbsDIetaFromPhoSeed", &hcalRechitAbsDIetaFromPhoSeed, &b_hcalRechitAbsDIetaFromPhoSeed);
   fChain->SetBranchAddress("hcalRechitAbsDIphiFromPhoSeed", &hcalRechitAbsDIphiFromPhoSeed, &b_hcalRechitAbsDIphiFromPhoSeed);
   fChain->SetBranchAddress("hcalRechitRawID", &hcalRechitRawID, &b_hcalRechitRawID);
   fChain->SetBranchAddress("hcalRechitDepth", &hcalRechitDepth, &b_hcalRechitDepth);
   fChain->SetBranchAddress("hcalRechitNoise", &hcalRechitNoise, &b_hcalRechitNoise);
   fChain->SetBranchAddress("hcalRechitEta", &hcalRechitEta, &b_hcalRechitEta);
   fChain->SetBranchAddress("diffFromPhi", &diffFromPhi, &b_diffFromPhi);
   fChain->SetBranchAddress("rechit_Noise", &rechit_Noise, &b_rechit_Noise);
   fChain->SetBranchAddress("hcalRechitPhi", &hcalRechitPhi, &b_hcalRechitPhi);
   fChain->SetBranchAddress("hcalRechitEta_noise", &hcalRechitEta_noise, &b_hcalRechitEta_noise);
   fChain->SetBranchAddress("hcalRechitPhi_noise", &hcalRechitPhi_noise, &b_hcalRechitPhi_noise);
   fChain->SetBranchAddress("hcal_diffFromPhi", &hcal_diffFromPhi, &b_hcal_diffFromPhi);
   Notify();
}

Bool_t HoE_Cone::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HoE_Cone::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HoE_Cone::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HoE_Cone_cxx
