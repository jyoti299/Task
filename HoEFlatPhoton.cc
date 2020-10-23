// -*- C++ -*-
// EleHoEAnalyzer.cc, but flattened out ie no vectors etc
// also, all hcal rechits are saved
// index pointing to ele is also saved if the min_dieta and min_diphi of hit and ele points to the same ele
// dieta and diphi kept signed

// system include files
#include <memory>
#include <cassert>

// user include files
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/plugins/CaloTopologyBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "TFile.h"
#include "TTree.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"  // added
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EGHcalRecHitSelector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
//#include "Code/HoEAnalyzer/plugins/HcalHaloAlgo.cc"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//class HoEFlatPhoton : public edm::one::EDAnalyzer<>  {
class HoEFlatPhoton : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HoEFlatPhoton(const edm::ParameterSet&);
  ~HoEFlatPhoton();
  
  static edm::ParameterSetDescription makePSetDescription();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 // TFile *file;
 // TTree *tree;
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");
  int run, lumi_block, event, bunch_crossing, orbit_number, store_number;

  int n_pho;
  std::vector<int>    pho_eb;
  std::vector<int>    pho_ee;
  std::vector<int>    pho_ecal;
  std::vector<int>    pho_gap_eb_ee;
  std::vector<int>    pho_gap_eb_eta;
  std::vector<int>    pho_gap_eb_phi;
  std::vector<int>    pho_gap_ee_dee;
  std::vector<int>    pho_gap_ee_ring;
  std::vector<float>  diffFromPhi;
  std::vector<float>  rechit_Noise;
  std::vector<int>    pho_golden;
  std::vector<int>    pho_unknown;
  std::vector<int>    pho_bigbrem;
  std::vector<int>    pho_badtrack;
  std::vector<int>    pho_showering;
  std::vector<int>    pho_gap;

  std::vector<float>  pho_track_fbrem;
  std::vector<float>  pho_sc_fbrem;
  std::vector<float>  pho_gen_ecal;
  std::vector<int>    pho_nbrem;
  std::vector<int>    pho_genmatch;
  std::vector<float>  pho_sc_energy;
  std::vector<float>  pho_dR_reco_gen;
  std::vector<float>  pho_pt_ratio_reco_gen;
  std::vector<float>  pho_sc_raw_energy;
  std::vector<float>  pho_ecal_energy;
  std::vector<float>  pho_seed_energy;
  std::vector<float>  pho_seed_corr_energy;
  std::vector<float>  pho_cmssw_hoe;
  std::vector<float>  pho_cmssw_hoe_tower;
  std::vector<float>  pho_cmssw_hoe_5x5;
  std::vector<float>  pho_sc_eta;
  std::vector<float>  pho_sc_phi;
  std::vector<float>  pho_pt;
  std::vector<float>  pho_eta;
  std::vector<float>  pho_phi;
  std::vector<float>  pho_sieie_5x5;
  std::vector<float>  pho_r9_5x5;

  std::vector<float>  pho_pfiso_pho;
  std::vector<float>  pho_pfiso_neu;
  std::vector<float>  pho_pfiso_cha;
  std::vector<float>  pho_pfiso_pu;
  std::vector<float>  pho_pfiso_hcal;
  std::vector<float>  pho_pfiso_ecal;

  std::vector<float>  pho_detiso03_ecalhit;
  std::vector<float>  pho_detiso03_hcaltower1;
  std::vector<float>  pho_detiso03_hcaltower2;
  std::vector<float>  pho_detiso03_trk;
  std::vector<float>  pho_detiso03_trk_heep;

  std::vector<int>    pho_seed_detid;
  std::vector<int>    pho_seed_subdetid;
  std::vector<int>    pho_seed_ieta;
  std::vector<int>    pho_seed_iphi;
  std::vector<float>  pho_seed_eta;
  std::vector<float>  pho_seed_phi;
  std::vector<int>    pho_seed_raw_id;
  std::vector<int>    pho_seed_hcal_ieta;
  std::vector<int>    pho_seed_hcal_iphi;

  int n_hcalhit;
  std::vector<int>    hcalhit_ieta;
  std::vector<int>    hcalhit_iphi;
  std::vector<float>  hcalhit_energy;
  std::vector<int>    hcalhit_seed_dieta;
  std::vector<int>    hcalhit_seed_diphi;
  std::vector<int>    hcalhit_raw_id;
  std::vector<int>    hcalhit_depth;
  std::vector<int>    hcalhit_pho_index;
  std::vector<float>  hcalhit_eta;
  std::vector<float>  hcalhit_phi;

  std::vector<std::vector<int>>    hcalRechitIeta;
  std::vector<std::vector<int>>    hcalRechitIphi;
  std::vector<std::vector<float>>  hcalRechitEnergy;
  std::vector<std::vector<float>>  hcalAllRechitEnergy;
  std::vector<std::vector<int>>    hcalRechitAbsDIetaFromPhoSeed;
  std::vector<std::vector<int>>    hcalRechitAbsDIphiFromPhoSeed;
  std::vector<std::vector<int>>    hcalRechitRawID;
  std::vector<std::vector<int>>    hcalRechitDepth; // mostly for Run 3 //
  std::vector<std::vector<int>>    hcalRechitDepth_noise;
  std::vector<std::vector<float>>  hcalRechitNoise; 
  std::vector<std::vector<float>>    hcalRechitEta;
  std::vector<std::vector<float>>    hcalRechitPhi;
  std::vector<std::vector<float>>    hcalRechitEta_noise;
  std::vector<std::vector<float>>    hcalRechitPhi_noise;
  std::vector<std::vector<float>>    hcalAllRechitEta;
  std::vector<std::vector<float>>    hcalAllRechitPhi;
  std::vector<std::vector<float>>    hcal_diffFromPhi;
  std::vector<int>    perPho_hcalRechitIeta;
  std::vector<int>    perPho_hcalRechitIphi;
  std::vector<float>  perPho_hcalRechitEnergy;
  std::vector<float>  perPho_hcalAllRechitEnergy;
  std::vector<int>    perPho_hcalRechitAbsDIetaFromPhoSeed;
  std::vector<int>    perPho_hcalRechitAbsDIphiFromPhoSeed;
  std::vector<int>    perPho_hcalRechitRawID;
  std::vector<int>    perPho_hcalRechitDepth; // mostly for Run 3 //
  std::vector<int>    perPho_hcalRechitDepth_noise;
  std::vector<float>  perPho_hcalRechitEta;
  std::vector<float>  perPho_hcalRechitPhi;
  std::vector<float>  perPho_hcalRechitEta_noise;
  std::vector<float>  perPho_hcalRechitPhi_noise;
  std::vector<float>  perPho_hcalAllRechitEta;
  std::vector<float>  perPho_hcalAllRechitPhi;
 
  int imin, min_dieta, min_diphi;
  float min_diR2;

  float pu_true;
  int pu_obs;
  float rho;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  static int calDIEta(int iEta1, int iEta2);
  static int calDIPhi(int iPhi1, int iPhi2);
  void reallocate_setaddress(int n_pho_ = 0, int n_hcalhit_ = 0);
  float getMinEnergyHCAL(HcalDetId id) const;
  bool EEGap,EBEtaGap,EBGap;
  int maxDIEta_ = 5;
  int maxDIPhi_ = 5;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<edm::View<reco::Photon> > phoToken_; // added
 // edm::EDGetTokenT<edm::View<pat::Photon> > phoToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> >        photonCollection_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<HBHERecHitCollection> hbhe_rechits_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;
  edm::ESHandle<CaloGeometry> theCaloGeometry;  
  edm::ESHandle<CaloTowerConstituentsMap> towerMap_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >     genParticlesCollection_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoPhoIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeuIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoChaIsoToken_; // added
  edm::Handle<std::vector<pat::PackedCandidate> > pfCands;
  std::string output;
  bool Run2_2018 ; // Now two options are supported, Run2_2018 and Run3
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HoEFlatPhoton::HoEFlatPhoton(const edm::ParameterSet& iConfig) :
 // phoToken_(consumes<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))), //added
 //phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),  
photonCollection_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonSrc"))),  
//  photonCollection_(consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonSrc"))),
  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  hbhe_rechits_(consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hbheInput"))),
  ebReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
  eeReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
  esReducedRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc"))),
  phoPhoIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("photonIsolation"))),
  phoNeuIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("neutralHadronIsolation"))),
  phoChaIsoToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("chargedIsolation"))), //added
  //output(iConfig.getParameter<std::string>("output_file")),
  Run2_2018(iConfig.getParameter<bool>("Run2_2018_"))
{
  //now do what ever initialization is needed
}


HoEFlatPhoton::~HoEFlatPhoton()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoEFlatPhoton::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
//  using mynamespace::HoEFlatPhoton;
  n_pho = 0;
  pho_eb.clear();
  pho_ee.clear();
  pho_gap_eb_ee.clear();
  pho_gap_eb_eta.clear();
  pho_gap_eb_phi.clear();
  pho_gap_ee_dee.clear();
  pho_gap_ee_ring.clear();
  diffFromPhi.clear();
  rechit_Noise.clear();
  pho_golden.clear();
  pho_unknown.clear();
  pho_badtrack.clear();
  pho_bigbrem.clear();
  pho_showering.clear();
  pho_gap.clear();
  pho_ecal.clear();
  pho_track_fbrem.clear();
  pho_sc_fbrem.clear();
  pho_nbrem.clear();
  pho_gen_ecal.clear();
  pho_genmatch.clear();
  pho_pt_ratio_reco_gen.clear();
  pho_dR_reco_gen.clear();
  pho_sc_energy.clear();
  pho_sc_raw_energy.clear();
  pho_ecal_energy.clear();
  pho_seed_energy.clear();
  pho_seed_corr_energy.clear();
  pho_cmssw_hoe.clear();
  pho_cmssw_hoe_tower.clear();
  pho_cmssw_hoe_5x5.clear();
  pho_sc_eta.clear();
  pho_sc_phi.clear();
  pho_pt.clear();
  pho_eta.clear();
  pho_phi.clear();
  pho_sieie_5x5.clear();
  pho_r9_5x5.clear();

  pho_pfiso_pho.clear();
  pho_pfiso_neu.clear();
  pho_pfiso_cha.clear();
  pho_pfiso_pu.clear();
  pho_pfiso_hcal.clear();
  pho_pfiso_ecal.clear();

  pho_detiso03_ecalhit.clear();
  pho_detiso03_hcaltower1.clear();
  pho_detiso03_hcaltower2.clear();
  pho_detiso03_trk.clear();
  pho_detiso03_trk_heep.clear();

  pho_seed_detid.clear();
  pho_seed_subdetid.clear();
  pho_seed_ieta.clear();
  pho_seed_iphi.clear();
  pho_seed_eta.clear();
  pho_seed_phi.clear();
  pho_seed_raw_id.clear();
  pho_seed_hcal_ieta.clear();
  pho_seed_hcal_iphi.clear();

  n_hcalhit = 0;
  hcalhit_ieta.clear();
  hcalhit_iphi.clear();
  hcalhit_energy.clear();
  hcalhit_seed_dieta.clear();
  hcalhit_seed_diphi.clear();
  hcalhit_raw_id.clear();
  hcalhit_depth.clear();
  hcalhit_pho_index.clear();
  hcalhit_eta.clear();
  hcalhit_phi.clear();

  hcalRechitIeta.clear();
  hcalRechitIphi.clear();
  hcalRechitEnergy.clear();
  hcalAllRechitEnergy.clear();
  hcalRechitAbsDIetaFromPhoSeed.clear();
  hcalRechitAbsDIphiFromPhoSeed.clear();
  hcalRechitRawID.clear();
  hcalRechitDepth.clear();
  hcalRechitDepth_noise.clear();
  hcalRechitNoise.clear(); 
  hcalRechitEta.clear();
  hcalRechitPhi.clear();
  hcalAllRechitEta.clear();
  hcalAllRechitPhi.clear();
  hcalRechitEta_noise.clear();
  hcalRechitPhi_noise.clear();
  hcal_diffFromPhi.clear();
  pu_true = -999999.f;
  pu_obs = -999999;
  rho = -999999.f;

  run = iEvent.eventAuxiliary().run();
  lumi_block = iEvent.eventAuxiliary().luminosityBlock();
  event = iEvent.eventAuxiliary().event();
  bunch_crossing = iEvent.eventAuxiliary().bunchCrossing();
  orbit_number = iEvent.eventAuxiliary().orbitNumber();
  store_number = iEvent.eventAuxiliary().storeNumber();
// adding from SAM
/*  std::vector<std::pair<float,float> > eleEtaPhi;
  eleEtaPhi.push_back(std::make_pair(ele->detEta(),ele->detPhi()));
  std::vector<const pat::Electron*> patEles;
  for(const auto& ele : *eleHandle){
    patEles.push_back(dynamic_cast<const pat::Electron*>(&ele));
  }
  for(size_t candNr=0;candNr<pfCands->size();candNr++){
    const pat::PackedCandidateRef pfCandRef(pfCands,candNr);
    const pat::PackedCandidate& pfCand = *pfCandRef;
    int scSeedCrysId=getSeedCrysIdOfPFCandSC(pfCandRef,patEles);
    bool accept =false;
    for(size_t eleNr=0;eleNr<eleEtaPhi.size();eleNr++){
      if(MathFuncs::calDeltaR2(eleEtaPhi[eleNr].first,eleEtaPhi[eleNr].second,
                               pfCand.eta(),pfCand.phi())<maxDR2){
        accept=true;
        break;
      }
    }
}*/
  edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);
 
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }
//  edm::Handle<edm::View<pat::Photon> > phoToken_;
//  iEvent.getByToken(photonCollection_, phoToken_);
 
  if (genPileupHandle.isValid()) {
    for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        pu_true = pu->getTrueNumInteractions();
        pu_obs = pu->getPU_NumInteractions();

        break;
      }
    }
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  if (!rhoHandle.failedToGet())
    rho = *(rhoHandle.product());
  else
    rho = -999999.f;

  edm::Handle<HBHERecHitCollection> hbheRechitsHandle;
  iEvent.getByToken(hbhe_rechits_, hbheRechitsHandle);
  iSetup.get<CaloGeometryRecord>().get(theCaloGeometry);
  iSetup.get<CaloGeometryRecord>().get(towerMap_);

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  edm::Handle<edm::ValueMap<float> > phoPhoIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoNeuIsoHandle;
  edm::Handle<edm::ValueMap<float> > phoChaIsoHandle;
  iEvent.getByToken(phoPhoIsoToken_, phoPhoIsoHandle);
  iEvent.getByToken(phoNeuIsoToken_, phoNeuIsoHandle);
  iEvent.getByToken(phoChaIsoToken_, phoChaIsoHandle);


//  if (iEvent.get(phoToken_).size() > pho_golden.capacity())
 //   reallocate_setaddress(iEvent.get(phoToken_).size(), 0);
for (edm::View<pat::Photon>::const_iterator pho = photonHandle->begin(); pho != photonHandle->end(); ++pho) {
 // for (const auto& pho : iEvent.get(phoToken_)) {
    int genmatch = 0;
    double min_dr2 = 999999.;
    double ptR = 999999.;
     
    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	//std::cout << " p->pdgId() " << p->pdgId() << std::endl;
	if ( std::abs(p->pdgId()) == 22 ) {  // ID changes from 11 to 22
	  //std::cout << "-----  p->status() " << p->status() << " p->pdgId() " << p->pdgId() << std::endl;
	}
	if ( (p->status() == 1) and (std::abs(p->pdgId()) == 22) ) {
	  double this_dr2 = reco::deltaR(*pho,*p);
          double ecal_gen = p->energy();
          pho_gen_ecal.push_back(ecal_gen);
	  if (this_dr2 < min_dr2) {
	    min_dr2 = this_dr2;
	    ptR = pho->pt() / p->pt();
	  }
	}  
      }
    }
  
    // these cuts were decided looking at min_dr and ptR distributions.
    if ( (min_dr2 < 0.0016) and (ptR > 0.7) && (ptR < 1.3) ) 
      genmatch = 1;

    perPho_hcalRechitIeta.clear();
    perPho_hcalRechitIphi.clear();
    perPho_hcalRechitEnergy.clear();
    perPho_hcalAllRechitEnergy.clear();
    perPho_hcalRechitAbsDIetaFromPhoSeed.clear();
    perPho_hcalRechitAbsDIphiFromPhoSeed.clear();
    perPho_hcalRechitRawID.clear();
    perPho_hcalRechitDepth.clear();
    perPho_hcalRechitDepth_noise.clear();
    perPho_hcalRechitEta.clear();
    perPho_hcalRechitPhi.clear();
    perPho_hcalRechitEta_noise.clear();
    perPho_hcalRechitPhi_noise.clear();
    perPho_hcalAllRechitEta.clear();
    perPho_hcalAllRechitPhi.clear();
    //rechit_Noise.clear();

    pho_dR_reco_gen.emplace_back( std::sqrt(min_dr2) );
    pho_pt_ratio_reco_gen.emplace_back(ptR);
    pho_genmatch.emplace_back(genmatch);

    pho_sc_eta.emplace_back((*pho).superCluster()->eta());
    pho_sc_phi.emplace_back((*pho).superCluster()->phi());
    pho_pt.emplace_back(pho->pt());
    pho_eta.emplace_back(pho->eta());
    pho_phi.emplace_back(pho->phi());
    pho_sieie_5x5.emplace_back(pho->full5x5_sigmaIetaIeta());
    pho_r9_5x5.emplace_back(pho->full5x5_r9());
    pho_ecal.emplace_back(pho->userFloat("ecalEnergyPreCorr")); 
   //added 
  /*  pho_pfiso_pho.push_back((*phoPhoIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
    pho_pfiso_neu.push_back((*phoNeuIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
    pho_pfiso_cha.push_back((*phoChaIsoHandle)[iEvent.get(phoToken_).ptrAt(n_pho)]);
*/
    /*pho_pfiso_pho.push_back((*phoPhoIsoHandle)[iEvent.get(photonCollection_).ptrAt(n_pho)]);
    pho_pfiso_neu.push_back((*phoNeuIsoHandle)[iEvent.get(photonCollection_).ptrAt(n_pho)]);
    pho_pfiso_cha.push_back((*phoChaIsoHandle)[iEvent.get(photonCollection_).ptrAt(n_pho)]);
*/
/*   pho_pfiso_pho       .push_back(pho->userFloat("phoChargedIsolation"));
   pho_pfiso_neu      .push_back(pho->userFloat("phoPhotonIsolation"));
   pho_pfiso_cha      .push_back(pho->userFloat("phoNeutralHadronIsolation"));
 
    reco::Photon::PflowIsolationVariables pfIso =pho.miniPFIsolation();// pho.pfIsolationVariables();
    pho_pfiso_pho.emplace_back(pfIso.sumPhotonEt);
    pho_pfiso_neu.emplace_back(pfIso.sumNeutralHadronEt);
    pho_pfiso_cha.emplace_back(pfIso.sumChargedHadronPt);
    pho_pfiso_pu.emplace_back(pfIso.sumPUPt);
    pho_pfiso_hcal.emplace_back(pho.hcalPFClusterIso());
    pho_pfiso_ecal.emplace_back(pho.ecalPFClusterIso());

    pho_detiso03_ecalhit.emplace_back(pho.dr03EcalRecHitSumEt());
    pho_detiso03_hcaltower1.emplace_back(pho.dr03HcalDepth1TowerSumEt());
    pho_detiso03_hcaltower2.emplace_back(pho.dr03HcalDepth2TowerSumEt());
    pho_detiso03_trk.emplace_back(pho.dr03TkSumPt());
    pho_detiso03_trk_heep.emplace_back(pho.dr03TkSumPtHEEP());
*/
    EcalClusterLazyTools lazyTool(iEvent, iSetup, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    //const reco::SuperCluster& superClus = *pho.superCluster();
    //const reco::CaloCluster &seedCluster = *superClus.seed();
    //DetId seedId = seedCluster.seed() ;
    DetId seedId = pho->superCluster()->seed()->seed();
    pho_seed_detid.emplace_back(seedId.det());
    pho_seed_subdetid.emplace_back(seedId.subdetId());

    float var_pho_seed_eta = -999999.f;
    float var_pho_seed_phi = -999999.f;

    DetId seed = (pho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      if ( (theSeedHit->id().rawId() != 0 ) ) {
	if (theCaloGeometry.product() != nullptr) {
	  const CaloSubdetectorGeometry *ecalgeo = theCaloGeometry->getSubdetectorGeometry(theSeedHit->id());
	  
	  if(ecalgeo->getGeometry(theSeedHit->id()) !=nullptr){
	    const GlobalPoint & ecalrechitPoint = (theCaloGeometry.product())->getPosition(theSeedHit->id());
	    var_pho_seed_eta=ecalrechitPoint.eta();
	    var_pho_seed_phi=ecalrechitPoint.phi();
	  }
	}
      }
    }  
    pho_seed_eta.emplace_back(var_pho_seed_eta);
    pho_seed_phi.emplace_back(var_pho_seed_phi);

    int var_pho_seed_ieta = -999999;
    int var_pho_seed_iphi = -999999;
    int var_pho_seed_raw_id = -999999;
    
    int var_pho_seed_hcal_ieta = -999999;
    int var_pho_seed_hcal_iphi = -999999;

    if ( seedId.det() == DetId::Ecal ) {
      if (seedId.subdetId() == EcalBarrel) {
	EBDetId ebId(seedId);
	var_pho_seed_ieta = ebId.ieta();
	var_pho_seed_iphi = ebId.iphi();
	var_pho_seed_raw_id = ebId.rawId();       
      }
      else if (seedId.subdetId() == EcalEndcap) {
	EEDetId eeId(seedId);
	var_pho_seed_ieta = eeId.ix();
        var_pho_seed_iphi = eeId.iy();
        var_pho_seed_raw_id = eeId.rawId();
      }

      // get hold of the seed hcal behind pho seed
      CaloTowerDetId towerId(towerMap_->towerOf(seedId));       
      var_pho_seed_hcal_ieta = towerId.ieta();
      var_pho_seed_hcal_iphi = towerId.iphi();
      bool hitsplushalfpi = true;
       std::vector<std::pair<float,float> > etaPhis;
      for (auto& hcalrh : iEvent.get(hbhe_rechits_) ) {
        int dIEtaAbs = std::abs(calDIEta(var_pho_seed_hcal_ieta, hcalrh.id().ieta()));
        int dIPhiAbs = std::abs(calDIPhi(var_pho_seed_hcal_iphi, hcalrh.id().iphi()));
          float rechitPhi_new=-99;
          float rechitEta_new=-99;
          float rechitEta_noise=-99;
          float rechitPhi_noise=-99;
          float diffFromSCPhi=-99;
          if ( (hcalrh.id().rawId() != 0 ) ) {
            if (theCaloGeometry.product() != nullptr) {
              const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());
              if(geo->getGeometry(hcalrh.id()) !=nullptr){
                const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());
std::cout<<" event "<<event<<std::endl;
                rechitPhi_new=rechitPoint.phi();
                rechitEta_new=rechitPoint.eta();
                diffFromSCPhi = deltaPhi(pho_sc_phi[n_pho],rechitPhi_new);
                  diffFromPhi.push_back(diffFromSCPhi);
                  std::cout<<" Doffr From SCPhi "<<diffFromSCPhi<<std::endl;
                  if( 1.25663 < fabs(diffFromSCPhi) && fabs(diffFromSCPhi) < 1.8849)
                     {   
                         rechit_Noise.push_back(hcalrh.energy()); 
                         rechitEta_noise=rechitPoint.eta();
                         rechitPhi_noise=rechitPoint.phi();
                         std::cout<<"rechit_Noise "<<hcalrh.energy()<<" "<<rechitEta_noise<<" "<<rechitPhi_noise<<std::endl; }
       std::cout<<" Check rechit_Noise "<<hcalrh.energy()<<" "<<rechitEta_noise<<" "<<rechitPhi_noise<<std::endl;
              }
            }
          }             
                        perPho_hcalRechitEta_noise.push_back(rechitEta_noise);
                        perPho_hcalRechitPhi_noise.push_back(rechitPhi_noise);
                        perPho_hcalRechitDepth_noise.push_back(hcalrh.id().depth());
                        perPho_hcalAllRechitEta.push_back(rechitEta_new);
                        perPho_hcalAllRechitPhi.push_back(rechitPhi_new);
                        perPho_hcalAllRechitEnergy.push_back(hcalrh.energy());
 /*for ( int i =0; i < 20 ;i++)
{std::cout<<"IEvent "<<i<<"Diff "<<diffFromSCPhi<<" "<<diffFromPhi[i]<<"rechit_Noise outisde"<<rechit_Noise[i]<<" "<<perPho_hcalRechitEta_noise[i]<<" "<<perPho_hcalRechitPhi_noise[i]<<std::endl;}*/
           if ( (dIEtaAbs <= maxDIEta_) && (dIPhiAbs <= maxDIPhi_)) { // &&  (hcalrh.energy()>getMinEnergyHCAL(hcalrh.id()) ) ) {
          perPho_hcalRechitIeta.push_back(hcalrh.id().ieta());
          perPho_hcalRechitIphi.push_back(hcalrh.id().iphi());
          perPho_hcalRechitEnergy.push_back(hcalrh.energy());
          perPho_hcalRechitAbsDIetaFromPhoSeed.push_back(dIEtaAbs);
          perPho_hcalRechitAbsDIphiFromPhoSeed.push_back(dIPhiAbs);

          perPho_hcalRechitRawID.push_back(hcalrh.id().rawId());
          perPho_hcalRechitDepth.push_back(hcalrh.id().depth());

          float rechitEta=-99;
          float rechitPhi=-99;
      //    float diffFromSCPhi=-99;
          if ( (hcalrh.id().rawId() != 0 ) ) {
            if (theCaloGeometry.product() != nullptr) {
              const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());
              if(geo->getGeometry(hcalrh.id()) !=nullptr){
                const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());

                rechitEta=rechitPoint.eta();
                rechitPhi=rechitPoint.phi(); 
              }
            }
          }
          perPho_hcalRechitEta.push_back(rechitEta);
          perPho_hcalRechitPhi.push_back(rechitPhi);
        }
      }
    } 


/*    pho_track_fbrem.emplace_back(pho.trackFbrem());
    pho_sc_fbrem.emplace_back(pho.superClusterFbrem());
    pho_nbrem.emplace_back(pho.numberOfBrems());

    int var_golden = 0;
    int var_unknown = 0;
    int var_gap = 0;
    int var_badtrack = 0;
    int var_showering = 0;
    int var_bigbrem = 0;

    if (pho.classification() == reco::GsfElectron::GOLDEN)
      var_golden = 1;

    if (pho.classification() == reco::GsfElectron::UNKNOWN)
      var_unknown = 1;

    if (pho.classification() == reco::GsfElectron::BIGBREM)
      var_bigbrem = 1;

    if (pho.classification() == reco::GsfElectron::BADTRACK)
      var_badtrack = 1;

    if (pho.classification() == reco::GsfElectron::SHOWERING)
      var_showering = 1;

    if (pho.classification() == reco::GsfElectron::GAP)
      var_gap = 1;

    pho_golden.emplace_back(var_golden);
    pho_unknown.emplace_back(var_unknown);
    pho_gap.emplace_back(var_gap);
    pho_badtrack.emplace_back(var_badtrack);
    pho_showering.emplace_back(var_showering);
    pho_bigbrem.emplace_back(var_bigbrem);
*/
    pho_sc_energy.emplace_back((*pho).superCluster()->energy());
    pho_sc_raw_energy.emplace_back((*pho).superCluster()->rawEnergy());
    pho_seed_energy.emplace_back((*theSeedHit).energy());
    //pho_seed_corr_energy.emplace_back(seedCluster->correctedEnergy());
    //pho_ecal_energy.emplace_back(pho.ecalEnergy());
    pho_cmssw_hoe.emplace_back(pho->hadTowOverEm());    //added HoverE diffrnt for Photon and commneyed above line
    //pho_cmssw_hoe_tower.emplace_back(pho.hcalOverEcalBc());
    //pho_cmssw_hoe_5x5.emplace_back(pho.full5x5_hcalOverEcal());
   
    pho_seed_ieta.emplace_back(var_pho_seed_ieta);
    pho_seed_iphi.emplace_back(var_pho_seed_iphi);
    pho_seed_raw_id.emplace_back(var_pho_seed_raw_id);

    pho_seed_hcal_ieta.emplace_back(var_pho_seed_hcal_ieta);
    pho_seed_hcal_iphi.emplace_back(var_pho_seed_hcal_iphi);

    pho_eb.emplace_back(pho->isEB());
    pho_ee.emplace_back(pho->isEE());
    pho_gap_eb_ee.emplace_back(pho->isEBEEGap());
    pho_gap_eb_eta.emplace_back(pho->isEBEtaGap());
    pho_gap_eb_phi.emplace_back(pho->isEBPhiGap());
    pho_gap_ee_dee.emplace_back(pho->isEEDeeGap());
    pho_gap_ee_ring.emplace_back(pho->isEERingGap());
    //diffFromPhi.emplace_back(diffFromSCPhi);
    //// Added by me by making this a vector of vectors 
    hcalRechitIeta.push_back(perPho_hcalRechitIeta);
    hcalRechitIphi.push_back(perPho_hcalRechitIphi);
    hcalRechitEnergy.push_back(perPho_hcalRechitEnergy);
    hcalAllRechitEnergy.push_back(perPho_hcalAllRechitEnergy);
    hcalRechitAbsDIetaFromPhoSeed.push_back(perPho_hcalRechitAbsDIetaFromPhoSeed);
    hcalRechitAbsDIphiFromPhoSeed.push_back(perPho_hcalRechitAbsDIphiFromPhoSeed);
    hcalRechitRawID.push_back(perPho_hcalRechitRawID);
    hcalRechitDepth.push_back(perPho_hcalRechitDepth);
    hcalRechitDepth_noise.push_back(perPho_hcalRechitDepth_noise);
    hcalRechitNoise.push_back(rechit_Noise);
    hcalRechitEta.push_back(perPho_hcalRechitEta);
    hcalRechitPhi.push_back(perPho_hcalRechitPhi);
    hcalAllRechitEta.push_back(perPho_hcalAllRechitEta);
    hcalAllRechitPhi.push_back(perPho_hcalAllRechitPhi);
    hcalRechitEta_noise.push_back(perPho_hcalRechitEta_noise);
    hcalRechitPhi_noise.push_back(perPho_hcalRechitPhi_noise);
    hcal_diffFromPhi.push_back(diffFromPhi);
    ++n_pho;
    //std::cout<<" No of photons after adding"<<n_pho<<std::endl;
  }

  // given the context, should be ok...
  if (n_pho == 0)
    return;
/*
  // just in case
  assert(((void) "ERROR: pho_eb size doesn't match n_pho!!!", int(pho_eb.size()) == n_pho));
  assert(((void) "ERROR: pho_ee size doesn't match n_pho!!!", int(pho_ee.size()) == n_pho));
  assert(((void) "ERROR: pho_gap_eb_ee size doesn't match n_pho!!!", int(pho_gap_eb_ee.size()) == n_pho));
  assert(((void) "ERROR: pho_gab_eb_eta size doesn't match n_pho!!!", int(pho_gap_eb_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_gap_eb_phi size doesn't match n_pho!!!", int(pho_gap_eb_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_gap_ee_dee size doesn't match n_pho!!!", int(pho_gap_ee_dee.size()) == n_pho));
  assert(((void) "ERROR: pho_gap_ee_ring size doesn't match n_pho!!!", int(pho_gap_ee_ring.size()) == n_pho));
*/
  /*assert(((void) "ERROR: pho_golden size doesn't match n_pho!!!", int(pho_golden.size()) == n_pho));
  assert(((void) "ERROR: pho_unknown size doesn't match n_pho!!!", int(pho_unknown.size()) == n_pho));
  assert(((void) "ERROR: pho_badtrack size doesn't match n_pho!!!", int(pho_badtrack.size()) == n_pho));
  assert(((void) "ERROR: pho_bigbrem size doesn't match n_pho!!!", int(pho_bigbrem.size()) == n_pho));
  assert(((void) "ERROR: pho_showering size doesn't match n_pho!!!", int(pho_showering.size()) == n_pho));
  // assert(((void) "ERROR: pho_gap size doesn't match n_pho!!!", int(pho_gap.size()) == n_pho));

  assert(((void) "ERROR: pho_track_fbrem size doesn't match n_pho!!!", int(pho_track_fbrem.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_fbrem size doesn't match n_pho!!!", int(pho_sc_fbrem.size()) == n_pho));
  assert(((void) "ERROR: pho_nbrem size doesn't match n_pho!!!", int(pho_nbrem.size()) == n_pho));*/
/*
  assert(((void) "ERROR: pho_genmatch size doesn't match n_pho!!!", int(pho_genmatch.size()) == n_pho));
  assert(((void) "ERROR: pho_pt_ratio_reco_gen size doesn't match n_pho!!!", int(pho_pt_ratio_reco_gen.size()) == n_pho));
  assert(((void) "ERROR: pho_dR_reco_gen size doesn't match n_pho!!!", int(pho_dR_reco_gen.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_energy size doesn't match n_pho!!!", int(pho_sc_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_raw_energy size doesn't match n_pho!!!", int(pho_sc_raw_energy.size()) == n_pho));
  //assert(((void) "ERROR: pho_ecal_energy size doesn't match n_pho!!!", int(pho_ecal_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_energy size doesn't match n_pho!!!", int(pho_seed_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_corr_energy size doesn't match n_pho!!!", int(pho_seed_corr_energy.size()) == n_pho));
  assert(((void) "ERROR: pho_cmssw_hoe size doesn't match n_pho!!!", int(pho_cmssw_hoe.size()) == n_pho));
 // assert(((void) "ERROR: pho_cmssw_hoe_tower size doesn't match n_pho!!!", int(pho_cmssw_hoe_tower.size()) == n_pho));
 // assert(((void) "ERROR: pho_cmssw_hoe_5x5 size doesn't match n_pho!!!", int(pho_cmssw_hoe_5x5.size()) == n_pho));*/
/*  assert(((void) "ERROR: pho_sc_eta size doesn't match n_pho!!!", int(pho_sc_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_sc_phi size doesn't match n_pho!!!", int(pho_sc_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_pt size doesn't match n_pho!!!", int(pho_pt.size()) == n_pho));
  assert(((void) "ERROR: pho_eta size doesn't match n_pho!!!", int(pho_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_phi size doesn't match n_pho!!!", int(pho_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_sieie_5x5 size doesn't match n_pho!!!", int(pho_sieie_5x5.size()) == n_pho));
  assert(((void) "ERROR: pho_r9_5x5 size doesn't match n_pho!!!", int(pho_r9_5x5.size()) == n_pho));

  assert(((void) "ERROR: pho_pfiso_pho size doesn't match n_pho!!!", int(pho_pfiso_pho.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_neu size doesn't match n_pho!!!", int(pho_pfiso_neu.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_cha size doesn't match n_pho!!!", int(pho_pfiso_cha.size()) == n_pho));
*/
 /* assert(((void) "ERROR: pho_pfiso_pu size doesn't match n_pho!!!", int(pho_pfiso_pu.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_hcal size doesn't match n_pho!!!", int(pho_pfiso_hcal.size()) == n_pho));
  assert(((void) "ERROR: pho_pfiso_ecal size doesn't match n_pho!!!", int(pho_pfiso_ecal.size()) == n_pho));

  assert(((void) "ERROR: pho_detiso03_ecalhit size doesn't match n_pho!!!", int(pho_detiso03_ecalhit.size()) == n_pho));
  assert(((void) "ERROR: pho_detiso03_hcaltower1 size doesn't match n_pho!!!", int(pho_detiso03_hcaltower1.size()) == n_pho));
  assert(((void) "ERROR: pho_detiso03_hcaltower2 size doesn't match n_pho!!!", int(pho_detiso03_hcaltower2.size()) == n_pho));
  assert(((void) "ERROR: pho_detiso03_trk size doesn't match n_pho!!!", int(pho_detiso03_trk.size()) == n_pho));
  assert(((void) "ERROR: pho_detiso03_trk_heep size doesn't match n_pho!!!", int(pho_detiso03_trk_heep.size()) == n_pho));
*/
 /* assert(((void) "ERROR: pho_seed_detid size doesn't match n_pho!!!", int(pho_seed_detid.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_subdetid size doesn't match n_pho!!!", int(pho_seed_subdetid.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_ieta size doesn't match n_pho!!!", int(pho_seed_ieta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_iphi size doesn't match n_pho!!!", int(pho_seed_iphi.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_eta size doesn't match n_pho!!!", int(pho_seed_eta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_phi size doesn't match n_pho!!!", int(pho_seed_phi.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_raw_id size doesn't match n_pho!!!", int(pho_seed_raw_id.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_hcal_ieta size doesn't match n_pho!!!", int(pho_seed_hcal_ieta.size()) == n_pho));
  assert(((void) "ERROR: pho_seed_hcal_iphi size doesn't match n_pho!!!", int(pho_seed_hcal_iphi.size()) == n_pho));
*/
  if (iEvent.get(hbhe_rechits_).size() > hcalhit_depth.capacity())  
   reallocate_setaddress(0, iEvent.get(hbhe_rechits_).size());
   
  std::cout<<" HCAL Depth Capacity "<<hcalhit_depth.capacity()<<std::endl;
  for (auto& hcalrh : iEvent.get(hbhe_rechits_)) {
   // if (hcalrh.energy() < getMinEnergyHCAL(hcalrh.id()))
   //   continue;
    if (hcalrh.id().depth() < 0 or hcalrh.id().depth() > 100)
      continue;

    hcalhit_ieta.emplace_back(hcalrh.id().ieta());
    hcalhit_iphi.emplace_back(hcalrh.id().iphi());
    hcalhit_energy.emplace_back(hcalrh.energy());

    hcalhit_raw_id.emplace_back(hcalrh.id().rawId());
    hcalhit_depth.emplace_back(hcalrh.id().depth());

    float rechitEta = -999999.f;
    float rechitPhi = -999999.f;
    if (hcalrh.id().rawId() != 0) {
      if (theCaloGeometry.product() != nullptr) {
        const CaloSubdetectorGeometry *geo = theCaloGeometry->getSubdetectorGeometry(hcalrh.id());

        if (geo->getGeometry(hcalrh.id()) != nullptr) {
          const GlobalPoint & rechitPoint = theCaloGeometry.product()->getPosition(hcalrh.id());

          rechitEta=rechitPoint.eta();
          rechitPhi=rechitPoint.phi();	
        }
      }
    }
    hcalhit_eta.emplace_back(rechitEta);
    hcalhit_phi.emplace_back(rechitPhi);

    imin = -1;
    min_dieta = 999999;
    min_diphi = 999999;
    min_diR2 = 999999.f;

    for (int iE = 0; iE < n_pho; ++iE) {
      int dieta = calDIEta(pho_seed_hcal_ieta[iE], hcalhit_ieta.back());
      int diphi = calDIPhi(pho_seed_hcal_iphi[iE], hcalhit_iphi.back());

      if (float(dieta * dieta) + float(diphi * diphi) < min_diR2) {
        min_diR2 = float(dieta * dieta) + float(diphi * diphi);
        min_dieta = dieta;
        min_diphi = diphi;
        imin = iE;
      }
    }

    hcalhit_seed_dieta.emplace_back(min_dieta);
    hcalhit_seed_diphi.emplace_back(min_diphi);
    hcalhit_pho_index.emplace_back(imin);

    ++n_hcalhit;
  }

  // again, just in case
  assert(((void) "ERROR: hcalhit_ieta size doesn't match n_hcalhit!!!", int(hcalhit_ieta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_iphi size doesn't match n_hcalhit!!!", int(hcalhit_iphi.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_energy size doesn't match n_hcalhit!!!", int(hcalhit_energy.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_seed_dieta size doesn't match n_hcalhit!!!", int(hcalhit_seed_dieta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_seed_diphi size doesn't match n_hcalhit!!!", int(hcalhit_seed_diphi.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_raw_id size doesn't match n_hcalhit!!!", int(hcalhit_raw_id.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_depth size doesn't match n_hcalhit!!!", int(hcalhit_depth.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_pho_index size doesn't match n_hcalhit!!!", int(hcalhit_pho_index.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_eta size doesn't match n_hcalhit!!!", int(hcalhit_eta.size()) == n_hcalhit));
  assert(((void) "ERROR: hcalhit_phi size doesn't match n_hcalhit!!!", int(hcalhit_phi.size()) == n_hcalhit));

  tree->Fill();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// doing some blatant copy paste from RecoEgamma/EgammaIsolationAlgos/src/EGHcalRecHitSphoctor.cc
int HoEFlatPhoton::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

int HoEFlatPhoton::calDIEta(int iEta1, int iEta2) {
  int dEta = iEta1 - iEta2;
  if (iEta1 * iEta2 < 0) {  //-ve to +ve transition and no crystal at zero
    if (dEta < 0)
      ++dEta;
    else
      --dEta;
  }

  return dEta;
}

// HCAL thresholds from here https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/CaloTowersCreator/python/calotowermaker_cfi.py?%21v=CMSSW_10_6_2
// Note: As far as I understood, 
// for 2018, HB threshold is 0.7, and for Run 3 it becomes 0.1 in depth1, 0.2 in depth2, 0.3 in other depths.
// In HE, 2018 and Run3 is same, and it is 0.1 in depth1, and 0.2 in other depths.
// Double check these HCAL thresholds from Sam.
float HoEFlatPhoton::getMinEnergyHCAL(HcalDetId id) const {
  if ( (id.subdetId() == HcalBarrel)  ) {
    if ( (Run2_2018 == 1) )
      return 0.7f;
    else if ( (Run2_2018 == 0) ) { // means Run3
      if (id.depth() == 1)
	return 0.1f;
      else if (id.depth() == 2)
	return 0.2f;
      else
	return 0.3f;
    }
    else // neither 2018, nor Run3, not supported
      return 999999.f;
  } 
  else if (id.subdetId() == HcalEndcap) {
    if (id.depth() == 1)
      return 0.1f;
    else
      return 0.2f;
  } else
    return 999999.f;
}

void HoEFlatPhoton::reallocate_setaddress(int n_pho_, int n_hcalhit_)
{
  static int cap_pho = 8;
  cap_pho = (n_pho_ == 0) ? cap_pho : n_pho_;

  pho_eb.reserve(cap_pho);
  pho_ee.reserve(cap_pho);
  pho_gap_eb_ee.reserve(cap_pho);
  pho_gap_eb_eta.reserve(cap_pho);
  pho_gap_eb_phi.reserve(cap_pho);
  pho_gap_ee_dee.reserve(cap_pho);
  pho_gap_ee_ring.reserve(cap_pho);

  pho_golden.reserve(cap_pho);
  pho_unknown.reserve(cap_pho);
  pho_badtrack.reserve(cap_pho);
  pho_bigbrem.reserve(cap_pho);
  pho_showering.reserve(cap_pho);
  pho_gap.reserve(cap_pho);

  pho_track_fbrem.reserve(cap_pho);
  pho_sc_fbrem.reserve(cap_pho);
  pho_nbrem.reserve(cap_pho);

  pho_genmatch.reserve(cap_pho);
  pho_pt_ratio_reco_gen.reserve(cap_pho);
  pho_dR_reco_gen.reserve(cap_pho);
  pho_sc_energy.reserve(cap_pho);
  pho_sc_raw_energy.reserve(cap_pho);
  pho_ecal_energy.reserve(cap_pho);
  pho_seed_energy.reserve(cap_pho);
  pho_seed_corr_energy.reserve(cap_pho);
  pho_cmssw_hoe.reserve(cap_pho);
  pho_cmssw_hoe_tower.reserve(cap_pho);
  pho_cmssw_hoe_5x5.reserve(cap_pho);
  pho_sc_eta.reserve(cap_pho);
  pho_sc_phi.reserve(cap_pho);
  pho_pt.reserve(cap_pho);
  pho_eta.reserve(cap_pho);
  pho_phi.reserve(cap_pho);
  pho_sieie_5x5.reserve(cap_pho);
  pho_r9_5x5.reserve(cap_pho);
  pho_ecal.reserve(cap_pho);

  pho_pfiso_pho.reserve(cap_pho);
  pho_pfiso_neu.reserve(cap_pho);
  pho_pfiso_cha.reserve(cap_pho);
  pho_pfiso_pu.reserve(cap_pho);
  pho_pfiso_hcal.reserve(cap_pho);
  pho_pfiso_ecal.reserve(cap_pho);

  pho_detiso03_ecalhit.reserve(cap_pho);
  pho_detiso03_hcaltower1.reserve(cap_pho);
  pho_detiso03_hcaltower2.reserve(cap_pho);
  pho_detiso03_trk.reserve(cap_pho);
  pho_detiso03_trk_heep.reserve(cap_pho);

  pho_seed_detid.reserve(cap_pho);
  pho_seed_subdetid.reserve(cap_pho);
  pho_seed_ieta.reserve(cap_pho);
  pho_seed_iphi.reserve(cap_pho);
  pho_seed_eta.reserve(cap_pho);
  pho_seed_phi.reserve(cap_pho);
  pho_seed_raw_id.reserve(cap_pho);
  pho_seed_hcal_ieta.reserve(cap_pho);
  pho_seed_hcal_iphi.reserve(cap_pho);

  static int cap_hcalhit = 10000;
  cap_hcalhit = (n_hcalhit_ == 0) ? cap_hcalhit : n_hcalhit_;
  hcalhit_ieta.reserve(cap_hcalhit);
  hcalhit_iphi.reserve(cap_hcalhit);
  hcalhit_energy.reserve(cap_hcalhit);
  hcalhit_seed_dieta.reserve(cap_hcalhit);
  hcalhit_seed_diphi.reserve(cap_hcalhit);
  hcalhit_raw_id.reserve(cap_hcalhit);
  hcalhit_depth.reserve(cap_hcalhit);
  hcalhit_pho_index.reserve(cap_hcalhit);
  hcalhit_eta.reserve(cap_hcalhit);
  hcalhit_phi.reserve(cap_hcalhit);

  if (n_pho_ == 0 and n_hcalhit_ == 0) {
    tree->Branch("run", &run, "run/I");
    tree->Branch("lumi_block", &lumi_block, "lumi_block/I");
    tree->Branch("event", &event, "event/I");
    tree->Branch("bunch_crossing", &bunch_crossing, "bunch_crossing/I");
    tree->Branch("orbit_number", &orbit_number, "orbit_number/I");
    tree->Branch("store_number", &store_number, "store_number/I");

    tree->Branch("n_pho", &n_pho, "n_pho/I");
    tree->Branch("n_hcalhit", &n_hcalhit, "n_hcalhit/I");
    tree->Branch("pu_true", &pu_true, "pu_true/F");
    tree->Branch("pu_obs", &pu_obs, "pu_obs/I");
    tree->Branch("rho", &rho, "rho/F");
  }

  TBranch *b_pho_eb = tree->Branch("pho_eb", &pho_eb);
  TBranch *b_pho_ee = tree->Branch("pho_ee", &pho_ee);
  TBranch *b_pho_gap_eb_ee = tree->Branch("pho_gap_eb_ee", &pho_gap_eb_ee);
  TBranch *b_pho_gap_eb_eta = tree->Branch("pho_gap_eb_eta", &pho_gap_eb_eta);
  TBranch *b_pho_gap_eb_phi = tree->Branch("pho_gap_eb_phi", &pho_gap_eb_phi);
  TBranch *b_pho_gap_ee_dee = tree->Branch("pho_gap_ee_dee", &pho_gap_ee_dee);
  TBranch *b_pho_gap_ee_ring = tree->Branch("pho_gap_ee_ring", &pho_gap_ee_ring);
  TBranch *b_pho_ecal = tree->Branch("pho_ecal", &pho_ecal);

  TBranch *b_pho_golden = tree->Branch("pho_golden", &pho_golden);
  TBranch *b_pho_unknown = tree->Branch("pho_unknown", &pho_unknown);
  TBranch *b_pho_bigbrem = tree->Branch("pho_bigbrem", &pho_bigbrem);
  TBranch *b_pho_gap = tree->Branch("pho_gap", &pho_gap);
  TBranch *b_pho_badtrack = tree->Branch("pho_badtrack", &pho_badtrack);
  TBranch *b_pho_showering = tree->Branch("pho_showering", &pho_showering);

  TBranch *b_pho_track_fbrem = tree->Branch("pho_track_fbrem", &pho_track_fbrem);
  TBranch *b_pho_sc_fbrem = tree->Branch("pho_sc_fbrem", &pho_sc_fbrem);
  TBranch *b_pho_gen_ecal = tree->Branch("pho_gen_ecal",&pho_gen_ecal);
  TBranch *b_pho_nbrem = tree->Branch("pho_nbrem", &pho_nbrem);
  TBranch *b_pho_genmatch = tree->Branch("pho_genmatch", &pho_genmatch);
  TBranch *b_pho_dR_reco_gen = tree->Branch("pho_dR_reco_gen", &pho_dR_reco_gen);
  TBranch *b_pho_pt_ratio_reco_gen = tree->Branch("pho_pt_ratio_reco_gen", &pho_pt_ratio_reco_gen);
  TBranch *b_pho_sc_energy = tree->Branch("pho_sc_energy", &pho_sc_energy);
  TBranch *b_pho_sc_raw_energy = tree->Branch("pho_sc_raw_energy", &pho_sc_raw_energy);
  TBranch *b_pho_ecal_energy = tree->Branch("pho_ecal_energy", &pho_ecal_energy);
  TBranch *b_pho_seed_energy = tree->Branch("pho_seed_energy", &pho_seed_energy);
  TBranch *b_pho_seed_corr_energy = tree->Branch("pho_seed_corr_energy", &pho_seed_corr_energy);
  TBranch *b_pho_cmssw_hoe = tree->Branch("pho_cmssw_hoe", &pho_cmssw_hoe);
  TBranch *b_pho_cmssw_hoe_tower = tree->Branch("pho_cmssw_hoe_tower", &pho_cmssw_hoe_tower);
  TBranch *b_pho_cmssw_hoe_5x5 = tree->Branch("pho_cmssw_hoe_5x5", &pho_cmssw_hoe_5x5);
  TBranch *b_pho_sc_eta = tree->Branch("pho_sc_eta", &pho_sc_eta);
  TBranch *b_pho_sc_phi = tree->Branch("pho_sc_phi", &pho_sc_phi);
  TBranch *b_pho_pt = tree->Branch("pho_pt", &pho_pt);
  TBranch *b_pho_eta = tree->Branch("pho_eta", &pho_eta);
  TBranch *b_pho_phi = tree->Branch("pho_phi", &pho_phi);
  TBranch *b_pho_sieie_5x5 = tree->Branch("pho_sieie_5x5", &pho_sieie_5x5);
  TBranch *b_pho_r9_5x5 = tree->Branch("pho_r9_5x5", &pho_r9_5x5);

  TBranch *b_pho_pfiso_pho = tree->Branch("pho_pfiso_pho", &pho_pfiso_pho);
  TBranch *b_pho_pfiso_neu = tree->Branch("pho_pfiso_neu", &pho_pfiso_neu);
  TBranch *b_pho_pfiso_cha = tree->Branch("pho_pfiso_cha", &pho_pfiso_cha);
  TBranch *b_pho_pfiso_pu = tree->Branch("pho_pfiso_pu", &pho_pfiso_pu);
  TBranch *b_pho_pfiso_hcal = tree->Branch("pho_pfiso_hcal", &pho_pfiso_hcal);
  TBranch *b_pho_pfiso_ecal = tree->Branch("pho_pfiso_ecal", &pho_pfiso_ecal);

  TBranch *b_pho_detiso03_ecalhit = tree->Branch("pho_detiso03_ecalhit", &pho_detiso03_ecalhit);
  TBranch *b_pho_detiso03_hcaltower1 = tree->Branch("pho_detiso03_hcaltower1", &pho_detiso03_hcaltower1);
  TBranch *b_pho_detiso03_hcaltower2 = tree->Branch("pho_detiso03_hcaltower2", &pho_detiso03_hcaltower2);
  TBranch *b_pho_detiso03_trk = tree->Branch("pho_detiso03_trk", &pho_detiso03_trk);
  TBranch *b_pho_detiso03_trk_heep = tree->Branch("pho_detiso03_trk_heep", &pho_detiso03_trk_heep);

  TBranch *b_pho_seed_detid = tree->Branch("pho_seed_detid", &pho_seed_detid);
  TBranch *b_pho_seed_subdetid = tree->Branch("pho_seed_subdetid", &pho_seed_subdetid);
  TBranch *b_pho_seed_ieta = tree->Branch("pho_seed_ieta", &pho_seed_ieta);
  TBranch *b_pho_seed_iphi = tree->Branch("pho_seed_iphi", &pho_seed_iphi);
  TBranch *b_pho_seed_eta = tree->Branch("pho_seed_eta", &pho_seed_eta);
  TBranch *b_pho_seed_phi = tree->Branch("pho_seed_phi", &pho_seed_phi);
  TBranch *b_pho_seed_raw_id = tree->Branch("pho_seed_raw_id", &pho_seed_raw_id);
  TBranch *b_pho_seed_hcal_ieta = tree->Branch("pho_seed_hcal_ieta", &pho_seed_hcal_ieta);
  TBranch *b_pho_seed_hcal_iphi = tree->Branch("pho_seed_hcal_iphi", &pho_seed_hcal_iphi);

  TBranch *b_hcalhit_ieta = tree->Branch("hcalhit_ieta", &hcalhit_ieta);
  TBranch *b_hcalhit_iphi = tree->Branch("hcalhit_iphi", &hcalhit_iphi);
  TBranch *b_hcalhit_energy = tree->Branch("hcalhit_energy", &hcalhit_energy);
  TBranch *b_hcalhit_seed_dieta = tree->Branch("hcalhit_seed_dieta", &hcalhit_seed_dieta);
  TBranch *b_hcalhit_seed_diphi = tree->Branch("hcalhit_seed_diphi", &hcalhit_seed_diphi);
  TBranch *b_hcalhit_raw_id = tree->Branch("hcalhit_raw_id", &hcalhit_raw_id);
  TBranch *b_hcalhit_depth = tree->Branch("hcalhit_depth", &hcalhit_depth);
  TBranch *b_hcalhit_pho_index = tree->Branch("hcalhit_pho_index", &hcalhit_pho_index);
  TBranch *b_hcalhit_eta = tree->Branch("hcalhit_eta", &hcalhit_eta);
  TBranch *b_hcalhit_phi = tree->Branch("hcalhit_phi", &hcalhit_phi);

  TBranch *b_hcalRechitIeta = tree->Branch("hcalRechitIeta", &hcalRechitIeta);
  TBranch *b_hcalRechitIphi = tree->Branch("hcalRechitIphi", &hcalRechitIphi);
  TBranch *b_hcalRechitEnergy = tree->Branch("hcalRechitEnergy", &hcalRechitEnergy);
  TBranch *b_hcalAllRechitEnergy = tree->Branch("hcalAllRechitEnergy", &hcalAllRechitEnergy);
  TBranch *b_hcalRechitAbsDIetaFromPhoSeed = tree->Branch("hcalRechitAbsDIetaFromPhoSeed", &hcalRechitAbsDIetaFromPhoSeed);
  TBranch *b_hcalRechitAbsDIphiFromPhoSeed = tree->Branch("hcalRechitAbsDIphiFromPhoSeed", &hcalRechitAbsDIphiFromPhoSeed);
  TBranch *b_hcalRechitRawID = tree->Branch("hcalRechitRawID", &hcalRechitRawID);
  TBranch *b_hcalRechitDepth = tree->Branch("hcalRechitDepth", &hcalRechitDepth);
  TBranch *b_hcalRechitDepth_noise = tree->Branch("hcalRechitDepth_noise", &hcalRechitDepth_noise);
  TBranch *b_hcalRechitNoise = tree->Branch("hcalRechitNoise", &hcalRechitNoise);
  TBranch *b_hcalRechitEta = tree->Branch("hcalRechitEta", &hcalRechitEta);
  TBranch *b_diffFromPhi  = tree->Branch("diffFromPhi", &diffFromPhi);
  TBranch *b_rechit_Noise  = tree->Branch("rechit_Noise", &rechit_Noise);
  TBranch *b_hcalRechitPhi = tree->Branch("hcalRechitPhi", &hcalRechitPhi);
  TBranch *b_hcalRechitEta_noise = tree->Branch("hcalRechitEta_noise", &hcalRechitEta_noise);
  TBranch *b_hcalRechitPhi_noise = tree->Branch("hcalRechitPhi_noise", &hcalRechitPhi_noise);
  TBranch *b_hcalAllRechitEta = tree->Branch("hcalAllRechitEta", &hcalAllRechitEta);
  TBranch *b_hcalAllRechitPhi = tree->Branch("hcalAllRechitPhi", &hcalAllRechitPhi);
  TBranch *b_hcal_diffFromPhi = tree->Branch("hcal_diffFromPhi",&hcal_diffFromPhi);
  if (n_pho_ != 0) {
    std::cout << "Electron block realloc to " << pho_golden.capacity() << "..." << std::endl;

    b_pho_eb->SetAddress(pho_eb.data());
    b_pho_ee->SetAddress(pho_ee.data());
    b_pho_gap_eb_ee->SetAddress(pho_gap_eb_ee.data());
    b_pho_gap_eb_eta->SetAddress(pho_gap_eb_eta.data());
    b_pho_gap_eb_phi->SetAddress(pho_gap_eb_phi.data());
    b_pho_gap_ee_dee->SetAddress(pho_gap_ee_dee.data());
    b_pho_gap_ee_ring->SetAddress(pho_gap_ee_ring.data());



    b_pho_golden->SetAddress(pho_golden.data());
    b_pho_unknown->SetAddress(pho_unknown.data());
    b_pho_bigbrem->SetAddress(pho_bigbrem.data());
    b_pho_gap->SetAddress(pho_gap.data());
    b_pho_badtrack->SetAddress(pho_badtrack.data());
    b_pho_showering->SetAddress(pho_showering.data());
    b_pho_gen_ecal->SetAddress(pho_gen_ecal.data());
    b_pho_track_fbrem->SetAddress(pho_track_fbrem.data());
    b_pho_sc_fbrem->SetAddress(pho_sc_fbrem.data());
    b_pho_nbrem->SetAddress(pho_nbrem.data());
    b_pho_genmatch->SetAddress(pho_genmatch.data());
    b_pho_dR_reco_gen->SetAddress(pho_dR_reco_gen.data());
    b_pho_pt_ratio_reco_gen->SetAddress(pho_pt_ratio_reco_gen.data());
    b_pho_sc_energy->SetAddress(pho_sc_energy.data());
    b_pho_sc_raw_energy->SetAddress(pho_sc_raw_energy.data());
    b_pho_ecal_energy->SetAddress(pho_ecal_energy.data());
    b_pho_seed_energy->SetAddress(pho_seed_energy.data());
    b_pho_seed_corr_energy->SetAddress(pho_seed_corr_energy.data());
    b_pho_cmssw_hoe->SetAddress(pho_cmssw_hoe.data());
    b_pho_cmssw_hoe_tower->SetAddress(pho_cmssw_hoe_tower.data());
    b_pho_cmssw_hoe_5x5->SetAddress(pho_cmssw_hoe_5x5.data());
    b_pho_sc_eta->SetAddress(pho_sc_eta.data());
    b_pho_sc_phi->SetAddress(pho_sc_phi.data());
    b_pho_pt->SetAddress(pho_pt.data());
    b_pho_eta->SetAddress(pho_eta.data());
    b_pho_phi->SetAddress(pho_phi.data());
    b_pho_sieie_5x5->SetAddress(pho_sieie_5x5.data());
    b_pho_r9_5x5->SetAddress(pho_r9_5x5.data());
    b_pho_ecal->SetAddress(pho_ecal.data());
   

    b_pho_pfiso_pho->SetAddress(pho_pfiso_pho.data());
    b_pho_pfiso_neu->SetAddress(pho_pfiso_neu.data());
    b_pho_pfiso_cha->SetAddress(pho_pfiso_cha.data());
    b_pho_pfiso_pu->SetAddress(pho_pfiso_pu.data());
    b_pho_pfiso_hcal->SetAddress(pho_pfiso_hcal.data());
    b_pho_pfiso_ecal->SetAddress(pho_pfiso_ecal.data());

    b_pho_detiso03_ecalhit->SetAddress(pho_detiso03_ecalhit.data());
    b_pho_detiso03_hcaltower1->SetAddress(pho_detiso03_hcaltower1.data());
    b_pho_detiso03_hcaltower2->SetAddress(pho_detiso03_hcaltower2.data());
    b_pho_detiso03_trk->SetAddress(pho_detiso03_trk.data());
    b_pho_detiso03_trk_heep->SetAddress(pho_detiso03_trk_heep.data());

    b_pho_seed_detid->SetAddress(pho_seed_detid.data());
    b_pho_seed_subdetid->SetAddress(pho_seed_subdetid.data());
    b_pho_seed_ieta->SetAddress(pho_seed_ieta.data());
    b_pho_seed_iphi->SetAddress(pho_seed_iphi.data());
    b_pho_seed_eta->SetAddress(pho_seed_eta.data());
    b_pho_seed_phi->SetAddress(pho_seed_phi.data());
    b_pho_seed_raw_id->SetAddress(pho_seed_raw_id.data());
    b_pho_seed_hcal_ieta->SetAddress(pho_seed_hcal_ieta.data());
    b_pho_seed_hcal_iphi->SetAddress(pho_seed_hcal_iphi.data());
  }

  if (n_hcalhit_ != 0) {
    std::cout << "Hcalhit block realloc to " << hcalhit_depth.capacity() << "..." << std::endl;

    b_hcalhit_ieta->SetAddress(hcalhit_ieta.data());
    b_hcalhit_iphi->SetAddress(hcalhit_iphi.data());
    b_hcalhit_energy->SetAddress(hcalhit_energy.data());
    b_hcalhit_seed_dieta->SetAddress(hcalhit_seed_dieta.data());
    b_hcalhit_seed_diphi->SetAddress(hcalhit_seed_diphi.data());
    b_hcalhit_raw_id->SetAddress(hcalhit_raw_id.data());
    b_hcalhit_depth->SetAddress(hcalhit_depth.data());
    b_hcalhit_pho_index->SetAddress(hcalhit_pho_index.data());
    b_hcalhit_eta->SetAddress(hcalhit_eta.data());
    b_hcalhit_phi->SetAddress(hcalhit_phi.data());

    b_hcalRechitIeta->SetAddress(hcalRechitIeta.data());
    b_hcalRechitIphi->SetAddress(hcalRechitIphi.data());
    b_hcalRechitEnergy->SetAddress(hcalRechitEnergy.data());
    b_hcalAllRechitEnergy->SetAddress(hcalAllRechitEnergy.data());
    b_hcalRechitAbsDIetaFromPhoSeed->SetAddress(hcalRechitAbsDIetaFromPhoSeed.data());
    b_hcalRechitAbsDIphiFromPhoSeed->SetAddress(hcalRechitAbsDIphiFromPhoSeed.data());
    b_hcalRechitEta->SetAddress(hcalRechitEta.data());
    b_diffFromPhi->SetAddress(diffFromPhi.data());
    b_rechit_Noise->SetAddress(rechit_Noise.data());
    b_hcalRechitPhi->SetAddress(hcalRechitPhi.data());
    b_hcalRechitEta_noise->SetAddress(hcalRechitEta_noise.data());
    b_hcalRechitPhi_noise->SetAddress(hcalRechitPhi_noise.data());
    b_hcalAllRechitEta->SetAddress(hcalAllRechitEta.data());
    b_hcalAllRechitPhi->SetAddress(hcalAllRechitPhi.data());
    b_hcal_diffFromPhi->SetAddress(hcal_diffFromPhi.data());
    b_hcalRechitRawID->SetAddress(hcalRechitRawID.data());
    b_hcalRechitDepth->SetAddress(hcalRechitDepth.data());
    b_hcalRechitDepth_noise->SetAddress(hcalRechitDepth_noise.data());
    b_hcalRechitNoise->SetAddress(hcalRechitNoise.data());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
HoEFlatPhoton::beginJob()
{
/*  file = new TFile(output.c_str(), "recreate");*/
  /*tree = new TTree("tree", "");
  tree->SetAutoSave(0);
  tree->SetImplicitMT(false);*/
  reallocate_setaddress();
}

// ------------ method called once each job just after ending the event loop  ------------
void
HoEFlatPhoton::endJob()
{
 // file->cd();
 // tree->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoEFlatPhoton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoEFlatPhoton);
