import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
process = cms.Process("HoE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load("Geometry.HcalEventSetup.CaloTowerTopology_cfi")
process.load("Configuration.Geometry.GeometryExtended2021Reco_cff")
#process.load("Configuration.Geometry.GeometryExtended2017_cff")
#process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
options = VarParsing.VarParsing ('analysis')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '') # 2018 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2023_realistic_v3', '') # 2023 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '') # 2021 MC
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#options.inputFiles_load = runover.list
options.inputFiles= 'root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_28.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_26.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_21.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_14.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_18.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-2.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-3.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-5.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-6.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-7.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-8.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-9.root','root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_1-10.root'
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring( options.inputFiles
                                # example 2021 file
                                #'/store/mc/Run3Summer19MiniAOD/QCD_Pt_1800to2400_TuneCP5_14TeV_pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/2C20FC5F-B670-4C46-BCC9-02EDF6EFA5F6.root'
#                                'root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_28.root'
                              #  'root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_26.root',
                               # 'root://cmsxrootd.fnal.gov//store/user/lpcqstar/Jyoti/EGamma_Studies/Photon/PhotonFlatPt0To200/crab_EGamma_Photon/200509_122216/0000/TSG-Run3Winter20DRMiniAOD-00007_21.root'
                                # 'file:TSG-Run3Winter20DRMiniAOD-00007_10000events.root' 
                            )
)
#root://eoscms.cern.ch//eos/cms
process.load("RecoEgamma.PhotonIdentification.photonIDValueMapProducer_cfi")
#/store/group/lpcqstar/Jyoti/EGamma_Studies/Photon/Ntuples/
process.demo = cms.EDAnalyzer('Rechits_check',
                              photons = cms.InputTag('slimmedPhotons'),
                              pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                              #hbheInput = cms.InputTag("reducedEgamma" ,  "reducedHBHEHits" ,  "PAT"),
                              hbheInput = cms.InputTag("hbhereco"),
                              ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                              eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                              esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                              genParticleSrc       = cms.InputTag("prunedGenParticles"),
                              #genParticleSrc       = cms.InputTag("genParticles"),
                              Run2_2018_ = cms.bool(False),
                              rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
                              photonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                              neutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                              chargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
)
process.TFileService = cms.Service("TFileService", fileName = cms.string('output_file_AllRechit_remReduced.root'))

process.p = cms.Path(process.photonIDValueMapProducer * process.demo)
