import FWCore.ParameterSet.Config as cms

susyNtuplizer = cms.EDAnalyzer('SusyNtuplizer',
                               lumiSummaryTag = cms.InputTag("lumiProducer"),
                               l1GTReadoutTag = cms.InputTag("gtDigis"),
                               l1GTObjectMapTag = cms.InputTag("hltL1GtObjectMap"),
                               hltCollectionTag = cms.InputTag("TriggerResults::HLT"),
                               vtxCollectionTag = cms.InputTag("offlinePrimaryVertices"),
                               trackCollectionTag = cms.InputTag("generalTracks"),
                               muonCollectionTag = cms.InputTag("muons"),
                               muonIDCollectionTags = cms.vstring("muidAllArbitrated",
                                                                  "muidGMStaChiCompatibility",
                                                                  "muidGMTkChiCompatibility",
                                                                  "muidGMTkKinkTight",
                                                                  "muidGlobalMuonPromptTight",
                                                                  "muidTM2DCompatibilityLoose",
                                                                  "muidTM2DCompatibilityTight",
                                                                  "muidTMLastStationAngLoose",
                                                                  "muidTMLastStationAngTight",
                                                                  "muidTMLastStationLoose",
                                                                  "muidTMLastStationOptimizedLowPtLoose",
                                                                  "muidTMLastStationOptimizedLowPtTight",
                                                                  "muidTMLastStationTight",
                                                                  "muidTMOneStationAngLoose",
                                                                  "muidTMOneStationAngTight",
                                                                  "muidTMOneStationLoose",
                                                                  "muidTMOneStationTight",
                                                                  "muidTrackerMuonArbitrated"),
                               electronCollectionTags = cms.vstring("gsfElectrons"),
                               electronIDCollectionTags = cms.vstring("eidLoose",
                                                                      "eidRobustHighEnergy",
                                                                      "eidRobustLoose",
                                                                      "eidRobustTight",
                                                                      "eidTight"),
                               photonCollectionTags = cms.vstring("photons","pfPhotonTranslator:pfphot"),
                               photonIDCollectionTags = cms.vstring("PhotonCutBasedIDLoose",
                                                                    "PhotonCutBasedIDLooseEM",
                                                                    "PhotonCutBasedIDTight"),
                               genCollectionTag = cms.InputTag("genParticles"),
                               simVertexCollectionTag = cms.InputTag("g4SimHits"),
                               caloJetCollectionTags = cms.vstring("ak5CaloJets"),
                               pfJetCollectionTags = cms.vstring("ak5PFJets"),
                               jptJetCollectionTags = cms.vstring("JetPlusTrackZSPCorJetAntiKt5"),
                               metCollectionTags = cms.vstring("pfMet","genMetTrue"),
                               pfCandidateCollectionTags = cms.vstring("particleFlow","particleFlow:electrons"),
                               puSummaryInfoTag = cms.InputTag("addPileupInfo"),
                               muonThreshold = cms.double(20.0),
                               electronThreshold = cms.double(20.0),
                               photonThreshold = cms.double(20.0),
                               jetThreshold = cms.double(20.0),
                               pfParticleThreshold = cms.double(1.0),
                               debugLevel = cms.int32(0),
                               storeGenInfos = cms.bool(True),
                               storeGeneralTracks = cms.bool(False),
                               recoMode = cms.bool(False),
                               outputFileName = cms.string("susyEvents.root")                               
)
