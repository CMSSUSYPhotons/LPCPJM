// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEvent.h
// 
/*

 Description: Objects definitions used for SusyNtuples

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEvent.h,v 1.1 2011/03/24 23:46:27 dwjang Exp $
//

#ifndef SusyEvent_h
#define SusyEvent_h

#include <vector>
#include <map>

#include <TLorentzVector.h>


namespace susy {

  const float etaGapBegin = 1.442;
  const float etaGapEnd = 1.556;
  const float etaGap = 1.499;
  const float etaMax = 2.5;


  class Particle {

  public:

    Particle()  { Init(); }
    ~Particle() { Init(); }
    void Init();

    int index;
    int motherIndex;
    int status;
    int pdgId;
    int charge;
    TVector3 vertex;
    TLorentzVector momentum;

  };


  class CorrMETData {

  public:
    CorrMETData() { Init(); }
    ~CorrMETData() { Init(); }

    void Init();

    Float_t  dmEx;             // for uncorrection, correctedEx - dmEx
    Float_t  dmEy;             // for uncorrection, correctedEy - dmEy
    Float_t  dsumEt;           // for uncorrection, correctedSumEt - dsumEt
    Float_t  dSignificance;    // for uncorrection, correctedSig - dSignificance

  };


  class MET {

  public:

    MET()  { Init(); }
    ~MET() { Init(); }
    void Init();

    Float_t met() const {  return mEt.Mod(); }
    Float_t metX() const { return mEt.X(); }
    Float_t metY() const { return mEt.Y(); }

    Float_t  sumEt;
    Float_t  significance;
    TVector2 mEt;
    TVector3 vertex;
    std::vector<susy::CorrMETData>  mEtCorr;

  };



  class Cluster {

  public:
    
    Cluster()  { Init(); }
    ~Cluster() { Init(); }
    void Init();

    Int_t    index;
    Int_t    nCrystals;
    Float_t  energy;
    TVector3 position;

  };


  class SuperCluster {
    
  public:
    
    SuperCluster()  { Init(); }
    ~SuperCluster() { Init(); }
    void Init();

    Int_t    index;
    Int_t    seedClusterIndex; // index in vector<Cluster> below
    Float_t  energy;
    Float_t  preshowerEnergy;
    Float_t  phiWidth;
    Float_t  etaWidth;
    TVector3 position;
    std::vector<Int_t> basicClusterIndices;
    
  };


  class Track {

  public:

    Track()  { Init(); }
    ~Track() { Init(); }
    void Init();

    bool isGsfTrack() { return (statusCode & (0x1 << 0) ); }
    bool InnerOk() {    return (statusCode & (0x1 << 1) ); }
    bool OuterOk() {    return (statusCode & (0x1 << 2) ); }

    // derived quantities
    float normChi2() const { return (ndof != 0) ? chi2/ndof : chi2*1e6; }
    float dxy() const { return (-vertex.X()*momentum.Py() + vertex.Y()*momentum.Px())/momentum.Pt(); }
    float d0() const { return -dxy(); }
    float phi() const { return momentum.Phi(); }
    bool loose() const {         return (quality & ( 0x1 << 0)); }
    bool tight() const {         return (quality & ( 0x1 << 1)); }
    bool highPurity() const {    return (quality & ( 0x1 << 2)); }
    bool confirmed() const {     return (quality & ( 0x1 << 3)); }
    bool goodIterative() const { return (confirmed() || highPurity()); }

    Int_t          index;
    Int_t          algorithm;
    Int_t          quality;
    Int_t          statusCode;
    Int_t          nHits;
    Float_t        chi2;
    Float_t        ndof;
    Float_t        charge;
    Float_t        error[5]; // qoverp, lambda, phi, dxy, dsz
    TVector3       vertex;
    TLorentzVector momentum;
    std::map<TString,TVector3> hitPositions;
    std::map<TString,TVector3>  extrapolatedPositions;
    
  };


  class Photon {

  public:

    Photon()  { Init(); }
    ~Photon() { Init(); }
    void Init();

    // fiducial bits
    bool isEB()           { return (fidBit & (0x1 << 0)); }
    bool isEE()           { return (fidBit & (0x1 << 1)); }
    bool isEBEtaGap()     { return (fidBit & (0x1 << 2)); }
    bool isEBPhiGap()     { return (fidBit & (0x1 << 3)); }
    bool isEERingGap()    { return (fidBit & (0x1 << 4)); }
    bool isEEDeeGap()     { return (fidBit & (0x1 << 5)); }
    bool isEBEEGap()      { return (fidBit & (0x1 << 6)); }

    Int_t          fidBit;
    Int_t          nPixelSeeds;
    Float_t        hadronicOverEm;
    Float_t        hadronicDepth1OverEm;
    Float_t        hadronicDepth2OverEm;
    Float_t        e1x2;
    Float_t        e1x5;
    Float_t        e2x5;
    Float_t        e3x3;
    Float_t        e5x5;
    Float_t        maxEnergyXtal;
    Float_t        sigmaEtaEta;
    Float_t        sigmaIetaIeta;
    Float_t        r1x5;
    Float_t        r2x5;
    Float_t        r9;

    Float_t        ecalRecHitSumEtConeDR04;
    Float_t        hcalTowerSumEtConeDR04;
    Float_t        hcalDepth1TowerSumEtConeDR04;
    Float_t        hcalDepth2TowerSumEtConeDR04;
    Float_t        trkSumPtSolidConeDR04;
    Float_t        trkSumPtHollowConeDR04;
    Int_t          nTrkSolidConeDR04;
    Int_t          nTrkHollowConeDR04;

    Float_t        ecalRecHitSumEtConeDR03;
    Float_t        hcalTowerSumEtConeDR03;
    Float_t        hcalDepth1TowerSumEtConeDR03;
    Float_t        hcalDepth2TowerSumEtConeDR03;
    Float_t        trkSumPtSolidConeDR03;
    Float_t        trkSumPtHollowConeDR03;
    Int_t          nTrkSolidConeDR03;
    Int_t          nTrkHollowConeDR03;

    Float_t        chargedHadronIso;
    Float_t        neutralHadronIso;
    Float_t        photonIso;

    Float_t        dist;
    Float_t        dcot;
    Float_t        radius;

    Int_t          superClusterIndex;
    Float_t        superClusterPreshowerEnergy;
    Float_t        superClusterPhiWidth;
    Float_t        superClusterEtaWidth;
    TVector3       caloPosition;
    TLorentzVector momentum;
    std::map<TString,UChar_t> idPairs;

  };


  class Electron {

  public:

    Electron()  { Init(); }
    ~Electron() { Init(); }
    void Init();

    // fiducial bits
    bool isEB() {        return (fidBit & (0x1 << 0)); }
    bool isEE() {        return (fidBit & (0x1 << 1)); }
    bool isEBEEGap() {   return (fidBit & (0x1 << 2)); }
    bool isEBEtaGap() {  return (fidBit & (0x1 << 3)); }
    bool isEBPhiGap() {  return (fidBit & (0x1 << 4)); }
    bool isEEDeeGap() {  return (fidBit & (0x1 << 5)); }
    bool isEERingGap() { return (fidBit & (0x1 << 6)); }
    bool isEBGap() {     return (isEBEtaGap() || isEBPhiGap()); }
    bool isEEGap() {     return (isEEDeeGap() || isEERingGap()); }
    bool isGap() {       return (isEBGap() || isEEGap() || isEBEEGap()); }

    // boolean variables packed in boolPack
    bool isGsfCtfScPixChargeConsistent() { return (boolPack & (0x1 << 0)); }
    bool isGsfScPixChargeConsistent() {    return (boolPack & (0x1 << 1)); }
    bool isGsfCtfChargeConsistent() {      return (boolPack & (0x1 << 2)); }
    bool ecalDrivenSeed() {                return (boolPack & (0x1 << 3)); }
    bool trackerDrivenSeed() {             return (boolPack & (0x1 << 4)); }
    bool ecalDriven() {                    return (boolPack & (0x1 << 5)); }
    bool passingCutBasedPreselection() {   return (boolPack & (0x1 << 6)); }
    bool passingMvaPreselection() {        return (boolPack & (0x1 << 7)); }
    bool ambiguous() {                     return (boolPack & (0x1 << 8)); }
    bool isEcalEnergyCorrected() {         return (boolPack & (0x1 << 9)); }
    bool isMomentumCorrected() {           return (boolPack & (0x1 << 10)); }
    bool convFlags() {                     return (boolPack & (0x1 << 11)); }

    Int_t          fidBit;
    Int_t          boolPack;
    Int_t          scPixCharge;

    Float_t        eSuperClusterOverP;
    Float_t        eSeedClusterOverP;
    Float_t        eSeedClusterOverPout;
    Float_t        eEleClusterOverPout;
    Float_t        deltaEtaSuperClusterTrackAtVtx;
    Float_t        deltaEtaSeedClusterTrackAtCalo;
    Float_t        deltaEtaEleClusterTrackAtCalo;
    Float_t        deltaPhiSuperClusterTrackAtVtx;
    Float_t        deltaPhiSeedClusterTrackAtCalo;
    Float_t        deltaPhiEleClusterTrackAtCalo;

    Float_t        shFracInnerHits;

    Float_t        sigmaEtaEta;
    Float_t        sigmaIetaIeta;
    Float_t        e1x5;
    Float_t        e2x5Max;
    Float_t        e5x5;
    Float_t        hcalDepth1OverEcal;          // hadronic energy on depth1 / em enrgy
    Float_t        hcalDepth2OverEcal;          // hadronic energy on depth2 / em enrgy
    Float_t        hcalOverEcal;                // hadronic energy / em energy

    Float_t        dr03TkSumPt;
    Float_t        dr03EcalRecHitSumEt;
    Float_t        dr03HcalDepth1TowerSumEt;
    Float_t        dr03HcalDepth2TowerSumEt;
    Float_t        dr03HcalTowerSumEt;

    Float_t        dr04TkSumPt;
    Float_t        dr04EcalRecHitSumEt;
    Float_t        dr04HcalDepth1TowerSumEt;
    Float_t        dr04HcalDepth2TowerSumEt;
    Float_t        dr04HcalTowerSumEt;

    Float_t        convDist;
    Float_t        convDcot;
    Float_t        convRadius;

    Float_t        mva;

    Int_t          bremClass;
    Float_t        fbrem;

    Float_t        ecalEnergy;                  // corrected
    Float_t        ecalEnergyError;             // correction error
    Float_t        trackMomentumError;
    Float_t        electronMomentumError;

    Int_t          gsfTrackIndex;
    Int_t          closestCtfTrackIndex;
    Int_t          electronClusterIndex;
    Int_t          superClusterIndex;

    // AtVtx, AtCalo
    std::map<TString,TVector3> trackPositions;
    // AtVtx, AtCalo, Out, AtEleClus, AtVtxWithConstraint
    std::map<TString,TLorentzVector> trackMomentums;

    TVector3       vertex;
    TLorentzVector momentum;
    std::map<TString,Float_t> idPairs;

  };



  class Muon {

  public:

    Muon()  { Init(); }
    ~Muon() { Init(); }
    void Init();

    // muon type
    bool isGlobalMuon() {     return (type & (0x1 << 1)); }
    bool isTrackerMuon() {    return (type & (0x1 << 2)); }
    bool isStandAloneMuon() { return (type & (0x1 << 3)); }
    bool isCaloMuon() {       return (type & (0x1 << 4)); }

    Int_t          type;
    Int_t          nMatches;
    Int_t          nValidHits;
    Int_t          nValidTrackerHits;
    Int_t          nValidMuonHits;
    Int_t          nChambers;
    Int_t          timeNDof;
    Int_t          timeDirection;
    Float_t        timeAtIp;
    Float_t        timeAtIpError;
    Float_t        caloCompatibility;
    Float_t        emEnergy;
    Float_t        hadEnergy;
    Float_t        trackIsoR03;
    Float_t        ecalIsoR03;
    Float_t        hcalIsoR03;
    Float_t        trackIsoR05;
    Float_t        ecalIsoR05;
    Float_t        hcalIsoR05;

    Int_t          trackIndex;             // tracker only
    Int_t          standAloneTrackIndex;   // muon detector only
    Int_t          combinedTrackIndex;     // combined
    TLorentzVector momentum;

    std::map<TString, UChar_t> idPairs;

  };



  class Tau {

  public:

    Tau()  { Init(); }
    ~Tau() { Init(); }
    void Init();

    bool IsCaloTau() {             return (status & (0x1 << 0)); }
    bool IsPFTau() {               return (status & (0x1 << 1)); }
    bool electronPreIDDecision() { return (status & (0x1 << 2)); }
    bool muonDecision() {          return (status & (0x1 << 3)); }

    Int_t status;
    Int_t decayMode;
    Int_t leadTrackIndex;
    Int_t leadParticleIndex;

    // from CaloTau
    Float_t leadTracksignedSipt;
    Float_t leadTrackHCAL3x3hitsEtSum;
    Float_t leadTrackHCAL3x3hottesthitDEta;
    Float_t signalTracksInvariantMass;
    Float_t TracksInvariantMass;
    Float_t isolationTracksPtSum;
    Float_t isolationECALhitsEtSum;
    Float_t maximumHCALhitEt;

    // from PFTau
    Float_t leadPFChargedHadrCandsignedSipt;
    Float_t isolationPFChargedHadrCandsPtSum;
    Float_t isolationPFGammaCandsEtSum;
    Float_t maximumHCALPFClusterEt;
    Float_t emFraction;
    Float_t hcalTotOverPLead;
    Float_t hcalMaxOverPLead;
    Float_t hcal3x3OverPLead;
    Float_t ecalStripSumEOverPLead;
    Float_t bremsRecoveryEOverPLead;
    Float_t electronPreIDOutput;
    Float_t caloComp;
    Float_t segComp;

    Float_t trackIso;
    Float_t ecalIso;
    Float_t hcalIso;

    TLorentzVector momentum;

    std::map<TString,Float_t> idPairs;

    std::vector<Int_t> sigTracks;
    std::vector<Int_t> isoTracks;
    std::vector<susy::Particle> sigParticles;
    std::vector<susy::Particle> isoParticles;

  };


  class CaloJet {

  public:
    
    CaloJet()  { Init(); }
    ~CaloJet() { Init(); }
    void Init();

    // Basic Jet Info
    Float_t        partonFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    Int_t          nPasses;
    Int_t          nConstituents;

    // CaloJet info
    Float_t        maxEInEmTowers;
    Float_t        maxEInHadTowers;
    Float_t        energyFractionHadronic;
    Float_t        emEnergyFraction;
    Float_t        hadEnergyInHB;
    Float_t        hadEnergyInHO;
    Float_t        hadEnergyInHE;
    Float_t        hadEnergyInHF;
    Float_t        emEnergyInEB;
    Float_t        emEnergyInEE;
    Float_t        emEnergyInHF;
    Float_t        towersArea;
    Int_t          n90;
    Int_t          n60;

    // Jet ID info
    Float_t        fHPD;
    Float_t        fRBX;
    Float_t        n90Hits;
    Float_t        fSubDetector1;
    Float_t        fSubDetector2;
    Float_t        fSubDetector3;
    Float_t        fSubDetector4;
    Float_t        restrictedEMF;
    Int_t          nHCALTowers;
    Int_t          nECALTowers;
    Float_t        approximatefHPD;
    Float_t        approximatefRBX;
    Int_t          hitsInN90;
    Int_t          numberOfHits2RPC;
    Int_t          numberOfHits3RPC;
    Int_t          numberOfHitsRPC;

    TVector3       vertex;
    TLorentzVector momentum;
    TLorentzVector detectorP4;

    // JES correction factor is stored the following order
    // For example, L4 stands for correction applied up to L4 from Raw
    // Raw, L1, L2, L3, L4, L5g, L5uds, L5c, L5b, L6g, L6uds, L6c, L6b, L7g, L7uds, L7c, L7b
    std::map<TString,Float_t> jesMap;

    // Btag discriminator name and value
    std::map<TString,Float_t> bTagMap;

  };


  class PFJet {

  public:
    
    PFJet()  { Init(); }
    ~PFJet() { Init(); }
    void Init();

    // Basic Jet Info
    Float_t        partonFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    Int_t          nPasses;
    Int_t          nConstituents;

    Float_t        chargedHadronEnergy;
    Float_t        neutralHadronEnergy;
    Float_t        photonEnergy;
    Float_t        electronEnergy;
    Float_t        muonEnergy;
    Float_t        HFHadronEnergy;
    Float_t        HFEMEnergy;
    Float_t        chargedEmEnergy;
    Float_t        chargedMuEnergy;
    Float_t        neutralEmEnergy;
    Int_t          chargedHadronMultiplicity;
    Int_t          neutralHadronMultiplicity;
    Int_t          photonMultiplicity;
    Int_t          electronMultiplicity;
    Int_t          muonMultiplicity;
    Int_t          HFHadronMultiplicity;
    Int_t          HFEMMultiplicity;
    Int_t          chargedMultiplicity;
    Int_t          neutralMultiplicity;

    // Jet ID info
    Float_t        fHPD;
    Float_t        fRBX;
    Float_t        n90Hits;
    Float_t        fSubDetector1;
    Float_t        fSubDetector2;
    Float_t        fSubDetector3;
    Float_t        fSubDetector4;
    Float_t        restrictedEMF;
    Int_t          nHCALTowers;
    Int_t          nECALTowers;
    Float_t        approximatefHPD;
    Float_t        approximatefRBX;
    Int_t          hitsInN90;
    Int_t          numberOfHits2RPC;
    Int_t          numberOfHits3RPC;
    Int_t          numberOfHitsRPC;

    TVector3       vertex;
    TLorentzVector momentum;

  };


  class JPTJet {

  public:
    
    JPTJet()  { Init(); }
    ~JPTJet() { Init(); }
    void Init();

    // Basic Jet Info
    Float_t        partonFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    Int_t          nPasses;
    Int_t          nConstituents;

    Float_t        chargedHadronEnergy;
    Float_t        neutralHadronEnergy;
    Float_t        chargedEmEnergy;
    Float_t        neutralEmEnergy;
    Int_t          chargedMultiplicity;
    Int_t          muonMultiplicity;
    Int_t          elecMultiplicity;
    Float_t        getZSPCor;

    TVector3       vertex;
    TLorentzVector momentum;

  };


  typedef std::vector<susy::CaloJet> CaloJetCollection;
  typedef std::vector<susy::PFJet> PFJetCollection;
  typedef std::vector<susy::JPTJet> JPTJetCollection;

  class Event {

  public:

    Event()  { Init(); }
    ~Event() { Init(); }

    // Initialize members
    void Init();

    // Members are made as public intentionally for easy access

    UChar_t                                     isRealData;
    Int_t                                       runNumber;
    ULong_t                                     eventNumber;
    Int_t                                       luminosityBlockNumber;
    Float_t                                     avgInsRecLumi;
    Float_t                                     intgRecLumi;
    Int_t                                       cosmicFlag;

    TVector3                                    beamSpot;

    std::map<TString,UChar_t>                   l1Map;
    std::map<TString,UChar_t>                   hltMap;
    std::map<TString,susy::MET>                 metMap;

    std::vector<TVector3>                       vertices;
    std::vector<susy::Track>                    tracks;          // only selected tracks associated with objects
    std::vector<susy::SuperCluster>             superClusters;   // only selected super clusters associated with objects
    std::vector<susy::Cluster>                  clusters;        // only selected basic clusters associated with super clusters
    std::vector<susy::Muon>                     muons;
    std::vector<susy::Photon>                   photons;
    std::vector<susy::Electron>                 electrons;
    std::map<TString,susy::CaloJetCollection>   caloJets;
    std::map<TString,susy::PFJetCollection>     pfJets;
    std::map<TString,susy::JPTJetCollection>    jptJets;

    // optional collections
    std::vector<susy::Tau>                      taus;            // not stored by default
    std::vector<susy::Track>                    generalTracks;   // not stored by default

    // generated information. Valid only for isRealData == 0, i.e. MC
    std::vector<TVector3>                       genVertices;
    std::vector<susy::Particle>                 genParticles;
    
  };




} // namespace susy

#endif