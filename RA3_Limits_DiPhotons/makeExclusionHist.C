#include <set>
#include <cmath>
#include <vector>
#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>

#include "GridPoint.C"


void InterpolatePotHoles(TH2D *h) {
  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 1e-10;

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      if(val > epsilon) continue;
      int ncnt = 0;
      double sum = 0;
      double sumErr = 0;
      double up    = h->GetBinContent(ix,iy+1);
      if(up > epsilon){
	sum += up;
	sumErr += h->GetBinError(ix,iy+1);
	ncnt++;
      }
      double down  = h->GetBinContent(ix,iy-1);
      if(down > epsilon){
	sum += down;
	sumErr += h->GetBinError(ix,iy-1);
	ncnt++;
      }
      double left  = h->GetBinContent(ix-1,iy);
      if(left > epsilon){
	sum += left;
	sumErr += h->GetBinError(ix-1,iy);
	ncnt++;
      }
      double right = h->GetBinContent(ix+1,iy);
      if(right > epsilon){
	sum += right;
	sumErr += h->GetBinError(ix+1,iy);
	ncnt++;
      }
      h->SetBinContent(ix,iy,sum/ncnt);
      h->SetBinError(ix,iy,std::sqrt(sumErr)/ncnt);
    } // for iy
  } // for ix

}



void makeExclusionHist(TString bino="bino", TString mNu="375", TString metCut="met100", TString jetCut="1jet", TString model="1") {

  TString label = bino + "_mN" + mNu + "_" + metCut + "_" + jetCut + "_model" + model;

  int debugLevel = 1;

  GridPoint* gp = new GridPoint;

  int mN = atoi(mNu.Data());

  // for full grids
//   const int nG = 21;
//   const int nS = 21;
//   float mGBins[nG+1] = {400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};
//   float mSBins[nS+1] = {400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};

  const int nG = 24;
  const int nS = 24;
  float mGBins[nG+1] = {160, 240, 320, 400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};
  float mSBins[nS+1] = {160, 240, 320, 400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};


  TFile* fin = new TFile("tree_"+label+".root","READ");
  TTree* tree = (TTree*) fin->Get("gridPoints");
  tree->SetBranchAddress("GridPoint",&gp);

  std::vector<GridPoint> grid;
  int nevent = tree->GetEntries();
  for(int ievt = 0; ievt < nevent; ievt++) {
    tree->GetEntry(ievt);
    grid.push_back(*gp);
  }

  char datafile[200];
  sprintf(datafile,"%s%sNLOxsec.dat",bino.Data(),mNu.Data());
  GetXSection(grid,TString(datafile));             // XS and RS errors

  TFile* fout = new TFile("exclusion_results_"+label+".root","RECREATE");
  TH2D* h_acc = new TH2D("acc",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_xsec = new TH2D("xsec",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_xsec_1L = new TH2D("xsec_1L",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_xsec_1H = new TH2D("xsec_1H",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_xsec_2L = new TH2D("xsec_2L",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_xsec_2H = new TH2D("xsec_2H",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);

  TH2D* h_limit = new TH2D("limit",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_exp = new TH2D("exp",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_exp_1L = new TH2D("exp_1L",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_exp_1H = new TH2D("exp_1H",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_exp_2L = new TH2D("exp_2L",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);
  TH2D* h_exp_2H = new TH2D("exp_2H",";M_{#tilde{q}} (GeV);M_{#tilde{g}} (GeV)",nS,mSBins,nG,mGBins);

  std::set<int> mGs;
  std::set<int> mSs;
  std::set<int> mNs;

  std::set<int> mGs_k0;
  std::set<int> mSs_k0;
  std::set<int> mNs_k0;

  std::set<int> mGs_x0;
  std::set<int> mSs_x0;
  std::set<int> mNs_x0;

  double minPDFerror = 999;
  double maxPDFerror = 0;
  double minRSerror = 999;
  double maxRSerror = 0;
  double minPDFAccError = 999;
  double maxPDFAccError = 0;

  nevent = grid.size();
  for(int ievt = 0; ievt < nevent; ievt++) {
    GridPoint* p = &(grid[ievt]);

    //    p->Print();

    mGs.insert(p->mGluino);
    mSs.insert(p->mSquark);
    mNs.insert(p->mNeutralino);

    // look only one Neutralino point
    if(p->mNeutralino != mN) continue;

    h_acc->Fill(p->mSquark,p->mGluino,p->accValue);
    h_xsec->Fill(p->mSquark,p->mGluino,p->xsecValue);

    double oneSigma_L = std::sqrt(p->xsecRSErrorNeg * p->xsecRSErrorNeg + p->xsecPDFError * p->xsecPDFError);
    double oneSigma_H = std::sqrt(p->xsecRSErrorPos * p->xsecRSErrorPos + p->xsecPDFError * p->xsecPDFError);
    h_xsec_1L->Fill(p->mSquark,p->mGluino,p->xsecValue - p->xsecValue*oneSigma_L);
    h_xsec_1H->Fill(p->mSquark,p->mGluino,p->xsecValue + p->xsecValue*oneSigma_H);
    h_xsec_2L->Fill(p->mSquark,p->mGluino,p->xsecValue - p->xsecValue*2*oneSigma_L);
    h_xsec_2H->Fill(p->mSquark,p->mGluino,p->xsecValue + p->xsecValue*2*oneSigma_H);


    h_limit->Fill(p->mSquark,p->mGluino,p->limit);
    h_exp->Fill(p->mSquark,p->mGluino,p->explimit);
    h_exp_1L->Fill(p->mSquark,p->mGluino,p->explimit_1L);
    h_exp_1H->Fill(p->mSquark,p->mGluino,p->explimit_1H);
    h_exp_2L->Fill(p->mSquark,p->mGluino,p->explimit_2L);
    h_exp_2H->Fill(p->mSquark,p->mGluino,p->explimit_2H);

    if(p->mGluino == 1920 && p->mSquark == 800 && p->mNeutralino == 375){
      std::cout << "last year's included point, ";
      std::cout << "xsec : " << p->xsecValue << ", obs : " << p->limit << std::endl;
      p->Print();
    }
    if(p->mGluino == 720 && p->mSquark == 720 && p->mNeutralino == 375){
      std::cout << "last year's excluded point, ";
      std::cout << "xsec : " << p->xsecValue << ", obs : " << p->limit << std::endl;
      p->Print();
    }

    if(p->accErrorPDF > maxPDFAccError) maxPDFAccError = p->accErrorPDF;
    if(p->accErrorPDF < minPDFAccError) minPDFAccError = p->accErrorPDF;

    if(p->xsecPDFError > maxPDFerror) maxPDFerror = p->xsecPDFError;
    if(p->xsecPDFError < minPDFerror) minPDFerror = p->xsecPDFError;

    double smaller = TMath::Min(p->xsecRSErrorNeg,p->xsecRSErrorPos);
    double bigger = TMath::Max(p->xsecRSErrorNeg,p->xsecRSErrorPos);
    if(bigger > maxRSerror) maxRSerror = bigger;
    if(smaller < minRSerror) minRSerror = smaller;

  }// for

  int nbinsX = h_xsec->GetXaxis()->GetNbins();
  int nbinsY = h_xsec->GetYaxis()->GetNbins();
  for(int i=1; i<=nbinsX; i++) {
    for(int j=1; j<=nbinsY; j++) {
      if(h_xsec->GetBinContent(i,j) < 1e-9) {
	std::cout << "mS, mG : " << h_xsec->GetXaxis()->GetBinLowEdge(i) << ", " << h_xsec->GetYaxis()->GetBinLowEdge(j) << std::endl;
      }
    }
  }

  std::cout << "mGs : " << mGs.size() << std::endl;
  std::cout << "{";
  for(std::set<int>::iterator it = mGs.begin(); it != mGs.end(); it++) {
    std::cout << *it << ", ";
  }
  std::cout << "};" << std::endl;

  std::cout << "mSs : " << mSs.size() << std::endl;
  std::cout << "{";
  for(std::set<int>::iterator it = mSs.begin(); it != mSs.end(); it++) {
    std::cout << *it << ", ";
  }
  std::cout << "};" << std::endl;

  std::cout << "mNs : " << mNs.size() << std::endl;
  std::cout << "{";
  for(std::set<int>::iterator it = mNs.begin(); it != mNs.end(); it++) {
    std::cout << *it << ", ";
  }
  std::cout << "};" << std::endl << std::endl;


  std::cout << "PDF errors on XS  : " << minPDFerror << " - " << maxPDFerror << std::endl;
  std::cout << "RS errors on XS   : " << minRSerror << " - " << maxRSerror << std::endl;
  std::cout << "PDF errors on Acc : " << minPDFAccError << " - " << maxPDFAccError << std::endl;

  InterpolatePotHoles(h_acc);
  InterpolatePotHoles(h_xsec);
  InterpolatePotHoles(h_xsec_1L);
  InterpolatePotHoles(h_xsec_1H);
  InterpolatePotHoles(h_xsec_2L);
  InterpolatePotHoles(h_xsec_2H);

  InterpolatePotHoles(h_limit);
  InterpolatePotHoles(h_exp);
  InterpolatePotHoles(h_exp_1L);
  InterpolatePotHoles(h_exp_1H);
  InterpolatePotHoles(h_exp_2L);
  InterpolatePotHoles(h_exp_2H);

  std::cout << "getentries : " << h_limit->GetEntries() << std::endl;

  fout->cd();
  fout->Write();
  fout->Close();

  fin->Close();

}

