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



void makeExclusionHistScan(TString label="bino_mNScan_met100_1jet_model1", TString fixed="q") {

  int debugLevel = 1;

  GridPoint* p = new GridPoint;

  const int nX = 10;
  const int nY = 24;
  float xBins[nY+1] = {150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150};
  float yBins[nY+1] = {160, 240, 320, 400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};

  TFile* fin = new TFile("tree_"+label+".root","READ");
  TTree* tree = (TTree*) fin->Get("gridPoints");
  tree->SetBranchAddress("GridPoint",&p);

  TString ylabel = "q";
  if(fixed.Contains("q")) ylabel = "g";

  TFile* fout = new TFile("exclusion_results_"+label+"_"+fixed+".root","RECREATE");
  TH2D* h_acc = new TH2D("acc",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_xsec = new TH2D("xsec",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_xsec_1L = new TH2D("xsec_1L",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_xsec_1H = new TH2D("xsec_1H",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_xsec_2L = new TH2D("xsec_2L",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_xsec_2H = new TH2D("xsec_2H",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);

  TH2D* h_limit = new TH2D("limit",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_exp = new TH2D("exp",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_exp_1L = new TH2D("exp_1L",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_exp_1H = new TH2D("exp_1H",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_exp_2L = new TH2D("exp_2L",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);
  TH2D* h_exp_2H = new TH2D("exp_2H",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{"+ylabel+"}} (GeV/c^{2})",nX,xBins,nY,yBins);

  std::set<int> mGs;
  std::set<int> mSs;
  std::set<int> mNs;

  double minPDFerror = 999;
  double maxPDFerror = 0;
  double minRSerror = 999;
  double maxRSerror = 0;
  double minPDFAccError = 999;
  double maxPDFAccError = 0;

  int nevent = tree->GetEntries();
  for(int ievt = 0; ievt < nevent; ievt++) {
    tree->GetEntry(ievt);
    //    p->Print();

    if(fixed.Contains("q") && p->mSquark != 2500) continue;
    if(fixed.Contains("g") && p->mGluino != 2500) continue;

    mGs.insert(p->mGluino);
    mSs.insert(p->mSquark);
    mNs.insert(p->mNeutralino);

    double x = p->mNeutralino;
    double y = p->mSquark;
    if(fixed.Contains("q")) y = p->mGluino;

    h_acc->Fill(x,y,p->accValue);
    h_xsec->Fill(x,y,p->xsecValue);

    double oneSigma_L = std::sqrt(p->xsecRSErrorNeg * p->xsecRSErrorNeg + p->xsecPDFError * p->xsecPDFError);
    double oneSigma_H = std::sqrt(p->xsecRSErrorPos * p->xsecRSErrorPos + p->xsecPDFError * p->xsecPDFError);
    h_xsec_1L->Fill(x,y,p->xsecValue - p->xsecValue*oneSigma_L);
    h_xsec_1H->Fill(x,y,p->xsecValue + p->xsecValue*oneSigma_H);
    h_xsec_2L->Fill(x,y,p->xsecValue - p->xsecValue*2*oneSigma_L);
    h_xsec_2H->Fill(x,y,p->xsecValue + p->xsecValue*2*oneSigma_H);


    h_limit->Fill(x,y,p->limit);
    h_exp->Fill(x,y,p->explimit);
    h_exp_1L->Fill(x,y,p->explimit_1L);
    h_exp_1H->Fill(x,y,p->explimit_1H);
    h_exp_2L->Fill(x,y,p->explimit_2L);
    h_exp_2H->Fill(x,y,p->explimit_2H);

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

