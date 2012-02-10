#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>

const int NCH = 5;
const double bins[NCH] = {50,60,70,80,100};
const double epsilon = 1e-10;
const double luminosity = 4680;
const double DataMCScale = 0.99; // 0.99 +/- 0.04

//const TString datacard_dir = "/uscms/home/dwjang/work/jobs/datacards/20120131";
const TString datacard_dir = "/uscms/home/dwjang/work/jobs/datacards/20120203";

class BinInfo {

public:

  BinInfo() { x = y = error = 0; }
  ~BinInfo() {}

  int x;
  double y;
  double error;

};

class GridPoint {

public:

  GridPoint() {Init();}
  ~GridPoint() {Init();}

  void Init();

  int mG;
  int mS;
  int mN;

  int ngen;
  double acc;

  double lumi;
  double xsecValue;       // NLO XS from Prospino in pb
  double xsecPDFError;    // PDF error
  double xsecRSErrorNeg;  // renormalization scale error
  double xsecRSErrorPos;  // renormalization scale error
  double accErrorPDF;     // acceptance errors on XS PDF
  double lumi_sysError;
  double qcd_sysError;
  double ew_sysError;
  double sig_sysError;

  std::vector<BinInfo> ggBins;
  std::vector<BinInfo> ewBins;
  std::vector<BinInfo> qcdBins;
  std::vector<BinInfo> qcdSysErrors;
  std::vector<BinInfo> sig_ggBins;
  std::vector<BinInfo> sig_ffBins;
  std::vector<BinInfo> sigBins;

  double limit;           // observed limit
  double explimit;        // expected limit
  double explimit_1L;     // expected limit -1 sigma
  double explimit_1H;     // expected limit +1 sigma
  double explimit_2L;     // expected limit -2 sigma
  double explimit_2H;     // expected limit +2 sigma

};

void GridPoint::Init() {
  mG = mS = mN = 0;
  ngen = 10000;
  acc = 0;
  lumi = luminosity; // int. luminosity
  xsecValue = 0;
  xsecPDFError = 0;
  xsecRSErrorNeg = 0;
  xsecRSErrorPos = 0;
  accErrorPDF = 0;

  lumi_sysError = 1.045;
  ew_sysError = 1.34; // fakerate error
  qcd_sysError = 1.25; // shape difference between ee and ff, will be calculated bin by bin later...

  // 0.5% (0.005) for electron/photon difference
  // 2.4% (0.02) for pileup
  // 1.9% (0.02) for signal fit over/underestimation
  // 1.8% (0.02) for signal/background shape assumption
  // total error of 0.04 (4%), and the scale factor is 0.99 +/- 0.04
  sig_sysError = 1.04;

  ggBins.clear();
  ewBins.clear();
  qcdBins.clear();
  qcdSysErrors.clear();
  sig_ggBins.clear();
  sig_ffBins.clear();
  sigBins.clear();

  limit = explimit = explimit_1L = explimit_1H = explimit_2L = explimit_2H = 0;
}


void readData(TFile* f, TString jet,
	      std::vector<BinInfo>& ggBins,
	      std::vector<BinInfo>& qcdBins,
	      std::vector<BinInfo>& ewBins,
	      std::vector<BinInfo>& qcdSysErrors);
void readSig(TFile* f, TString jet, int mS, int mG, int mN,
	     std::vector<BinInfo>& sig_ggBins,
	     std::vector<BinInfo>& sig_ffBins,
	     std::vector<BinInfo>& sigBins);
void getBins(TH1F* h, std::vector<BinInfo>& binInfos);
void printBins(std::vector<BinInfo>& binInfos);
void printBins(std::fstream& of, std::vector<BinInfo>& binInfos);
void GetNGen(std::vector<GridPoint>& grid, TString datafile);
void GetXSection(std::vector<GridPoint>& grid, TString datafile);
void GetPDFErrorsXSection(std::vector<GridPoint>& grid, TString datafile);
void GetPDFErrorsAcceptance(std::vector<GridPoint>& grid, TString datafile);
void makeSignalGains(std::vector<GridPoint>& grid);
void makeDataCard(std::vector<GridPoint>& grid, TString bino, TString jet);
void makeDataCardSingleChannel(std::vector<GridPoint>& grid, TString bino, TString jet); // datacard for MET > 100 GeV
void makeDataCardEachChannel(std::vector<GridPoint>& grid, TString bino, TString jet); // datacard for individual bins

void makeTemplate(TString bino = "bino", TString jet="1jet") {

  TString hist_dir = "inputHists.20120130";
  TString dataHist = hist_dir + "/" + "limit_setting_0130.root";
  TString sigHist = hist_dir + "/" + "signal_contamination_" + bino + "_chi0375.root";
  if(bino.Contains("mNScan")) sigHist = hist_dir + "/" + "signal_contamination_bino_chi0.root";

  TFile* fData = new TFile(dataHist,"READ");
  std::vector<BinInfo> ggBins;
  std::vector<BinInfo> qcdBins;
  std::vector<BinInfo> ewBins;
  std::vector<BinInfo> qcdSysErrors;
  readData(fData,jet,ggBins,qcdBins,ewBins,qcdSysErrors);

  TFile* fSig = new TFile(sigHist,"READ");

  std::vector<GridPoint> grids;

  int npoint = 0;

  if(bino.Contains("mNScan")) {
    for(int iN=150; iN<=1050; iN+=100) {
      for(int iG=160; iG<=2000; iG+=80) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = 2500;
	grid.mN = iN;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	std::vector<BinInfo> sig_ggBins;
	std::vector<BinInfo> sig_ffBins;
	std::vector<BinInfo> sigBins;
	readSig(fSig,jet,2500,iG,iN,sig_ggBins,sig_ffBins,sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sigBins = sigBins;
	grids.push_back(grid);
      }// iG
    }// iN
    GetXSection(grids,"xsecdat/binoNLOxsec_mNScan.dat"); // get XS and RS error
    GetPDFErrorsXSection(grids,"xsecdat/xsectionPDFErrors_mNScan.dat");    // PDF errors on XS
    GetPDFErrorsAcceptance(grids,"xsecdat/acceptancePDFErrors_mNScan.dat");  // PDF errors on acceptance
  }
  else {
    for(int iS=400; iS<=2000; iS+=80) {
      for(int iG=400; iG<=2000; iG+=80) {
	npoint++;
	GridPoint grid;
	grid.mG = iG;
	grid.mS = iS;
	grid.mN = 375;
	grid.ngen = 10000;
	grid.lumi = luminosity;

	std::vector<BinInfo> sig_ggBins;
	std::vector<BinInfo> sig_ffBins;
	std::vector<BinInfo> sigBins;
	readSig(fSig,jet,iS,iG,375,sig_ggBins,sig_ffBins,sigBins);

	grid.ggBins = ggBins;
	grid.qcdBins = qcdBins;
	grid.ewBins = ewBins;
	grid.qcdSysErrors = qcdSysErrors;
	grid.sig_ggBins = sig_ggBins;
	grid.sig_ffBins = sig_ffBins;
	grid.sigBins = sigBins;
	grids.push_back(grid);
      }// iG
    }// iS
    GetXSection(grids,"xsecdat/"+bino+"375NLOxsec.dat"); // get XS and RS error
    GetPDFErrorsXSection(grids,"xsecdat/xsectionPDFErrors.dat");    // PDF errors on XS
    GetPDFErrorsAcceptance(grids,"xsecdat/acceptancePDFErrors.dat");  // PDF errors on acceptance
  }

  if(bino.Contains("wino")) GetNGen(grids,"xsecdat/mcAccMap_wino_mN375.dat");

  makeSignalGains(grids);

  makeDataCardEachChannel(grids,bino,jet);
  makeDataCard(grids,bino,jet);

  std::cout << "npoint : " << npoint << std::endl;

}



void readData(TFile* f, TString jet,
	      std::vector<BinInfo>& ggBins,
	      std::vector<BinInfo>& qcdBins,
	      std::vector<BinInfo>& ewBins,
	      std::vector<BinInfo>& qcdSysErrors) {

  TH1F* gg = (TH1F*) f->Get("met_gg_"+jet);
  ggBins.clear();
  getBins(gg,ggBins);
  std::cout << "gg  events ----------" << std::endl;
  printBins(ggBins);

  //  TH1F* qcd = (TH1F*) f->Get("met_qcd_avg_"+jet);
  TH1F* qcd = (TH1F*) f->Get("met_qcd_ff_"+jet);
  qcdBins.clear();
  getBins(qcd,qcdBins);
  std::cout << "qcd events ----------" << std::endl;
  printBins(qcdBins);

  TH1F* ew = (TH1F*) f->Get("met_eg_"+jet);
  ewBins.clear();
  getBins(ew,ewBins);
  std::cout << "ew  events ----------" << std::endl;
  printBins(ewBins);

  TH1F* qcd_ee = (TH1F*) f->Get("met_qcd_ee_"+jet);
  std::vector<BinInfo> qcd_eeBins;
  getBins(qcd_ee,qcd_eeBins);
  std::cout << "qcd_ee events ----------" << std::endl;
  printBins(qcd_eeBins);

  int N = int(qcd_eeBins.size());
  for(int i=0; i<N; i++){
    double diff = qcd_eeBins[i].y - qcdBins[i].y;
    BinInfo bin;
    bin.x = qcdBins[i].x;
    bin.y = fabs(diff);
    qcdSysErrors.push_back(bin);
  }

  std::cout << "qcd_sysErrors ---------" << std::endl;
  printBins(qcdSysErrors);

}


void readSig(TFile* f, TString jet, int mS, int mG, int mN,
	     std::vector<BinInfo>& sig_ggBins,
	     std::vector<BinInfo>& sig_ffBins,
	     std::vector<BinInfo>& sigBins) {

  std::stringstream ggname;
  ggname << "gg_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;
  TH1F* gg = (TH1F*) f->Get(ggname.str().c_str());
  if(!gg) {
    std::cout << "histogram " << ggname.str() << " is not available!!!" << std::endl;
    return;
  }
  sig_ggBins.clear();
  getBins(gg,sig_ggBins);

  std::stringstream ffname;
  ffname << "ff_met_" << jet << "_mS" << mS << "_mG" << mG << "_mN" << mN;
  TH1F* ff = (TH1F*) f->Get(ffname.str().c_str());
  sig_ffBins.clear();
  getBins(ff,sig_ffBins);

  sigBins.clear();
  int n = sig_ggBins.size();
  for(int k=0; k<n; k++) {
    BinInfo b;
    b.x = sig_ggBins[k].x;
    b.y = sig_ggBins[k].y - sig_ffBins[k].y;
    if(b.y < 1e-6) b.y = 0.0;
    b.error = std::sqrt(sig_ggBins[k].error * sig_ggBins[k].error + sig_ffBins[k].error * sig_ffBins[k].error);
    sigBins.push_back(b);
  }//k
  //  printBins(sigBins);

}


void getBins(TH1F* h, std::vector<BinInfo>& binInfos){

  for(int i=0; i<NCH; i++) {
    int ibin = h->GetXaxis()->FindBin(bins[i]);
    int jbin = -1;
    if(i<NCH-1) jbin = h->GetXaxis()->FindBin(bins[i+1]) - 1;
    BinInfo binInfo;
    binInfo.x = bins[i];
    //    binInfo.y = h->IntegralAndError(ibin,jbin,binInfo.error,"width");
    binInfo.y = h->IntegralAndError(ibin,jbin,binInfo.error);
    binInfos.push_back(binInfo);
  }

}

void printBins(std::vector<BinInfo>& binInfos){

  for(std::vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    printf("%3d(%6.2f +/- %5.2f) ",it->x,it->y,it->error);
  }
  printf("\n");

}

void printBins(std::fstream& of, std::vector<BinInfo>& binInfos) {
  for(std::vector<BinInfo>::iterator it = binInfos.begin();
      it != binInfos.end(); it++) {
    of << "(" << it->x << ", " << it->y << " +/- " << it->error << ") ";
  }
  of << std::endl;
}


void GetNGen(std::vector<GridPoint>& grid, TString datafile) {
  unsigned int ngrid = grid.size();
  std::ifstream fin;
  fin.open(datafile.Data());
  while(1){
    int mG, mS, mChi, nGen, nAcc;
    double dummy;
    fin >> mG >> mS >> mChi >> nGen >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<ngrid; ig++){
      if( (mG ==grid[ig].mG) && (mS ==grid[ig].mS) ){
	grid[ig].ngen = nGen;
      }
    }
  }
  fin.close();
}


void GetXSection(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
  // RS errors are absolute.
  fin.open(datafile.Data());
  while(1){
    int mG, mS, mN;
    double xsec, rsp, rsm, dummy;
    std::string dummyS;
    // printf("%6.0f %6.0f %6.0f LO: %9.3e + %9.3e - %9.3e NLO: %9.3e + %9.3e - %9.3e\n", $ms, $mg, $mchi, $loxsec{$ms}{$mg}{$mchi},$lohi,$lolo,$nloxsec{$ms}{$mg}{$mchi},$nlohi,$nlolo);
    fin >> mS >> mG >> mN >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> dummy >> dummyS >> xsec >> dummyS >> rsp >> dummyS >> rsm;
    if(!fin.good()) break;
    for(unsigned int ig=0; ig<ngrid; ig++){
      if( (mG ==grid[ig].mG) && (mS ==grid[ig].mS) ){
	grid[ig].xsecValue = std::abs(xsec);
	grid[ig].xsecRSErrorPos = std::abs(rsp/xsec);  // convert to relative error
	grid[ig].xsecRSErrorNeg = std::abs(rsm/xsec);  // convert to relative error
      }
    }
  }
  fin.close();
}


void GetPDFErrorsXSection(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
  fin.open(datafile.Data());
  if(datafile.Contains("mNScan")) {
    while(1){
      int mG, mS, mN;
      double exs;
      fin >> mG >> mS >> mN >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) && (mN == grid[ig].mN) ) {
	  // errors are in %. convert to relative error
	  grid[ig].xsecPDFError = 0.01*exs;
	}
      } 
    }
  }
  else {
    while(1){
      int mG, mS;
      double exs;
      fin >> mG >> mS >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) ) {
	  // errors are in %. convert to relative error
	  grid[ig].xsecPDFError = 0.01*exs;
	}
      } 
    }
  }
  fin.close();
}


void GetPDFErrorsAcceptance(std::vector<GridPoint>& grid, TString datafile)
{
  unsigned int ngrid = grid.size();
  std::ifstream fin;
  fin.open(datafile.Data());
  if(datafile.Contains("mNScan")) {
    while(1){
      int mG, mS, mN;
      double exs;
      fin >> mG >> mS >> mN >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) && (mN == grid[ig].mN) ) {
	  // errors are in %. convert to relative error
	  grid[ig].accErrorPDF = 0.01*exs;
	}
      } 
    }
  }
  else {
    while(1){
      int mG, mS;
      double exs;
      fin >> mG >> mS >> exs;
      if(!fin.good()) break;
      for(unsigned int ig=0; ig<ngrid; ig++){
	if( (mG == grid[ig].mG) && (mS == grid[ig].mS) ) {
	  // errors are in %. convert to relative error
	  grid[ig].accErrorPDF = 0.01*exs;
	}
      } 
    }
  }
  fin.close();
}


void makeSignalGains(std::vector<GridPoint>& grid) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {
    double nevt = 0;
    for(std::vector<BinInfo>::iterator bit = it->sigBins.begin();
	bit != it->sigBins.end(); bit++) {
      nevt += bit->y;
      double acc = bit->y * DataMCScale * DataMCScale / it->ngen;
      double nexpected = it->xsecValue * it->lumi * acc;
      double error = bit->error * it->xsecValue * it->lumi * DataMCScale * DataMCScale / it->ngen;
      bit->y = nexpected;
      bit->error = error;
      //      std::cout << "y, xsec, acc : " << bit->y << ", " << it->xsecValue << ", " << acc << std::endl;
    }// bit
    it->sig_sysError = std::sqrt(it->sig_sysError*it->sig_sysError + it->accErrorPDF*it->accErrorPDF);
    it->acc = nevt/it->ngen;
    //    printBins(it->sigBins);
  }// it

}


void makeDataCard(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    if(it->sigBins.size() == 0) continue;

    std::stringstream outname;
    outname << datacard_dir.Data() << "/multiChannel/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
    std::fstream outfile(outname.str().c_str(),std::ios::out);

    outfile << "# mG : " << it->mG << std::endl;
    outfile << "# mS : " << it->mS << std::endl;
    outfile << "# mN : " << it->mN << std::endl;
    outfile << "# ngen : " << it->ngen << std::endl;
    outfile << "# acc : " << it->acc << std::endl;
    outfile << "# lumi : " << it->lumi << std::endl;
    outfile << "# xsecValue : " << it->xsecValue << std::endl;
    outfile << "# xsecPDFError : " << it->xsecPDFError << std::endl;
    outfile << "# xsecRSErrorNeg : " << it->xsecRSErrorNeg << std::endl;
    outfile << "# xsecRSErrorPos : " << it->xsecRSErrorPos << std::endl;
    outfile << "# accErrorPDF : " << it->accErrorPDF << std::endl;
    outfile << "# lumi_sysError : " << it->lumi_sysError << std::endl;
    outfile << "# qcd_sysError : " << it->qcd_sysError << std::endl;
    outfile << "# ew_sysError : " << it->ew_sysError << std::endl;
    outfile << "# sig_sysError : " << it->sig_sysError << std::endl;
    outfile << "# ggBins : "; printBins(outfile,it->ggBins);
    outfile << "# ewBins : "; printBins(outfile,it->ewBins);
    outfile << "# qcdBins : "; printBins(outfile,it->qcdBins);
    outfile << "# qcdSysErrors : "; printBins(outfile,it->qcdSysErrors);
    outfile << "# sig_ggBins : "; printBins(outfile,it->sig_ggBins);
    outfile << "# sig_ffBins : "; printBins(outfile,it->sig_ffBins);
    outfile << "# sigBins : "; printBins(outfile,it->sigBins);

    std::vector<int> sensitive_bins;
    for(int i=0; i<it->sigBins.size(); i++){
      if(it->sigBins[i].y > epsilon) sensitive_bins.push_back(i);
    }

    int nch = sensitive_bins.size();

    outfile << "imax " << nch << " number of channels" << std::endl;
    outfile << "jmax 2 number of backgrounds" << std::endl;
    outfile << "kmax " << (5 + nch*3) << " number of nuisance parameters" << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i;
    outfile << std::endl;

    outfile << "observation ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->ggBins[bin].y;
    }
    outfile << std::endl;

    outfile << "--------------" << std::endl;
    outfile << "bin         ";
    for(int i=0; i<nch; i++) outfile << "\t" << i << " " << i << " " << i;
    outfile << std::endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t susy qcd ew";
    outfile << std::endl;

    outfile << "process       ";
    for(int i=0; i<nch; i++) outfile << "\t 0 1 2";
    outfile << std::endl;

    outfile << "rate          ";
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "\t" << it->sigBins[bin].y << " " << it->qcdBins[bin].y << " " << it->ewBins[bin].y;
    }
    outfile << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "u_lumi lnN    ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->lumi_sysError << " - - ";
    outfile << std::endl;

    outfile << "u_id lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << it->sig_sysError << " - - ";
    outfile << std::endl;

    outfile << "u_JES lnN     ";
    for(int i=0; i<nch; i++) outfile << "\t" << "1.02 - - ";
    outfile << std::endl;

    outfile << "u_qcd lnN     ";
    for(int i=0; i<nch; i++){
      int bin = sensitive_bins[i];
      outfile << "\t" << "- " << 1+(it->qcdSysErrors[bin].y/it->qcdBins[bin].y) << " - ";
    }
    outfile << std::endl;
    
    outfile << "u_ew lnN      ";
    for(int i=0; i<nch; i++) outfile << "\t" << "- - " << it->ew_sysError;
    outfile << std::endl;
    
    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_susy_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << 1+(it->sigBins[bin].error/it->sigBins[bin].y) << " - - ";
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_qcd_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - " << 1 + (it->qcdBins[bin].error/it->qcdBins[bin].y) << " - ";
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

    for(int i=0; i<nch; i++) {
      int bin = sensitive_bins[i];
      outfile << "stat_ew_bin" << i << " lnN ";
      for(int j=0; j<nch; j++) {
	if(i==j) outfile << " - - " << 1 + (it->ewBins[bin].error/it->ewBins[bin].y);
	else outfile << " - - - ";
      }//j
      outfile << std::endl;
    }// i

  }// it


}



void makeDataCardSingleChannel(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    int nbins = int(it->sigBins.size());
    if(nbins == 0) continue;

    double nData = 0;
    double nQcd = 0;
    double nEw = 0;
    double nSig = 0;
    double qcd_err = 0;
    double ew_err = 0;
    double sig_err = 0;
    double qcdSysError = 0;

    for(int i=2; i<nbins; i++) { // starting from 100 GeV bin
      nData += it->ggBins[i].y;
      nQcd += it->qcdBins[i].y;
      nEw += it->ewBins[i].y;
      nSig += it->sigBins[i].y;
      qcd_err += it->qcdBins[i].error * it->qcdBins[i].error;
      ew_err += it->ewBins[i].error * it->ewBins[i].error;
      sig_err += it->sigBins[i].error * it->sigBins[i].error;
      qcdSysError += it->qcdSysErrors[i].y;
    }
    qcd_err = std::sqrt(qcd_err);
    ew_err = std::sqrt(ew_err);
    sig_err = std::sqrt(sig_err);
    qcdSysError = 1 + qcdSysError/nQcd;

    if(nSig < epsilon) continue;

    std::stringstream outname;
    outname << datacard_dir.Data() << "/singleChannel/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
    std::fstream outfile(outname.str().c_str(),std::ios::out);

    outfile << "# mG : " << it->mG << std::endl;
    outfile << "# mS : " << it->mS << std::endl;
    outfile << "# mN : " << it->mN << std::endl;
    outfile << "# ngen : " << it->ngen << std::endl;
    outfile << "# acc : " << it->acc << std::endl;
    outfile << "# lumi : " << it->lumi << std::endl;
    outfile << "# xsecValue : " << it->xsecValue << std::endl;
    outfile << "# xsecPDFError : " << it->xsecPDFError << std::endl;
    outfile << "# xsecRSErrorNeg : " << it->xsecRSErrorNeg << std::endl;
    outfile << "# xsecRSErrorPos : " << it->xsecRSErrorPos << std::endl;
    outfile << "# accErrorPDF : " << it->accErrorPDF << std::endl;
    outfile << "# lumi_sysError : " << it->lumi_sysError << std::endl;
    outfile << "# qcd_sysError : " << it->qcd_sysError << std::endl;
    outfile << "# ew_sysError : " << it->ew_sysError << std::endl;
    outfile << "# sig_sysError : " << it->sig_sysError << std::endl;
    outfile << "# ggBins : "; printBins(outfile,it->ggBins);
    outfile << "# ewBins : "; printBins(outfile,it->ewBins);
    outfile << "# qcdBins : "; printBins(outfile,it->qcdBins);
    outfile << "# qcdSysErrors : "; printBins(outfile,it->qcdSysErrors);
    outfile << "# sig_ggBins : "; printBins(outfile,it->sig_ggBins);
    outfile << "# sig_ffBins : "; printBins(outfile,it->sig_ffBins);
    outfile << "# sigBins : "; printBins(outfile,it->sigBins);

    outfile << "imax 1 number of channels" << std::endl;
    outfile << "jmax 2 number of backgrounds" << std::endl;
    outfile << "kmax 8 number of nuisance parameters" << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "bin 0" << std::endl;
    outfile << "observation " << nData << std::endl;
    outfile << "--------------" << std::endl;
    outfile << "bin 0 0 0" << std::endl;
    outfile << "process susy qcd ew" << std::endl;
    outfile << "process 0 1 2" << std::endl;
    outfile << "rate " << nSig << " " << nQcd << " " << nEw << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "u_lumi lnN    " << it->lumi_sysError << " - - " << std::endl;
    outfile << "u_id lnN      " << it->sig_sysError << " - - " << std::endl;
    outfile << "u_JES lnN     1.02 - - " << std::endl;
    outfile << "u_qcd lnN     - " << qcdSysError << " - " << std::endl;
    outfile << "u_ew lnN      - - " << it->ew_sysError << std::endl;
    
    outfile << "stat_susy lnN " << (1 + sig_err/nSig) << " - -" << std::endl;
    outfile << "stat_qcd  lnN - " << (1 + qcd_err/nQcd) << " -" << std::endl;
    outfile << "stat_ew   lnN - - " << (1 + ew_err/nEw) << std::endl;

  }// it

}



// The numbers in the histogram do not match with event counts. So here I put the numbers manually
// by separating QCD/EW and stat/sys errors.
// This should be fixed when the proper histograms are provided.
/*
void makeDataCardSingleChannel(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    int nbins = int(it->sigBins.size());
    if(nbins == 0) continue;

    double nData = 11;
    double nQcd = 11.392;
    double nEw = 2.8985;
    double qcd_err = 2.933;
    double ew_err = 0.2065;

    if(jet.Contains("nojet")){
      nData = 17;
      nQcd = 20.188;
      nEw = 3.4724;
      qcd_err = 6.659;
      ew_err = 0.226;
    }

    double nSig = 0;
    double sig_err = 0;
    for(int i=2; i<nbins; i++) { // starting from 100 GeV bin
      nSig += it->sigBins[i].y;
      sig_err += it->sigBins[i].error * it->sigBins[i].error;
    }
    sig_err = std::sqrt(sig_err);

    if(nSig < epsilon) continue;

    std::stringstream outname;
    outname << datacard_dir.Data() << "/singleChannel/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
    std::fstream outfile(outname.str().c_str(),std::ios::out);

    outfile << "# mG : " << it->mG << std::endl;
    outfile << "# mS : " << it->mS << std::endl;
    outfile << "# mN : " << it->mN << std::endl;
    outfile << "# ngen : " << it->ngen << std::endl;
    outfile << "# acc : " << it->acc << std::endl;
    outfile << "# lumi : " << it->lumi << std::endl;
    outfile << "# xsecValue : " << it->xsecValue << std::endl;
    outfile << "# xsecPDFError : " << it->xsecPDFError << std::endl;
    outfile << "# xsecRSErrorNeg : " << it->xsecRSErrorNeg << std::endl;
    outfile << "# xsecRSErrorPos : " << it->xsecRSErrorPos << std::endl;
    outfile << "# accErrorPDF : " << it->accErrorPDF << std::endl;
    outfile << "# lumi_sysError : " << it->lumi_sysError << std::endl;
    outfile << "# qcd_sysError : " << it->qcd_sysError << std::endl;
    outfile << "# ew_sysError : " << it->ew_sysError << std::endl;
    outfile << "# sig_sysError : " << it->sig_sysError << std::endl;
    outfile << "# ggBins : "; printBins(outfile,it->ggBins);
    outfile << "# ewBins : "; printBins(outfile,it->ewBins);
    outfile << "# qcdBins : "; printBins(outfile,it->qcdBins);
    outfile << "# sig_ggBins : "; printBins(outfile,it->sig_ggBins);
    outfile << "# sig_ffBins : "; printBins(outfile,it->sig_ffBins);
    outfile << "# sigBins : "; printBins(outfile,it->sigBins);

    outfile << "imax 1 number of channels" << std::endl;
    outfile << "jmax 2 number of backgrounds" << std::endl;
    outfile << "kmax 8 number of nuisance parameters" << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "bin 0" << std::endl;
    outfile << "observation " << nData << std::endl;
    outfile << "--------------" << std::endl;
    outfile << "bin 0 0 0" << std::endl;
    outfile << "process susy qcd ew" << std::endl;
    outfile << "process 0 1 2" << std::endl;
    outfile << "rate " << nSig << " " << nQcd << " " << nEw << std::endl;
    outfile << "--------------" << std::endl;

    outfile << "u_lumi lnN    " << it->lumi_sysError << " - - " << std::endl;
    outfile << "u_id lnN      " << it->sig_sysError << " - - " << std::endl;
    outfile << "u_JES lnN     1.02 - - " << std::endl;
    outfile << "u_qcd lnN     - " << it->qcd_sysError << " - " << std::endl;
    outfile << "u_ew lnN      - - " << it->ew_sysError << std::endl;
    
    outfile << "stat_susy lnN " << (1 + sig_err/nSig) << " - -" << std::endl;
    outfile << "stat_qcd  lnN - " << (1 + qcd_err/nQcd) << " -" << std::endl;
    outfile << "stat_ew   lnN - - " << (1 + ew_err/nEw) << std::endl;

  }// it

}
*/


void makeDataCardEachChannel(std::vector<GridPoint>& grid, TString bino, TString jet) {

  for(std::vector<GridPoint>::iterator it = grid.begin();
      it != grid.end(); it++) {

    int nbins = int(it->sigBins.size());
    if(nbins == 0) continue;

    for(int i=0; i<nbins; i++) {

      double nData = it->ggBins[i].y;
      double nQcd = it->qcdBins[i].y;
      double nEw = it->ewBins[i].y;
      double nSig = it->sigBins[i].y;
      double qcd_err = it->qcdBins[i].error;
      double ew_err = it->ewBins[i].error;
      double sig_err = it->sigBins[i].error;

      double qcdSysError = 1 + it->qcdSysErrors[i].y/it->qcdBins[i].y;

      if(nSig < epsilon) continue;

      std::stringstream outname;
      outname << datacard_dir.Data() << "/bin" << i << "/" << bino.Data() << "_mS" << it->mS << "_mG" << it->mG << "_mN" << it->mN << "_" << jet << ".dat";
      std::fstream outfile(outname.str().c_str(),std::ios::out);

      outfile << "# mG : " << it->mG << std::endl;
      outfile << "# mS : " << it->mS << std::endl;
      outfile << "# mN : " << it->mN << std::endl;
      outfile << "# ngen : " << it->ngen << std::endl;
      outfile << "# acc : " << it->acc << std::endl;
      outfile << "# lumi : " << it->lumi << std::endl;
      outfile << "# xsecValue : " << it->xsecValue << std::endl;
      outfile << "# xsecPDFError : " << it->xsecPDFError << std::endl;
      outfile << "# xsecRSErrorNeg : " << it->xsecRSErrorNeg << std::endl;
      outfile << "# xsecRSErrorPos : " << it->xsecRSErrorPos << std::endl;
      outfile << "# accErrorPDF : " << it->accErrorPDF << std::endl;
      outfile << "# lumi_sysError : " << it->lumi_sysError << std::endl;
      outfile << "# qcd_sysError : " << it->qcd_sysError << std::endl;
      outfile << "# ew_sysError : " << it->ew_sysError << std::endl;
      outfile << "# sig_sysError : " << it->sig_sysError << std::endl;
      outfile << "# ggBins : "; printBins(outfile,it->ggBins);
      outfile << "# ewBins : "; printBins(outfile,it->ewBins);
      outfile << "# qcdBins : "; printBins(outfile,it->qcdBins);
      outfile << "# qcdSysErrors : "; printBins(outfile,it->qcdSysErrors);
      outfile << "# sig_ggBins : "; printBins(outfile,it->sig_ggBins);
      outfile << "# sig_ffBins : "; printBins(outfile,it->sig_ffBins);
      outfile << "# sigBins : "; printBins(outfile,it->sigBins);

      outfile << "imax 1 number of channels" << std::endl;
      outfile << "jmax 2 number of backgrounds" << std::endl;
      outfile << "kmax 8 number of nuisance parameters" << std::endl;
      outfile << "--------------" << std::endl;

      outfile << "bin " << i << std::endl;
      outfile << "observation " << nData << std::endl;
      outfile << "--------------" << std::endl;
      outfile << "bin " << i << " " << i << " " << i << std::endl;
      outfile << "process susy qcd ew" << std::endl;
      outfile << "process 0 1 2" << std::endl;
      outfile << "rate " << nSig << " " << nQcd << " " << nEw << std::endl;
      outfile << "--------------" << std::endl;
      
      outfile << "u_lumi lnN    " << it->lumi_sysError << " - - " << std::endl;
      outfile << "u_id lnN      " << it->sig_sysError << " - - " << std::endl;
      outfile << "u_JES lnN     1.02 - - " << std::endl;
      outfile << "u_qcd lnN     - " << qcdSysError << " - " << std::endl;
      outfile << "u_ew lnN      - - " << it->ew_sysError << std::endl;
      
      outfile << "stat_susy lnN " << (1 + sig_err/nSig) << " - -" << std::endl;
      outfile << "stat_qcd  lnN - " << (1 + qcd_err/nQcd) << " -" << std::endl;
      outfile << "stat_ew   lnN - - " << (1 + ew_err/nEw) << std::endl;

    }// for i
  }// it


}
