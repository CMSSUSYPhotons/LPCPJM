#include <vector>
#include <iostream>
#include <fstream>

#include <TSystem.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>

#include "GridPoint.C"
#include "roostats_cl95.C"

void ra3(int mG = 1040, int mS = 1040, int mN = 375, int nuisanceModel = 1, TString jetCut="1jet", TString metCut="met100", TString bino="bino") {

  // mG : 400 - 2000 GeV incremented by 80 GeV step
  // mS : 400 - 2000 GeV incremented by 80 GeV step
  // mN : 375, 150 only two points
  // nuisanceModel : 0(gauss), 1(lognormal), 2(gamma) : lognormal is recommended by stat committee
  // jetCut : 1jet, nojet
  // metCut : met50, met70, met100
  // bino : bino, wino

  TStopwatch ts;
  ts.Start();

  std::vector<GridPoint> grid;
  grid.clear();

  char datafile[100];
  sprintf(datafile,"mcAccMap_%s_mN%d.dat",bino.Data(),mN);

  CreateGrid(grid,TString(datafile),metCut,jetCut);              // signal acceptance

  sprintf(datafile,"%s%dNLOxsec.dat",bino.Data(),mN);
  GetXSection(grid,TString(datafile));             // XS and RS errors
  GetPDFErrorsXSection(grid,"xsectionPDFErrors.dat");    // PDF errors on XS
  GetPDFErrorsAcceptance(grid,"acceptancePDFErrors.dat");  // PDF errors on acceptance
  CalculateAcceptanceCorrectedJuly2011(grid); // Data/MC scale, JES scale, etc.

  unsigned int nUL = grid.size();

  // Met 100 GeV version 1+jet
  int n = 0;                   // candidates
  double ilum = 1143;          // int. luminosity
  double slum = 0.06*ilum;     // 6% on lumi
  double eff;                  // taken from the text file
  double seff;                 // taken from the grid points
  double bck  = 1.45821;       // background (ee, ff combined)
  double sbck = 0.977438;       // absolute error on background
  bool gauss = false;
  int seed = 89275;
  std::string prob = "cls";

  // Met 100 GeV version nojet
  if(metCut.Contains("met100") && jetCut.Contains("nojet")) {
    n = 1;
    bck  = 2.15612;
    sbck = 0.666149;
  }
  else if(metCut.Contains("met100") && jetCut.Contains("ff")) {
    n = 0;
    bck  = 2.52149;
    sbck = 2.19458;
  }
  else if(metCut.Contains("met100") && jetCut.Contains("ee")) {
    n = 0;
    bck  = 1.31018;
    sbck = 0.818844;
  }
  else if(metCut.Contains("met50") && jetCut.Contains("1jet")) {
    n = 9;
    bck  = 11.3275;
    sbck = 2.06947;
  }
  else if(metCut.Contains("met50") && jetCut.Contains("nojet")) {
    n = 12;
    bck  = 17.6904;
    sbck = 1.94204;
  }


  // ------------ Upper limit in One Grid Point -------
  GridPoint* gridPoint = 0;
  char tree_name[200];
  sprintf(tree_name,"tree_%s_%s_%s_mG%d_mS%d_mN%d_model%d.root",bino.Data(),metCut.Data(),jetCut.Data(),mG,mS,mN,nuisanceModel);

  TFile* f = new TFile(tree_name,"RECREATE");
  TTree* tree = new TTree("gridPoints","GridPoint");
  tree->Branch("GridPoint",&gridPoint,32000,99);

  for(unsigned int i=0; i<nUL; i++) {

    if( grid[i].mGluino == mG && grid[i].mSquark == mS && grid[i].mNeutralino == mN ) {

      eff  = grid[i].accValue;
      seff = grid[i].accErrorTot;
      seed = int(eff*seff*1e8);

      std::cout << "Input parameters =====>" << std::endl;
      std::cout << "ilum : " << ilum << " +/- " << slum << std::endl;
      std::cout << "eff : " << eff << " +/- " << seff << std::endl;
      std::cout << "bck : " << bck << " +/- " << sbck << std::endl;
      std::cout << "n : " << n << std::endl;
      std::cout << "gauss : " << gauss << std::endl;
      std::cout << "seed : " << seed << std::endl;

      LimitResult limit = roostats_limit(ilum, slum, eff, seff, bck, sbck, n, gauss, nuisanceModel, prob, "test.gif", seed);

      std::cout << "ObservedLimit     : " << limit.GetObservedLimit() << std::endl;
      std::cout << "ExpectedLimit     : " << limit.GetExpectedLimit() << std::endl;
      std::cout << "OneSigmaLowRange  : " << limit.GetOneSigmaLowRange() << std::endl;
      std::cout << "OneSigmaHighRange : " << limit.GetOneSigmaHighRange() << std::endl;
      std::cout << "TwoSigmaLowRange  : " << limit.GetTwoSigmaLowRange() << std::endl;
      std::cout << "TwoSigmaHighRange : " << limit.GetTwoSigmaHighRange() << std::endl;

      grid[i].nuisanceModel = nuisanceModel;
      grid[i].limit       = limit.GetObservedLimit();
      grid[i].explimit    = limit.GetExpectedLimit();
      grid[i].explimit_1L = limit.GetOneSigmaLowRange();
      grid[i].explimit_1H = limit.GetOneSigmaHighRange();
      grid[i].explimit_2L = limit.GetTwoSigmaLowRange();
      grid[i].explimit_2H = limit.GetTwoSigmaHighRange();

      gridPoint = &grid[i];
      tree->Fill();

      grid[i].Print();

    } // if
  } // for

  f->cd();
  tree->Write();
  f->Close();
  

  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
