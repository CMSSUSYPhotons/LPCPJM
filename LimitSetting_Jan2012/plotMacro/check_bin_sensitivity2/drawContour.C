
#include "util.C"

void drawContour(TString bino="bino", TString jet="nojet", TString channel = "multiChannel", bool print=false) {

  bool useCustomGetContour = false;

  TString data_dir = "/uscms/home/dwjang/rel/428/src/ra3/plotMacro/20120131/"+channel;

  gStyle->SetOptStat(0);

  TString label = bino + "_mN375_met100_" + jet + "_" + channel;
  if(bino.Contains("mNScan")) label = bino + "_met100_" + jet + "_" + channel;

  gStyle->SetPalette(1);
  gStyle->SetPadLeftMargin(0.15);

  // for bino & wino
  const int nG = 21;
  const int nS = 21;
  float mGVals[nG+1] = {400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};
  float mSVals[nS+1] = {400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};

  // to make valuse on the center of each bin
  float* mGBins = new float[nG+1];
  float* mSBins = new float[nS+1];
  for(int i=0; i<nG+1; i++) mGBins[i] = mGVals[i]-40;
  for(int i=0; i<nS+1; i++) mSBins[i] = mSVals[i]-40;


  // for mNScan
  const int nX = 10;
  const int nY = 24;
  float xVals[nX+1] = {150, 250, 350, 450, 550, 650, 750, 850, 950, 1050, 1150};
  float yVals[nY+1] = {160, 240, 320, 400, 480, 560, 640, 720, 800, 880, 960, 1040, 1120, 1200, 1280, 1360, 1440, 1520, 1600, 1680, 1760, 1840, 1920, 2000, 2100};

  float* xBins = new float[nX+1];
  float* yBins = new float[nY+1];

  for(int i=0; i<nX+1; i++) xBins[i] = xVals[i]-50;
  for(int i=0; i<nY+1; i++) yBins[i] = yVals[i]-40;


  TFile* fout = new TFile("hist_exclusion_"+label+".root","RECREATE");

  const int nxs = 6;
  TString xsname[nxs] = {"xsec","xsec_1L","xsec_1H","xsec_2L","xsec_2H","acc"};
  TH2D* h_xs[nxs];
  for(int i=0; i<nxs; i++) {
    if(bino.Contains("mNScan")) h_xs[i] = new TH2D(xsname[i],xsname[i],nX,xBins,nY,yBins);
    else h_xs[i] = new TH2D(xsname[i],xsname[i],nS,mSBins,nG,mGBins);
  }

  const int nlimit = 6;
  TString limitname[nlimit] = {"limit", "exp", "exp_1L","exp_1H","exp_2L","exp_2H"};
  TH2D* h_limit[nlimit];
  for(int i=0; i<nlimit; i++) {
    if(bino.Contains("mNScan")) h_limit[i] = new TH2D(limitname[i],limitname[i],nX,xBins,nY,yBins);
    else h_limit[i] = new TH2D(limitname[i],limitname[i],nS,mSBins,nG,mGBins);
  }

  TString datafile = data_dir + "/" + bino + "_" + jet + ".table";

  std::ifstream fin;
  fin.open(datafile.Data());
  while(1){
    // #echo "mS mG mN acc xsec xsecPDFError xsecRSErrorNeg xsecRSErrorPos obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s"
    int mS, mG, mN;
    double acc, xsec, xsecPDFError, xsecRSErrorNeg, xsecRSErrorPos, obsLimit, expLimit, exp_m1s, exp_m2s, exp_p1s, exp_p2s;
    fin >> mS >> mG >> mN >> acc >> xsec >> xsecPDFError >> xsecRSErrorNeg >> xsecRSErrorPos >> obsLimit >> expLimit >> exp_m1s >> exp_m2s >> exp_p1s >> exp_p2s;
    if(!fin.good()) break;

//     std::cout << mS << ", " << mG << ", " << mN << ", " << acc << ", " << xsec << ", " << xsecPDFError << ", "
// 	      << xsecRSErrorNeg << ", " << xsecRSErrorPos << ", " << obsLimit << ", " << expLimit << ", "
// 	      << exp_m1s << ", " << exp_m2s << ", " << exp_p1s << ", " << exp_p2s << std::endl;

    double oneSigma_L = std::sqrt(xsecRSErrorNeg * xsecRSErrorNeg + xsecPDFError * xsecPDFError);
    double oneSigma_H = std::sqrt(xsecRSErrorPos * xsecRSErrorPos + xsecPDFError * xsecPDFError);

    if(bino.Contains("mNScan")) {
      if(mS != 2500) continue;
      h_xs[5]->Fill(mN,mG,acc);
      h_xs[0]->Fill(mN,mG,xsec);
      h_xs[1]->Fill(mN,mG,xsec - xsec*oneSigma_L);
      h_xs[2]->Fill(mN,mG,xsec + xsec*oneSigma_H);
      h_xs[3]->Fill(mN,mG,xsec - xsec*2*oneSigma_L);
      h_xs[4]->Fill(mN,mG,xsec + xsec*2*oneSigma_H);

      h_limit[0]->Fill(mN,mG,obsLimit*xsec);
      h_limit[1]->Fill(mN,mG,expLimit*xsec);
      h_limit[2]->Fill(mN,mG,exp_m1s*xsec);
      h_limit[3]->Fill(mN,mG,exp_p1s*xsec);
      h_limit[4]->Fill(mN,mG,exp_m2s*xsec);
      h_limit[5]->Fill(mN,mG,exp_p2s*xsec);
    }
    else {
      if(mN != 375) continue;
      h_xs[5]->Fill(mS,mG,acc);
      h_xs[0]->Fill(mS,mG,xsec);
      h_xs[1]->Fill(mS,mG,xsec - xsec*oneSigma_L);
      h_xs[2]->Fill(mS,mG,xsec + xsec*oneSigma_H);
      h_xs[3]->Fill(mS,mG,xsec - xsec*2*oneSigma_L);
      h_xs[4]->Fill(mS,mG,xsec + xsec*2*oneSigma_H);

      h_limit[0]->Fill(mS,mG,obsLimit*xsec);
      h_limit[1]->Fill(mS,mG,expLimit*xsec);
      h_limit[2]->Fill(mS,mG,exp_m1s*xsec);
      h_limit[3]->Fill(mS,mG,exp_p1s*xsec);
      h_limit[4]->Fill(mS,mG,exp_m2s*xsec);
      h_limit[5]->Fill(mS,mG,exp_p2s*xsec);
    }// if - else

  }// while
  fin.close();

  for(int i=0; i<nxs; i++) fillPotHoles(h_xs[i]);
  for(int i=0; i<nlimit; i++) fillPotHoles(h_limit[i]);

  TString option2D = "COL Z";
  //  TString option2D = "TEXT";

  TCanvas* can_acc = new TCanvas("can_acc_"+label,"can_acc_"+label,1000,800);
  can_acc->SetRightMargin(0.15);
  h_xs[5]->GetXaxis()->SetNdivisions(505);
  h_xs[5]->GetYaxis()->SetNdivisions(505);
  h_xs[5]->Draw(option2D);
  if(print) {
    can_acc->Print("",".gif");
    can_acc->Print("",".pdf");
  }

  TCanvas* can_xs = new TCanvas("can_xsec_"+label,"can_xsec_"+label,1000,800);
  can_xs->SetRightMargin(0.15);
  h_xs[0]->GetXaxis()->SetNdivisions(505);
  h_xs[0]->GetYaxis()->SetNdivisions(505);
  h_xs[0]->Draw(option2D);
  if(print) {
    can_xs->Print("",".gif");
    can_xs->Print("",".pdf");
  }

  TCanvas* can_limit = new TCanvas("can_limit_"+label,"can_limit_"+label,1000,800);
  can_limit->SetRightMargin(0.15);
  h_limit[0]->GetXaxis()->SetNdivisions(505);
  h_limit[0]->GetYaxis()->SetNdivisions(505);
  h_limit[0]->Draw(option2D);
  if(print) {
    can_limit->Print("",".gif");
    can_limit->Print("",".pdf");
  }

  // now find exclusion curves
  TCanvas* can_diff = new TCanvas("can_diff_"+label,"can_diff_"+label,1200,800);
  //  can_diff->Divide(nlimit,nxs);
  TH2D* h_excl[nlimit][nxs-1];
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs-1; j++) {

      h_excl[i][j] = (TH2D*) h_limit[i]->Clone("exclusion_"+limitname[i]+"_"+xsname[j]);
      int nbinsx = h_excl[i][j]->GetXaxis()->GetNbins();
      int nbinsy = h_excl[i][j]->GetYaxis()->GetNbins();

      for( int ibx=1; ibx<=nbinsx; ++ibx){
	for( int iby=1; iby<=nbinsy; ++iby){
	  double x1 = h_limit[i]->GetBinContent(ibx,iby);
	  double x2 = h_xs[j]->GetBinContent(ibx,iby);
	  h_excl[i][j]->SetBinContent(ibx,iby,x2-x1);
	  x1 = h_limit[i]->GetBinError(ibx,iby);
	  x2 = h_xs[j]->GetBinError(ibx,iby);
	  h_excl[i][j]->SetBinError(ibx,iby,std::sqrt(x1*x1+x2*x2));
	}// for iby
      }// for ibx
      fixBadCells(h_excl[i][j]);
      if(i==0 && j==0) h_excl[i][j]->Draw("TEXT");
    }// for j
  }// for i


  float xmin = 400;
  float ymin = 400;
  float xmax = 2000;
  float ymax = 2000;
  if(bino.Contains("mNScan")){
    xmin = 200;
    xmax = 1500;
    ymin = 500;
    ymax = 2000;
  }


  TGraph* curv[nlimit][nxs-1];

  TH2D* h_back;
  if(bino.Contains("mNScan")) h_back = new TH2D("h_back",";M_{#tilde{#chi^{0}}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,xmin,xmax,100,ymin,ymax);
  else h_back = new TH2D("h_back",";M_{#tilde{q}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,xmin,xmax,100,ymin,ymax);

  h_back->GetXaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetTitleOffset(1.2);


  double contours[2]={ 0.0, 1.0 };
  TCanvas *can_excl01 = new TCanvas("can_excl01_"+label, "can_excl01_"+label,1200,800);
  can_excl01->Divide(nlimit,nxs-1);
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs-1; j++) {

      can_excl01->cd(j*nlimit + i + 1);
      h_back->Draw();

      if(useCustomGetContour) {
        curv[i][j] = getContour(h_excl[i][j],"excl_curv_"+limitname[i]+"_"+xsname[j]);
        curv[i][j]->Draw("SAME L");
      }
      else {
	h_excl[i][j]->SetContour(2,contours);
	h_excl[i][j]->Draw("SAME CONT LIST");
	gPad->Update();

	TObjArray *contsM = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = (TList*)contsM->At(0);
	curv[i][j] = (TGraph*)contLevel->First()->Clone("excl_curv_"+limitname[i]+"_"+xsname[j]);
      }
      //      PrintPoints(curv[i][j]);
      //      RemovePoints(curv[i][j]);

    }// for j
  }// for i

  if(bino.Contains("mNScan")) {
    for(int i=0; i<nlimit; i++) {
      for(int j=0; j<nxs-1; j++) {
	RemovePoints(curv[i][j]);
      }
    }

    for(int i=0; i<nlimit; i++) {
      for(int j=0; j<nxs-1; j++) {
	double x,y;
	int whichone = curv[i][j]->GetN()-1;
	curv[i][j]->GetPoint(whichone-1,x,y);
	curv[i][j]->SetPoint(whichone,y,y);
      }
    }
  }


  TGraph* excludedRegion = new TGraph(curv[0][0]->GetN()+3);
  int nbins = curv[0][0]->GetN();
  for(int i=0; i<nbins; i++){
    double x,y;
    curv[0][0]->GetPoint(i,x,y);
    excludedRegion->SetPoint(i,x,y);
  }

  excludedRegion->SetPoint(nbins,xmax,ymin);
  excludedRegion->SetPoint(nbins+1,xmin,ymin);
  excludedRegion->SetPoint(nbins+2,xmin,ymax);

  TGraph* oneSigma = new TGraph(curv[0][1]->GetN()+curv[0][2]->GetN());
  TGraph* twoSigma = new TGraph(curv[0][3]->GetN()+curv[0][4]->GetN());

  for(int i=1; i<nxs-1; i++) {
    int nbins = curv[0][i]->GetN();
    for(int j=0; j<nbins; j++) {
      double x,y;
      curv[0][i]->GetPoint(j,x,y);
      if(i==1) oneSigma->SetPoint(j,x,y);
      if(i==2) oneSigma->SetPoint(curv[0][1]->GetN()+nbins-j-1,x,y);
      if(i==3) twoSigma->SetPoint(j,x,y);
      if(i==4) twoSigma->SetPoint(curv[0][3]->GetN()+nbins-j-1,x,y);
    }
  }


  TCanvas* can_excl02 = new TCanvas("can_excl02_"+label, "can_excl02_"+label,1000,800);
  h_back->Draw();
  can_excl02->SetGrid(1,1);

  // ecluded region
  excludedRegion->SetFillStyle(3004);
  excludedRegion->Draw("same f");

//   // 2 sigma band
//   twoSigma->SetFillColor(42);
//   twoSigma->Draw("same f");

  // 1 sigma band
  oneSigma->SetFillColor(42);
  oneSigma->Draw("same f");

  // expected limit
  curv[1][0]->SetLineStyle(2);
  curv[1][0]->SetLineWidth(4);
  curv[1][0]->Draw("SAME L");

  // observed limit
  curv[0][0]->SetLineWidth(4);
  curv[0][0]->Draw("SAME L");

  PrintPoints(curv[0][0]);


  float leg_xmin = 0.18;
  float leg_xmax = 0.53;
  float leg_ymin = 0.45;
  float leg_ymax = 0.7;
  if(bino.Contains("mNScan")){
    leg_xmin = 0.18;
    leg_xmax = 0.5;
    leg_ymin = 0.5;
    leg_ymax = 0.7;
  }
  else if(bino.Contains("wino")){
    leg_xmin = 0.55;
    leg_xmax = 0.85;
    leg_ymin = 0.5;
    leg_ymax = 0.7;
  }


  TLegend* leg;
  if(bino.Contains("mNScan")) leg = new TLegend(leg_xmin+0.05,leg_ymin,leg_xmax,leg_ymax,"m_{#tilde{g}} = 2500 (GeV/c^{2})");
  else leg = new TLegend(leg_xmin+0.05,leg_ymin,leg_xmax,leg_ymax,bino+"-like,  m_{#tilde{#Chi}^{0}} = 375 (GeV/c^{2})");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(curv[0][0],"Observed, NLO","L");
  leg->AddEntry(curv[1][0],"Expected, NLO","L");
  leg->AddEntry(oneSigma,"#pm1#sigma, NLO","F");
  //  leg->AddEntry(twoSigma,"#pm2#sigma, NLO","F");
  leg->Draw("same");

  TLatex* lat = new TLatex(leg_xmin,0.85,"CMS Preliminary");
  lat->SetNDC(true);
  lat->SetTextSize(0.04);
  lat->Draw("same");

  TLatex* lat2 = new TLatex(leg_xmin,0.77,"#sqrt{s} = 7 TeV,  #int#it{L}dt = 4.7 fb^{-1}");
  lat2->SetNDC(true);
  lat2->SetTextSize(0.04);
  lat2->Draw("same");

  float yv = 0.25;
  if(bino.Contains("wino")) yv = 0.16;
  TLatex* lat2 = new TLatex(0.25,yv,"Excluded");
  lat2->SetNDC(true);
  lat2->SetTextSize(0.04);
  lat2->Draw("same");

  TGraph* noRegion = new TGraph(3);
  noRegion->SetPoint(0,TMath::Min(xmin,ymin),TMath::Min(xmin,ymin));
  noRegion->SetPoint(1,xmax,ymin);
  noRegion->SetPoint(2,TMath::Min(xmax,ymax),TMath::Min(xmax,ymax));
  noRegion->SetFillColor(16);

  TLatex* lat3 = new TLatex(0.7,0.25,"#tilde{g} NLSP");
  lat3->SetNDC(true);
  lat3->SetTextSize(0.04);

  if(bino.Contains("mNScan")){
    noRegion->Draw("same f");
    lat3->Draw("same");
  }

  if(print) {
    can_excl02->Print("",".gif");
    can_excl02->Print("",".pdf");
  }

  fout->cd();
  fout->Write();

  can_acc->Write();
  can_xs->Write();
  can_limit->Write();
  can_excl01->Write();
  can_excl02->Write();

  for(int i=0; i<nlimit; i++)
    for(int j=0; j<nxs-1; j++)
      curv[i][j]->Write();



}



