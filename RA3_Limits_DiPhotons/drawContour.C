TGraph* getContour(TH2D* h, TString name) {

  TGraph* graph = new TGraph(100);
  graph->SetName(name);
  int ip = 0;
  int nx = h->GetXaxis()->GetNbins();
  int ny = h->GetYaxis()->GetNbins();

  // for y>x
  int ix = -1;
  for(int j=ny;true; j--) {
    int k = -1;
    for(int i=2; i<nx-1; i++) {
      if(h->GetBinContent(i,j) < 0) {
	std::cout << "i,j,z : " << i << ", " << j << ", " << h->GetBinContent(i,j) << std::endl;
	k = i;
	break;
      }
    }// for i
    if(k<0) continue;
    double y = h->GetYaxis()->GetBinLowEdge(j);
    double x1 = h->GetXaxis()->GetBinLowEdge(k-1);
    double x2 = h->GetXaxis()->GetBinLowEdge(k);
    double z1 = h->GetBinContent(k-1,j);
    double z2 = h->GetBinContent(k,j);
    double x = x1 + (x2-x1)*fabs(z1)/fabs(z2-z1);
    std::cout << "y, x1, x2, z1, z2, x : " << y << ", " << x1 << ", " << x2 << ", " << x << ", " << z1 << ", " << z2 << std::endl;
    graph->SetPoint(ip++,x,y);

    if(h->GetYaxis()->GetBinLowEdge(j) < 1100) {
      ix = k;
      break;
    }
  }// for j

  if(ix < 0) std::cout << "Something wrong...." << std::endl;

  // for y<x
  for(int i=ix; i<=nx; i++) {
    int k = -1;
    for(int j=2; j<ny-1; j++) {
      if(h->GetBinContent(i,j) < 0) {
	k = j;
	break;
      }
    }// for j
    if(k<0) continue;
    double x = h->GetXaxis()->GetBinLowEdge(i);
    double y1 = h->GetYaxis()->GetBinLowEdge(k-1);
    double y2 = h->GetYaxis()->GetBinLowEdge(k);
    double z1 = h->GetBinContent(i,k-1);
    double z2 = h->GetBinContent(i,k);
    double y = y1 + (y2-y1)*fabs(z1)/fabs(z2-z1);
    std::cout << "x, y1, y2, z1, z2, y : " << x << ", " << y1 << ", " << y2 << ", " << y << ", " << z1 << ", " << z2 << std::endl;
    graph->SetPoint(ip++,x,y);
  }// for i

  ip = graph->GetN()-1;
  while(1) {
    double x, y;
    graph->GetPoint(ip,x,y);
    if(x>1) break;
    else graph->RemovePoint(ip);
    ip--;
  }
  

  return graph;
}

void fillPotHoles(TH2D* h) {
  int nbinsX = h->GetXaxis()->GetNbins();
  int nbinsY = h->GetYaxis()->GetNbins();

  double epsilon = 0;

  for(int ix=1; ix <= nbinsX; ix++) {
    for(int iy=1; iy <= nbinsY; iy++) {
      double val = h->GetBinContent(ix,iy);
      if(val < epsilon) {
	int ncnt = 0;
	double up    = h->GetBinContent(ix,iy+1);
	if(up > epsilon) ncnt++;
	double down  = h->GetBinContent(ix,iy-1);
	if(down > epsilon) ncnt++;
	double left  = h->GetBinContent(ix-1,iy);
	if(left > epsilon) ncnt++;
	double right = h->GetBinContent(ix+1,iy);
	if(right > epsilon) ncnt++;
	if(ncnt == 4){
	  val = (up+down+left+right)/ncnt;
	  h->SetBinContent(ix,iy,val);
	  up    = h->GetBinError(ix,iy+1);
	  down  = h->GetBinError(ix,iy-1);
	  left  = h->GetBinError(ix-1,iy);
	  right = h->GetBinError(ix+1,iy);
	  val = std::sqrt(up*up + down*down + left*left + right*right)/ncnt;
	  h->SetBinError(ix,iy,val);
	}
      }
      else {
	int ncnt = 0;
	double up    = h->GetBinContent(ix,iy+1);
	if(up < epsilon) ncnt++;
	double down  = h->GetBinContent(ix,iy-1);
	if(down < epsilon) ncnt++;
	double left  = h->GetBinContent(ix-1,iy);
	if(left < epsilon) ncnt++;
	double right = h->GetBinContent(ix+1,iy);
	if(right < epsilon) ncnt++;
	if(ncnt == 4){
	  val = (up+down+left+right)/ncnt;
	  h->SetBinContent(ix,iy,val);
	  up    = h->GetBinError(ix,iy+1);
	  down  = h->GetBinError(ix,iy-1);
	  left  = h->GetBinError(ix-1,iy);
	  right = h->GetBinError(ix+1,iy);
	  val = std::sqrt(up*up + down*down + left*left + right*right)/ncnt;
	  h->SetBinError(ix,iy,val);
	}
      }
    } // for iy
  } // for ix
}


void RemovePoints(TGraph* g) {
  int N = g->GetN();
  int last = -1;
  for(int i=0; i<N; i++) {
    double x,y;
    g->GetPoint(i,x,y);
    //    if(abs(x-y) < 0.1) last = i;
    if(x > 810) {
      last = i+1;
      break;
    }
  }
  for(int i=N-1; i>=last; i--) g->RemovePoint(i);
}


void drawContour(TString bino="bino", TString mNu="375", TString metCut="met100", TString jetCut="1jet", TString model="1") {

  TString label = bino + "_mN" + mNu + "_" + metCut + "_" + jetCut + "_model" + model;

  float xymin = 400;
  if(bino.Contains("wino") && mNu.Contains("150")) xymin = 150;

  gStyle->SetPalette(1);
  gStyle->SetPadLeftMargin(0.15);

  TFile* fin = new TFile("exclusion_results_"+label+".root","READ");
  TH2D* h_acc = (TH2D*) fin->Get("acc");

  const int nxs = 5;
  TString xsname[nxs] = {"xsec","xsec_1L","xsec_1H","xsec_2L","xsec_2H"};
  TH2D* h_xs[nxs];
  for(int i=0; i<nxs; i++) {
    h_xs[i] = (TH2D*) fin->Get(xsname[i]);
  }

  const int nlimit = 6;
  TString limitname[nlimit] = {"limit", "exp", "exp_1L","exp_1H","exp_2L","exp_2H"};
  TH2D* h_limit[nlimit];
  for(int i=0; i<nlimit; i++) {
    h_limit[i] = (TH2D*) fin->Get(limitname[i]);
  }


  TFile* fout = new TFile("exclusion_plots_"+label+".root","RECREATE");

  TCanvas* can_acc = new TCanvas("can_acc_"+label,"can_acc_"+label,1000,800);
  can_acc->SetRightMargin(0.15);
  h_acc->GetXaxis()->SetNdivisions(505);
  h_acc->GetYaxis()->SetNdivisions(505);
  h_acc->Draw("COL Z");
  can_acc->Print("",".gif");
  can_acc->Print("",".eps");

  TCanvas* can_xs = new TCanvas("can_xsec_"+label,"can_xsec_"+label,1000,800);
  can_xs->SetRightMargin(0.15);
  h_xs[0]->GetXaxis()->SetNdivisions(505);
  h_xs[0]->GetYaxis()->SetNdivisions(505);
  h_xs[0]->Draw("COL Z");
  can_xs->Print("",".gif");
  can_xs->Print("",".eps");

  TCanvas* can_limit = new TCanvas("can_limit_"+label,"can_limit_"+label,1000,800);
  can_limit->SetRightMargin(0.15);
  h_limit[0]->GetXaxis()->SetNdivisions(505);
  h_limit[0]->GetYaxis()->SetNdivisions(505);
  h_limit[0]->Draw("COL Z");
  can_limit->Print("",".gif");
  can_limit->Print("",".eps");


  // now find exclusion curves
  TCanvas* can_diff = new TCanvas("can_diff_"+label,"can_diff_"+label,1200,800);
  //  can_diff->Divide(nlimit,nxs);
  TH2D* h_excl[nlimit][nxs];
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs; j++) {

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
// 	  if(i==0 && j==0 && x2-x1 < 0) {
// 	    std::cout << ibx << "," << iby << "=("
// 		      << h_excl[i][j]->GetXaxis()->GetBinLowEdge(ibx) << ","
// 		      << h_excl[i][j]->GetYaxis()->GetBinLowEdge(iby) << ") : "
// 		      << "(" << x2 << "-" << x1 << ") = " << x2-x1 << std::endl;
// 	  }// if
	}// for iby
      }// for ibx
      fillPotHoles(h_excl[i][j]);
      //      can_diff->cd(j*nlimit + i+1);
      if(i==0 && j==0) h_excl[i][j]->Draw("TEXT");
    }// for j
  }// for i

  TGraph* curv[nlimit][nxs];

  TH2D* h_back = new TH2D("h_back",";M_{#tilde{q}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,xymin,2000,100,xymin,2000);
  h_back->GetXaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetTitleOffset(1.2);

  double contours[2]={ 0.0, 1.1 };
  TCanvas *can_excl01 = new TCanvas("can_excl01_"+label, "can_excl01_"+label,1200,800);
  can_excl01->Divide(nlimit,nxs);
  for(int i=0; i<nlimit; i++) {
    for(int j=0; j<nxs; j++) {
      can_excl01->cd(j*nlimit + i + 1);
      h_back->Draw();
      curv[i][j] = getContour(h_excl[i][j],"excl_curv_"+limitname[i]+"_"+xsname[j]);
      curv[i][j]->Draw("SAME L");
//       h_excl[i][j]->SetContour(2,contours);
//       h_excl[i][j]->Draw("SAME CONT LIST");
//       gPad->Update();

//       TObjArray *contsM = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
//       TList* contLevel = (TList*)contsM->At(0);
//       curv[i][j] = (TGraph*)contLevel->First()->Clone("excl_curv_"+limitname[i]+"_"+xsname[j]);
//       //      RemovePoints(curv[i][j]);
    }// for j
  }// for i


  TGraph* excludedRegion = new TGraph(curv[0][0]->GetN()+3);
  int nbins = curv[0][0]->GetN();
  for(int i=0; i<nbins; i++){
    double x,y;
    curv[0][0]->GetPoint(i,x,y);
    excludedRegion->SetPoint(i,x,y);
  }
  if(bino.Contains("wino") && mNu.Contains("375")){
    excludedRegion->SetPoint(nbins,2000,xymin);
    excludedRegion->SetPoint(nbins+1,xymin,xymin);
    excludedRegion->RemovePoint(nbins+2);
  }
  else {
    excludedRegion->SetPoint(nbins,2000,xymin);
    excludedRegion->SetPoint(nbins+1,xymin,xymin);
    excludedRegion->SetPoint(nbins+2,xymin,2000);
  }

  TGraph* oneSigma = new TGraph(curv[0][1]->GetN()+curv[0][2]->GetN());
  TGraph* twoSigma = new TGraph(curv[0][3]->GetN()+curv[0][4]->GetN());

  for(int i=1; i<nxs; i++) {
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


  TLegend* leg = new TLegend(0.57,0.5,0.95,0.75,bino+"-like,  m_{#tilde{#Chi}^{0}} = "+mNu+" (GeV/c^{2})");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(curv[0][0],"Observed, NLO","L");
  leg->AddEntry(curv[1][0],"Expected, NLO","L");
  leg->AddEntry(oneSigma,"#pm1#sigma, NLO","F");
  //  leg->AddEntry(twoSigma,"#pm2#sigma, NLO","F");
  leg->Draw("same");

  TLatex* lat = new TLatex(0.57,0.85,"CMS Preliminary");
  lat->SetNDC(true);
  lat->SetTextSize(0.04);
  lat->Draw("same");

  TLatex* lat2 = new TLatex(0.57,0.77,"#sqrt{s} = 7 TeV,  #int#it{L}dt = 1.14 fb^{-1}");
  lat2->SetNDC(true);
  lat2->SetTextSize(0.04);
  lat2->Draw("same");

  float yv = 0.25;
  if(bino.Contains("wino")) yv = 0.16;
  TLatex* lat2 = new TLatex(0.25,yv,"Excluded");
  lat2->SetNDC(true);
  lat2->SetTextSize(0.04);
  lat2->Draw("same");

  can_excl02->Print("",".gif");
  can_excl02->Print("",".eps");

  fout->cd();
  fout->Write();

  can_acc->Write();
  can_xs->Write();
  can_limit->Write();
  can_excl01->Write();
  can_excl02->Write();

  for(int i=0; i<nlimit; i++)
    for(int j=0; j<nxs; j++)
      curv[i][j]->Write();

}

