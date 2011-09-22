
void comp_limit(TString model="gauss") {

  gStyle->SetPadLeftMargin(0.15);

  const int N = 4;
  TString label[N] = {"1jet","nojet","ff","ee"};
  int lcolor[N] = {1,3,4,2};
  int lstyle[N] = {1,2,5,9};
  TFile* f[N];
  TGraph* g[N];
  for(int i=0; i<N; i++) {
    f[i] = new TFile("exclusion_plots_"+model+"_met100_"+label[i]+".root","READ");
    g[i] = (TGraph*) f[i]->Get("excl_curv_limit_xsec");
    g[i]->SetLineColor(lcolor[i]);
    g[i]->SetLineStyle(lstyle[i]);
    g[i]->SetLineWidth(4);
  }

  TH2D* h_back = new TH2D("h_back",";M_{#tilde{q}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,400,2000,100,400,2000);
  h_back->GetXaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetNdivisions(505);
  h_back->GetYaxis()->SetTitleOffset(1.2);

  TLegend* leg = new TLegend(0.65,0.55,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetLineColor(0);

  TCanvas* can = new TCanvas("can_comp_limit","can_comp_limit",1000,800);
  h_back->Draw();
  for(int i=0; i<N; i++) {
    g[i]->Draw("SAME L");
    leg->AddEntry(g[i],label[i],"L");
  }
  leg->Draw("SAME");
  can->Print("",".gif");
  can->Print("",".eps");

}
