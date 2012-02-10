
void comp() {

  gStyle->SetOptStat(0);

  const int N = 6;
  TString label[N] = {"bin0","bin1","bin2","bin3","bin4","multiChannel"};
  TString bin_label[N] = {"50-60","60-70","70-80","80-100","100-","combined"};
  int color[N] = {2,3,4,5,6,1};

  TGraph* graph[N];

  for(int i=0; i<N; i++) {
    TFile* f = new TFile("hist_exclusion_bino_mN375_met100_1jet_"+label[i]+".root","READ");
    graph[i] = (TGraph*) f->Get("excl_curv_limit_xsec");
    //graph[i] = (TGraph*) f->Get("excl_curv_exp_xsec");
    graph[i]->SetLineColor(color[i]);
  }

  TH2D* h_back = new TH2D("h_back",";M_{#tilde{q}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,400,2000,100,400,2000);

  TLegend* leg = new TLegend(0.7,0.6,0.9,0.9,"Bin Ranges");
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  TCanvas* can = new TCanvas("can_binComp","can_binComp",800,600);
  h_back->GetYaxis()->SetTitleOffset(1.2);
  h_back->Draw();
  for(int i=0; i<N; i++){
    graph[i]->Draw("L same");
    leg->AddEntry(graph[i],bin_label[i],"L");
  }
  leg->Draw("same");

  can->Print("",".gif");
  can->Print("",".pdf");

}

