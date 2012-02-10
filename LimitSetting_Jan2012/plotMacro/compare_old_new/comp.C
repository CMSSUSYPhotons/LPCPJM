

void comp() {

  TFile* f_new = new TFile("hist_exclusion_bino_mN375_met100_1jet_new.root","READ");
  TFile* f_old = new TFile("hist_exclusion_bino_mN375_met100_1jet_old.root","READ");

  TGraph* h_new_obs = (TGraph*) f_new->Get("excl_curv_limit_xsec");
  TGraph* h_new_exp = (TGraph*) f_new->Get("excl_curv_exp_xsec");

  h_new_obs->SetLineColor(2);
  h_new_exp->SetLineColor(2);

  TGraph* h_old_obs = (TGraph*) f_old->Get("excl_curv_limit_xsec");
  TGraph* h_old_exp = (TGraph*) f_old->Get("excl_curv_exp_xsec");

  TH2D* h_back = new TH2D("h_back",";M_{#tilde{q}} (GeV/c^{2});M_{#tilde{g}} (GeV/c^{2})",100,400,2000,100,400,2000);

  TLegend* leg_obs = new TLegend(0.2,0.3,0.5,0.5,"Observed Limit");
  leg_obs->SetLineColor(0);
  leg_obs->SetFillColor(0);

  TCanvas* can_obs = new TCanvas("can_obs","can_obs",800,600);
  h_back->Draw();
  h_new_obs->Draw("L same");
  h_old_obs->Draw("L same");
  leg_obs->AddEntry(h_new_obs,"NewTool","L");
  leg_obs->AddEntry(h_old_obs,"OldTool","L");
  leg_obs->Draw("same");
  can_obs->Print("",".gif");
  can_obs->Print("",".pdf");

  TLegend* leg_exp = new TLegend(0.2,0.3,0.5,0.5,"Expected Limit");
  leg_exp->SetLineColor(0);
  leg_exp->SetFillColor(0);

  TCanvas* can_exp = new TCanvas("can_exp","can_exp",800,600);
  h_back->Draw();
  h_new_exp->Draw("L same");
  h_old_exp->Draw("L same");
  leg_exp->AddEntry(h_new_exp,"NewTool","L");
  leg_exp->AddEntry(h_old_exp,"OldTool","L");
  leg_exp->Draw("same");
  can_exp->Print("",".gif");
  can_exp->Print("",".pdf");

}
