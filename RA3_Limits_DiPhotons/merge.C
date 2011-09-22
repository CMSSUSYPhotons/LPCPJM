
void merge() {

  const int N = 6;
  TString label[N] = {"bino_mN150_met100_1jet",
		      "bino_mN375_met100_1jet",
		      "bino_mN375_met100_nojet",
		      "bino_mN375_met50_1jet",
		      "wino_mN150_met100_1jet",
		      "wino_mN375_met100_1jet"};
  const int M = 2;
  TString model[M] = {"0","1"};

  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      TChain ch("gridPoints");
      ch.Add("/uscmst1b_scratch/lpc1/lpctau/dwjang/work/jobs/limits/"+label[i]+"/tree_*model"+model[j]+".root");
      ch.Merge("tree_"+label[i]+"_model"+model[j]+".root");
    }
  }

  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      TFile* f = new TFile("tree_"+label[i]+"_model"+model[j]+".root","READ");
      TTree* t = (TTree*) f->Get("gridPoints");
      std::cout << label[i] << ", " << model[j] << " : " << t->GetEntries() << std::endl;
    }
  }

}
