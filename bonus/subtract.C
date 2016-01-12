{
   TFile *f1 = TFile::Open("bonus_default_gemc_out.root");
   TCanvas *c1 = (TCanvas*)f1->Get("c1");
   TH1D *h1 = (TH1D*)c1->FindObject("cutg");
   TFile *f2 = TFile::Open("bonus_rohacell_gemc_out.root");
   TCanvas *c2 = (TCanvas*)f1->Get("c1");
   TH1D *h2 = (TH1D*)c2->FindObject("cutg");
   h1->Add(h2,-1);
   
}
