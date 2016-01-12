void comp(){

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(20);
  
  c1 = new TCanvas("c1","Helix vs. KF fits",600,400);
  
  TFile *f2 = TFile::Open("kfilt.root");
  TH1F *h2 = new TH1F("hframe2","KF fit",500,-10.,5.);
  ntuple2->Draw("pTresol_kf>>hframe2");
  h2->Draw();
  h2->SetLineColorAlpha(kRed, 0.6);
  h2->GetXaxis()->SetTitle("p_{T} resol.");

  TFile *f1 = TFile::Open("helix.root");  
  TH1F *h1 = new TH1F("hframe1","Helix fit",500,-10.,5.);
  ntuple1->Draw("pTresol_helix>>hframe1","","same");
  h1->Draw("same");
  h1->SetLineColorAlpha(kBlue, 0.6);
  gPad->SetLogy(1);

  
  Tl.DrawLatex(0, 30,   "red: p_{T} KF fit");
  Tl.DrawLatex(0, 20,   "blue: p_{T} Helix fit");



  c2 = new TCanvas("c2","Helix vs. KF fits",600,400);
  
  TFile *f2 = TFile::Open("kfilt.root");
  TH1F *h2 = new TH1F("hframe2","KF fit",500,-15.,15.);
  ntuple2->Draw("thetaDiff_kf>>hframe2");
  h2->Draw();
  h2->SetLineColorAlpha(kRed, 0.2);
  h2->GetXaxis()->SetTitle("#Delta #theta (mrad)");

  TFile *f1 = TFile::Open("helix.root");  
  TH1F *h1 = new TH1F("hframe1","Helix fit",500,-15.,15.);
  ntuple1->Draw("thetaDiff_helix>>hframe1","","same");
  h1->Draw("same");
  h1->SetLineColorAlpha(kBlue, 0.2);
  gPad->SetLogy(1);

  
  Tl.DrawLatex(0, 30,   "red: #theta KF fit");
  Tl.DrawLatex(0, 20,   "blue: #theta Helix fit");









  
  c3 = new TCanvas("c3","Helix vs. KF fits",600,400);
  
  TFile *f2 = TFile::Open("kfilt.root");
  TH1F *h2 = new TH1F("hframe2","KF fit",500,-150.,150.);
  ntuple2->Draw("phiDiff_kf>>hframe2");
  h2->Draw();
  h2->SetLineColorAlpha(kRed, 0.6);
  h2->GetXaxis()->SetTitle("#Delta #phi (mrad)");

  TFile *f1 = TFile::Open("helix.root");  
  TH1F *h1 = new TH1F("hframe1","Helix fit",500,-150.,150.);
  ntuple1->Draw("phiDiff_helix>>hframe1","","same");
  h1->Draw("same");
  h1->SetLineColorAlpha(kBlue, 0.6);
  gPad->SetLogy(1);

  
  Tl.DrawLatex(0, 30,   "red: #phi KF fit");
  Tl.DrawLatex(0, 20,   "blue: #phi Helix fit");



  
}
