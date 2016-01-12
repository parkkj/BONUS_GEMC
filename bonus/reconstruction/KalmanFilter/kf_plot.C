void kf_plot(){

  TH1F *hframe = new TH1F("hframe","momentum resolution",100,-10.,10.);
  hframe->Draw(); //you can set the axis att via hframe->GetXaxis()..
  gPad->SetLogx(0);
  hframe->GetXaxis()->SetTitle(" #Delta p_T /p_T [%]");
  hframe->GetYaxis()->SetTitle("Number of event");
  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(20);

  double pTresol,pXresol,pYresol;
  int nline =0;
  TTree *MyTree = new TTree("MyTree", "MyTree");
  MyTree->ReadFile("helix_out.txt","pTresol:pXresol:pYresol");
  MyTree->SetMarkerStyle(kCircle);
  MyTree->SetMarkerColorAlpha(kRed, 0.3);
  while(1){
    hframe->Fill(pTresol);
    nline++;
  }
  hframe->Draw("same");

  Tl.DrawLatex(10, 142,   "open: p_T resolution(helix)");

}
