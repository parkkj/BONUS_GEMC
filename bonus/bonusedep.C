void bonusedep()
{

 TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(20);
  
  //  TH2F *hframe = new TH2F("hframe","Step length vs. track length",100,10.,80,100,0.00001,0.1);
  TH2F *hframe = new TH2F("hframe","GEMC BONUS12 simulation",100,25.,75.,100,0.,10.);
  hframe->Draw(); //you can set the axis att via hframe->GetXaxis()..
  TTree *MyTree = new TTree("MyTree", "MyTree");
  //  MyTree->ReadFile("R70mevstraighttrack2L.txt","track:edep");
  MyTree->ReadFile("R70mevstraighttrack2L.txt","edep:lstep");
  MyTree->SetMarkerStyle(kFullCircle);
  MyTree->SetMarkerColorAlpha(kRed, 0.2);
  MyTree->Draw("edep*1000:lstep","","same");
  hframe->GetXaxis()->SetTitle("Track Length [mm]");
  hframe->GetYaxis()->SetTitle("Edep [keV]");
  gPad->SetLogy(1);
  
  TTree *MyTree2 = new TTree("MyTree2", "MyTree2");
  MyTree2->ReadFile("R70mevstraighttrack3L.txt","edep2:lstep2");
  MyTree2->SetMarkerStyle(kOpenCircle);
  MyTree2->SetMarkerColorAlpha(kBlue, 0.2);
  MyTree2->Draw("edep2*1000:lstep2","","same");

  Tl.DrawLatex(10, 17500,   "red open circle:TimeWindow=0ns");

  }
