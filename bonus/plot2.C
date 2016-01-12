void plot2(){

  TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length for various momenta",50,0,50,1000,0.5,180.0);
  hframe->Draw(); //you can set the axis att via hframe->GetXaxis()..
  gPad->SetLogx(0);
  hframe->GetXaxis()->SetTitle("Track Length [mm]");
  hframe->GetYaxis()->SetTitle("Total Energy [MeV]");
  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(20);
  
  TTree *MyTree = new TTree("MyTree", "MyTree");
  MyTree->ReadFile("Rout_proton65mevOnSolenoid3LIV.txt","track:edep");
  MyTree->SetMarkerStyle(kCircle);
  MyTree->SetMarkerColorAlpha(kRed, 0.3);
  MyTree->Draw("edep:track","","same");
  TTree *MyTreeB = new TTree("MyTreeB", "MyTreeB");
  MyTreeB->ReadFile("Rout_proton65mevOnSolenoid3STD.txt","trackB:edepB");
  MyTreeB->SetMarkerStyle(kFullCircle);
  MyTreeB->SetMarkerColorAlpha(kRed, 0.3);
  MyTreeB->Draw("edepB:trackB","","same");
  Tl.DrawLatex(5, 25,   "p_{s} = 65 MeV/c @ vtx=0.0");
 
  
  TTree *MyTree1 = new TTree("MyTree1", "MyTree1");
  MyTree1->ReadFile("Rout_proton72mevOnSolenoid3LIV.txt","track1:edep1");
  MyTree1->SetMarkerStyle(kOpenSquare);
  MyTree1->SetMarkerColorAlpha(kBlue, 0.3);
  MyTree1->Draw("edep1:track1","","same");
  TTree *MyTree1B = new TTree("MyTree1B", "MyTree1B");
  MyTree1B->ReadFile("Rout_proton72mevOnSolenoid3STD.txt","track1B:edep1B");
  MyTree1B->SetMarkerStyle(kFullSquare);
  MyTree1B->SetMarkerColorAlpha(kBlue, 0.3);
  MyTree1B->Draw("edep1B:track1B","","same");
  Tl.DrawLatex(10, 43,   "p_{s} = 72 MeV/c");

  TTree *MyTree2 = new TTree("MyTree2", "MyTree2");
  MyTree2->ReadFile("Rout_proton90mevOnSolenoid3LIV.txt","track2:edep2");
  MyTree2->SetMarkerStyle(kOpenTriangleUp);
  MyTree2->SetMarkerColorAlpha(kGreen, 0.3);
  MyTree2->Draw("edep2:track2","","same");
  TTree *MyTree2B = new TTree("MyTree2B", "MyTree2B");
  MyTree2B->ReadFile("Rout_proton90mevOnSolenoid3STD.txt","track2B:edep2B");
  MyTree2B->SetMarkerStyle(kFullTriangleUp);
  MyTree2B->SetMarkerColorAlpha(kGreen, 0.3);
  MyTree2B->Draw("edep2B:track2B","","same");
  Tl.DrawLatex(15, 70,   "p_{s} = 90 MeV/c");

  TTree *MyTree3 = new TTree("MyTree3", "MyTree3");
  MyTree3->ReadFile("Rout_proton120mevOnSolenoid3LIV.txt","track3:edep3");
  MyTree3->SetMarkerStyle(kOpenDiamond);
  MyTree3->SetMarkerColorAlpha(kMagenta, 0.3);
  MyTree3->Draw("edep3:track3","","same");
  TTree *MyTree3B = new TTree("MyTree3B", "MyTree3B");
  MyTree3B->ReadFile("Rout_proton120mevOnSolenoid3STD.txt","track3B:edep3B");
  MyTree3B->SetMarkerStyle(kFullDiamond);
  MyTree3B->SetMarkerColorAlpha(kMagenta, 0.3);
  MyTree3B->Draw("edep3B:track3B","","same");
  Tl.DrawLatex(20, 105,   "p_{s} = 120 MeV/c");

  TTree *MyTree4 = new TTree("MyTree4", "MyTree4");
  MyTree4->ReadFile("Rout_proton180mevOnSolenoid3LIV.txt","track4:edep4");
  MyTree4->SetMarkerStyle(kOpenTriangleDown);
  MyTree4->SetMarkerColorAlpha(kBlack, 0.3);
  MyTree4->Draw("edep4:track4","","same");
  TTree *MyTree4B = new TTree("MyTree4B", "MyTree4B");
  MyTree4B->ReadFile("Rout_proton180mevOnSolenoid3STD.txt","track4B:edep4B");
  MyTree4B->SetMarkerStyle(kFullTriangleDown);
  MyTree4B->SetMarkerColorAlpha(kBlack, 0.3);
  MyTree4B->Draw("edep4B:track4B","","same");
  Tl.DrawLatex(25, 170,   "p_{s} = 180 MeV/c");

  Tl.DrawLatex(10, 142,   "open:LIV, full:STD,(0.1% lower)");

}
