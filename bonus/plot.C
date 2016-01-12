void plot(){

  //TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length p=180MeV/C",50,0,50,100,178.0,179.0);
  //TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length p=120MeV/c",50,0,50,100,115.2,116.8);
  //  TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length p=90MeV/c",500,0,50,100,79.5,83.0);
  TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length p=72MeV/c",500,0,50,100,45.5,56.1);
  //  TH2F *hframe = new TH2F("hframe","Edep_tot vs. track length p=65MeV/c",500,0,50,100,0.5,46.1);
  TTree *MyTree = new TTree("MyTree", "MyTree");
  MyTree->ReadFile("Rout_proton72mevOnSolenoid3LIV.txt","track:edep");
  MyTree->SetMarkerStyle(kCircle);
  MyTree->SetMarkerColorAlpha(kRed, 0.35);
  hframe->Draw(); //you can set the axis att via hframe->GetXaxis()..
  MyTree->Draw("edep:track","","same");


  /*
  TTree *MyTree1 = new TTree("MyTree1", "MyTree1");
  MyTree1->ReadFile("Rout_proton120mevOnSolenoid3STD.txt","track1:edep1");
  MyTree1->SetMarkerStyle(kOpenSquare);
  MyTree1->SetMarkerColorAlpha(kBlue, 0.35);
  MyTree1->Draw("edep1:track1","","same");
  */
}
