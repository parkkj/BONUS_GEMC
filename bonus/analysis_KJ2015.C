#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void analysis_KJ2015(string input_filename)
{
    // gROOT->Reset();
    gStyle->SetPalette(1);
    gStyle->SetOptStat(1);
  
    const double DEG=180./3.1415926;

    // bool Is_PVDIS=false,Is_SIDIS_3he=false,Is_SIDIS_proton=false,Is_JPsi=false;
    // if (input_filename.find("PVDIS",0) != string::npos) Is_PVDIS=true;
    // else if (input_filename.find("SIDIS_3he",0) != string::npos) Is_SIDIS_3he=true;
    // else if (input_filename.find("SIDIS_proton",0) != string::npos) Is_SIDIS_proton=true;
    // else if (input_filename.find("JPsi",0) != string::npos) Is_JPsi=true;
    // else {cout << "not PVDIS or SIDIS or JPsi " << endl; return;}

    char the_filename[200];
    sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());

    char output_filename[200];
    sprintf(output_filename, "%s_output.root",the_filename);
    TFile *outputfile=new TFile(output_filename, "recreate");

    TH2F *hacceptance_forwardangle,*hacceptance_largeangle,*hacceptance_overall;


    // # == id name ==============================================================
    // # digit     beamline                side          magnet      number          window
    // #        ion    1        upstream    1    dipole    1           n     front      1
    // #        ele    2       downstream   2  quadrupole  2                 back       2
    // # =========================================================================

    char *name_beamline[2]={"ion","ele"};
    char *name_side[2]={"upstream","downstream"};
    char *name_magnet[2]={"dipole","quadrupole"};
    char *name_number[8]={"1","2","3","4","5","6","7","8"};
    char *name_window[2]={"front","back"};

    TH2F *hhit_alphaXY[2][2][2][8][2];
    TH2F *hhit_ThetaP[2][2][2][8][2],*hhit_PhiP[2][2][2][8][2],*hhit_ThetaPhi[2][2][2][8][2],*hhit_xyl[2][2][2][8][2];
    for (int h=0;h<2;h++){    
	for (int i=0;i<2;i++){
	    for (int j=0;j<2;j++){
		for (int k=0;k<8;k++){
		    for (int l=0;l<2;l++){ 
			hhit_alphaXY[h][i][j][k][l]=new TH2F(Form("hit_alphaXY_%i_%i_%i_%i_%i",h,i,j,k,l),Form("hit of %s %s %s %s %s; #alpha_{x} (mrad);#alpha_{y} (mrad)",name_beamline[h],name_side[i],name_magnet[j],name_number[k],name_window[l]),1000,-120,-80,1000,-30,30);
			hhit_ThetaPhi[h][i][j][k][l]=new TH2F(Form("hit_ThetaPhi_%i_%i_%i_%i_%i",h,i,j,k,l),Form("hit of %s %s %s %s %s; vertex #phi (deg);vertex #theta (deg)",name_beamline[h],name_side[i],name_magnet[j],name_number[k],name_window[l]),360,-180,180,180,0,20);

		    }}}}}

    TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","gen_ThetaP;vertex P (GeV);vertex #theta (deg)",1000,0,100,1000,0,180);
    TH2F *hgen_ThetaPhi=new TH2F("gen_ThetaPhi","gen_ThetaPhi;vertex #phi (deg);vertex #theta (deg)",360,-180,180,180,0,180);
    TH2F *hgen_PhiP=new TH2F("gen_PhiP","gen_PhiP;vertex P (GeV);vertex #phi (deg)",1000,0,100,1000,-180,180);

    TH2F *hgen_ThetaP_ion=new TH2F("gen_ThetaP_ion","gen_ThetaP_ion;vertex P (GeV);vertex #theta (mrad)",1000,0,100,1000,0,100);
    TH2F *hgen_ThetaPhi_ion=new TH2F("gen_ThetaPhi_ion","gen_ThetaPhi_ion;vertex #phi (mrad);vertex #theta (mrad)",1000,-3142,3142,1000,0,100);
    TH2F *hgen_PhiP_ion=new TH2F("gen_PhiP_ion","gen_PhiP_ion;vertex P (GeV);vertex #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hgen_Thetat_ion=new TH2F("gen_Thetat_ion","gen_Thetat_ion;-t (GeV^{2});vertex #theta (mrad)",100,0,5,100,0,100);

    //    Ion QUADRUPOL 1
    TH2F *hacc_ThetaP_ionquad1back=new TH2F("acc_ThetaP_ionquad1back","acc_ThetaP_ionquad1back;vertex P (GeV);vertex #theta (mrad)",1000,0,100,1000,0,100);
    TH2F *hacc_ThetaPhi_ionquad1back=new TH2F("acc_ThetaPhi_ionquad1back","acc_ThetaPhi_ionquad1back;vertex #phi (mrad);vertex #theta (mrad)",1000,-3142,3142,1000,0,100);
    TH2F *hacc_PhiP_ionquad1back=new TH2F("acc_PhiP_ionquad1back","acc_PhiP_ionquad1back;vertex P (GeV);vertex #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hacc_Thetat_ionquad1back=new TH2F("acc_Thetat_ionquad1back","acc_Thetat_ionquad1back;-t (GeV^{2});vertex #theta (mrad)",1000,0,5,1000,0,1000);

    TH2F *hhit_ThetaP_ionquad1back=new TH2F("hit_ThetaP_ionquad1back","hit_ThetaP_ionquad1back; P (GeV); #theta (mrad)",1000,0,100,1000,-157,157);
    TH2F *hhit_ThetaPhi_ionquad1back=new TH2F("hit_ThetaPhi_ionquad1back","hit_ThetaPhi_ionquad1back; #phi (mrad); #theta (mrad)",1000,-3142,3142,1000,40,60);
    TH2F *hhit_PhiP_ionquad1back=new TH2F("hit_PhiP_ionquad1back","hit_PhiP_ionquad1back; P (GeV); #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hhit_Thetat_ionquad1back=new TH2F("hit_Thetat_ionquad1back","hit_Thetat_ionquad1back;-t (GeV^{2}); #theta (mrad)",100,0,5,100,0,100);


    //    Ion QUADRUPOL 2
    TH2F *hacc_ThetaP_ionquad2back=new TH2F("acc_ThetaP_ionquad2back","acc_ThetaP_ionquad2back;vertex P (GeV);vertex #theta (mrad)",1000,0,100,1000,0,100);
    TH2F *hacc_ThetaPhi_ionquad2back=new TH2F("acc_ThetaPhi_ionquad2back","acc_ThetaPhi_ionquad2back;vertex #phi (mrad);vertex #theta (mrad)",1000,-3142,3142,1000,0,100);
    TH2F *hacc_PhiP_ionquad2back=new TH2F("acc_PhiP_ionquad2back","acc_PhiP_ionquad2back;vertex P (GeV);vertex #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hacc_Thetat_ionquad2back=new TH2F("acc_Thetat_ionquad2back","acc_Thetat_ionquad2back;-t (GeV^{2});vertex #theta (mrad)",1000,0,5,1000,0,1000);

    TH2F *hhit_ThetaP_ionquad2back=new TH2F("hit_ThetaP_ionquad2back","hit_ThetaP_ionquad2back; P (GeV); #theta (mrad)",1000,0,100,1000,-157,157);
    TH2F *hhit_ThetaPhi_ionquad2back=new TH2F("hit_ThetaPhi_ionquad2back","hit_ThetaPhi_ionquad2back; #phi (mrad); #theta (mrad)",1000,-3142,3142,1000,40,60);
    TH2F *hhit_PhiP_ionquad2back=new TH2F("hit_PhiP_ionquad2back","hit_PhiP_ionquad2back; P (GeV); #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hhit_Thetat_ionquad2back=new TH2F("hit_Thetat_ionquad2back","hit_Thetat_ionquad2back;-t (GeV^{2}); #theta (mrad)",100,0,5,100,0,100);


    //    Ion QUADRUPOL 3
    TH2F *hacc_ThetaP_ionquad3back=new TH2F("acc_ThetaP_ionquad3back","acc_ThetaP_ionquad3back;vertex P (GeV);vertex #theta (mrad)",1000,0,100,1000,0,100);
    TH2F *hacc_ThetaPhi_ionquad3back=new TH2F("acc_ThetaPhi_ionquad3back","acc_ThetaPhi_ionquad3back;vertex #phi (mrad);vertex #theta (mrad)",1000,-3142,3142,1000,0,100);
    TH2F *hacc_PhiP_ionquad3back=new TH2F("acc_PhiP_ionquad3back","acc_PhiP_ionquad3back;vertex P (GeV);vertex #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hacc_Thetat_ionquad3back=new TH2F("acc_Thetat_ionquad3back","acc_Thetat_ionquad3back;-t (GeV^{2});vertex #theta (mrad)",1000,0,5,1000,0,1000);

    TH2F *hhit_ThetaP_ionquad3back=new TH2F("hit_ThetaP_ionquad3back","hit_ThetaP_ionquad3back; P (GeV); #theta (mrad)",1000,0,100,1000,-157,157);
    TH2F *hhit_ThetaPhi_ionquad3back=new TH2F("hit_ThetaPhi_ionquad3back","hit_ThetaPhi_ionquad3back; #phi (mrad); #theta (mrad)",1000,-3142,3142,1000,40,60);
    TH2F *hhit_PhiP_ionquad3back=new TH2F("hit_PhiP_ionquad3back","hit_PhiP_ionquad3back; P (GeV); #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hhit_Thetat_ionquad3back=new TH2F("hit_Thetat_ionquad3back","hit_Thetat_ionquad3back;-t (GeV^{2}); #theta (mrad)",100,0,5,100,0,100);





    TH2F *hacc_ThetaP_iondipole2back=new TH2F("acc_ThetaP_iondipole2back","acc_ThetaP_iondipole2back;vertex P (GeV);vertex #theta (mrad)",1000,0,100,1000,0,100);
    TH2F *hacc_ThetaPhi_iondipole2back=new TH2F("acc_ThetaPhi_iondipole2back","acc_ThetaPhi_iondipole2back;vertex #phi (mrad);vertex #theta (mrad)",1000,-3142,3142,1000,0,100);
    TH2F *hacc_PhiP_iondipole2back=new TH2F("acc_PhiP_iondipole2back","acc_PhiP_iondipole2back;vertex P (GeV);vertex #phi (mrad)",1000,0,100,1000,-3142,3142);
    TH2F *hacc_Thetat_iondipole2back=new TH2F("acc_Thetat_iondipole2back","acc_Thetat_iondipole2back;-t (GeV^{2});vertex #theta (mrad)",1000,0,5,1000,0,100);


    TH3F *hgen_ThetaPhiP_ion=new TH3F("gen_ThetaPhiP_ion","gen_ThetaPhiP_ion;vertex P (GeV);vertex #phi (mrad);vertex #theta (mrad)",100,0,100,1000,-3142,3142,100,0,100);

    TH3F *hacc_ThetaPhiP_iondipole2back=new TH3F("acc_ThetaPhiP_iondipole2back","acc_ThetaPhiP_iondipole2back;vertex P (GeV);vertex #phi (mrad);vertex #theta (mrad)",100,0,100,1000,-3142,3142,100,0,100);


    TH2F *hhit_ThetaP_iondipole2back=new TH2F("hit_ThetaP_iondipole2back","hit_ThetaP_iondipole2back; P (GeV); #theta (mrad)",100,0,100,1000,0,314);
    TH2F *hhit_ThetaPhi_iondipole2back=new TH2F("hit_ThetaPhi_iondipole2back","hit_ThetaPhi_iondipole2back; #phi (mrad); #theta (mrad)",100,-314,314,1000,0,314);
    TH2F *hhit_PhiP_iondipole2back=new TH2F("hit_PhiP_iondipole2back","hit_PhiP_iondipole2back; P (GeV); #phi (mrad)",100,0,100,1000,-314,314);
    TH2F *hhit_Thetat_iondipole2back=new TH2F("hit_Thetat_iondipole2back","hit_Thetat_iondipole2back;-t (GeV^{2}); #theta (mrad)",100,0,5,100,0,100);



    TH2F *hgen_ThetaXP_ion=new TH2F("gen_ThetaXP_ion","gen_ThetaXP_ion;P (GeV); #theta_x (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hgen_ThetaYP_ion=new TH2F("gen_ThetaYP_ion","gen_ThetaYP_ion;P (GeV); #theta_y (deg)",1000,0,100,1000,-180,180);
    TH2F *hgen_ThetaXY_ion=new TH2F("gen_ThetaXY_ion","gen_ThetaXY_ion;#theta_x (deg); #theta_y (deg)",1000,-180,180,1000,-180,180);

    TH2F *hacc_ThetaXP_iondipole2back=new TH2F("acc_ThetaXP_iondipole2back","acc_ThetaXP_iondipole2back;P (GeV); #theta_x (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaYP_iondipole2back=new TH2F("acc_ThetaYP_iondipole2back","acc_ThetaYP_iondipole2back;P (GeV); #theta_y (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaXY_iondipole2back=new TH2F("acc_ThetaXY_iondipole2back","acc_ThetaXY_iondipole2back;#theta_x (deg); #theta_y (deg)",1000,-180.,180.,1000,-180.,180.);


    //    Ion QUADRUPOL 1
    TH2F *hacc_ThetaXP_ionquad1back=new TH2F("acc_ThetaXP_ionquad1back","acc_ThetaXP_ionquad1back;P (GeV); #theta_x (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaYP_ionquad1back=new TH2F("acc_ThetaYP_ionquad1back","acc_ThetaYP_ionquad1back;P (GeV); #theta_y (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaXY_ionquad1back=new TH2F("acc_ThetaXY_ionquad1back","acc_ThetaXY_ionquad1back;#theta_x (deg); #theta_y (deg)",1000,-180.,180.,1000,-180.,180.);

    //    Ion QUADRUPOL 2
    TH2F *hacc_ThetaXP_ionquad2back=new TH2F("acc_ThetaXP_ionquad2back","acc_ThetaXP_ionquad2back;P (GeV); #theta_x (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaYP_ionquad2back=new TH2F("acc_ThetaYP_ionquad2back","acc_ThetaYP_ionquad2back;P (GeV); #theta_y (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaXY_ionquad2back=new TH2F("acc_ThetaXY_ionquad2back","acc_ThetaXY_ionquad2back;#theta_x (deg); #theta_y (deg)",1000,-180.,180.,1000,-180.,180.);


    //    Ion QUADRUPOL 3
    TH2F *hacc_ThetaXP_ionquad3back=new TH2F("acc_ThetaXP_ionquad3back","acc_ThetaXP_ionquad3back;P (GeV); #theta_x (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaYP_ionquad3back=new TH2F("acc_ThetaYP_ionquad3back","acc_ThetaYP_ionquad3back;P (GeV); #theta_y (deg)",1000,0,100,1000,-180.,180.);
    TH2F *hacc_ThetaXY_ionquad3back=new TH2F("acc_ThetaXY_ionquad3back","acc_ThetaXY_ionquad3back;#theta_x (deg); #theta_y (deg)",1000,-180.,180.,1000,-180.,180.);

    TH2F *hgen_alpha_pRT=new TH2F("gen_alpha_pRT","gen_alpha_pRT;#alpha; P_{T}(GeV)",50,0.7,1.35,50,0,0.1);
    /* TH2F *hhit_alpha_pRT_iondipole2back=new TH2F("hit_alpha_pRT_iondipole2back","hit_alpha_pRT_iondipole2back;#alpha; P_{T}(GeV)",1000,0,2,1000,0,0.3); */
    TH2F *hhit_alpha_pRT_iondipole2back=new TH2F("hit_alpha_pRT_iondipole2back","hit_alpha_pRT_iondipole2back;#alpha; P_{T}(GeV)",50,0.74,1.36,50,0,0.1);

    TFile *file=new TFile(input_filename.c_str());
    if (file->IsZombie()) {
	cout << "Error opening file" << input_filename << endl;
	exit(-1);
    }
    else cout << "open file " << input_filename << endl;

    // variable in lund file and variable in header bank
    // nptl	ntgt	nproton	x	xinv_q2	alpha	xinv_pRT	tprime_KJ	crs1	dummy
    // ngen	var1	var2	var3	var4	var5	var6		var7		var8	var9

    TTree *header = (TTree*) file->Get("header");
    // int evn,evn_typ;
    // double beamPol;
    // double var1,var2,var3,var4,var5,var6,var7,var8,var9;
    vector <int> *header_evn=0,*header_evn_type=0;
    vector <double> *header_beamPol=0;
    /* vector <double> *header_var1=0,*header_var2=0,*header_var3=0,*header_var4=0,*header_var5=0,*header_var6=0,*header_var7=0,*header_var8=0,*header_var9=0; */
    vector <double> *header_var1=0,*header_var2=0,*header_var3=0,*header_var4=0,*header_var5=0,*header_var6=0,*header_var7=0,*header_var8=0;
    header->SetBranchAddress("evn",&header_evn);
    header->SetBranchAddress("evn_type",&header_evn_type);
    header->SetBranchAddress("beamPol",&header_beamPol);
    header->SetBranchAddress("var1",&header_var1); // xbj
    header->SetBranchAddress("var2",&header_var2); // Q2
    header->SetBranchAddress("var3",&header_var3); // xinv_se
    header->SetBranchAddress("var4",&header_var4); // alpha_R
    header->SetBranchAddress("var5",&header_var5); // xinv_pRT
    header->SetBranchAddress("var6",&header_var6); // xinv_tPrime
    header->SetBranchAddress("var7",&header_var7); // cross-section
    header->SetBranchAddress("var8",&header_var8); // ion beam pol

    TTree *generated = (TTree*) file->Get("generated");
    vector <double> *gen_pid=0,*gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
    generated->SetBranchAddress("pid",&gen_pid);
    generated->SetBranchAddress("px",&gen_px);
    generated->SetBranchAddress("py",&gen_py);
    generated->SetBranchAddress("pz",&gen_pz);
    generated->SetBranchAddress("vx",&gen_vx);
    generated->SetBranchAddress("vy",&gen_vy);
    generated->SetBranchAddress("vz",&gen_vz);

    TTree *flux = (TTree*) file->Get("flux");
    vector<int> *flux_id=0;
    vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0,*flux_avg_trackE=0,*flux_avg_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
    flux->SetBranchAddress("id",&flux_id);
    flux->SetBranchAddress("pid",&flux_pid);
    flux->SetBranchAddress("mpid",&flux_mpid);
    flux->SetBranchAddress("tid",&flux_tid);
    flux->SetBranchAddress("mtid",&flux_mtid);
    flux->SetBranchAddress("otid",&flux_otid);
    flux->SetBranchAddress("trackE",&flux_avg_trackE);
    flux->SetBranchAddress("totEdep",&flux_avg_totEdep);
    flux->SetBranchAddress("avg_x",&flux_avg_x);
    flux->SetBranchAddress("avg_y",&flux_avg_y);
    flux->SetBranchAddress("avg_z",&flux_avg_z);
    flux->SetBranchAddress("avg_lx",&flux_avg_lx);
    flux->SetBranchAddress("avg_ly",&flux_avg_ly);
    flux->SetBranchAddress("avg_lz",&flux_avg_lz);
    flux->SetBranchAddress("px",&flux_px);
    flux->SetBranchAddress("py",&flux_py);
    flux->SetBranchAddress("pz",&flux_pz);
    flux->SetBranchAddress("vx",&flux_vx);
    flux->SetBranchAddress("vy",&flux_vy);
    flux->SetBranchAddress("vz",&flux_vz);
    flux->SetBranchAddress("mvx",&flux_mvx);
    flux->SetBranchAddress("mvy",&flux_mvy);
    flux->SetBranchAddress("mvz",&flux_mvz);
    flux->SetBranchAddress("avg_t",&flux_avg_t);

    int nevent = (int)generated->GetEntries();
    // int nevent = (int)flux->GetEntries();
    int nselected = 0;
    cout << "nevent " << nevent << endl;

    int counter=0;

    double alpha,xinv_pRT,xinv_tPrime,crs;

    for (Long64_t i=0;i<nevent;i++) {
	//   cout << i << "\r";
	//   cout << i << "\n";

	int pid_gen;
	double p_gen,theta_gen,phi_gen,px_gen,py_gen,pz_gen,vx_gen,vy_gen,vz_gen;
	double t=0,Q2=0;
	double theta_x,theta_y;  
  
	header->GetEntry(i);
	for (int j=0;j<header_evn->size();j++) {
	    alpha       = header_var4->at(j);
	    xinv_pRT    = header_var5->at(j);
	    xinv_tPrime = header_var6->at(j);
	    crs         = header_var7->at(j);    
	}  
      
	generated->GetEntry(i);
	for (int j=0;j<gen_pid->size();j++) {

	    pid_gen=gen_pid->at(j);
	    px_gen=gen_px->at(j)/1e3;    	
	    py_gen=gen_py->at(j)/1e3;	
	    pz_gen=gen_pz->at(j)/1e3;      
	    vx_gen=gen_vx->at(j);    	
	    vy_gen=gen_vy->at(j);		
	    vz_gen=gen_vz->at(j);     

	    double mass=0;
	    if (pid_gen==2212) mass=0.938272;
	    else if (abs(pid_gen)==11) mass=0.5e-3;
	    else if (pid_gen==22) mass=0;
	    else if (pid_gen==2112) mass=93955;
	    else {cout << "unknown pid in gen " << pid_gen << endl;  return;}

	    p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);    	//in GeV
	    TLorentzVector p(px_gen,py_gen,pz_gen,sqrt(p_gen*p_gen+mass*mass));
	    theta_gen=p.Theta()*DEG;  	//in deg
	    phi_gen=p.Phi()*DEG;     	//in deg
	    //phi_gen=p.Theta()*DEG;     	//in deg

	    //      cout << " theta=" << theta_gen << " phi="<< phi_gen<< endl;

	    hgen_ThetaP->Fill(p_gen,theta_gen,crs);
	    hgen_ThetaPhi->Fill(phi_gen,theta_gen,crs);    
	    hgen_PhiP->Fill(p_gen,phi_gen,crs);
    
	    TLorentzVector p0(100*sin(50e-3),0,100*cos(50e-3),sqrt(100*100+0.938*0.938));
	    TLorentzVector p1(p_gen*sin(theta_gen/DEG)*cos(phi_gen/DEG),p_gen*sin(theta_gen/DEG)*sin(phi_gen/DEG),p_gen*cos(theta_gen/DEG),sqrt(p_gen*p_gen+0.938*0.938));
	    t=(p1-p0).M2();  //assume it's proton
	    Q2=-(p1-p0).M2(); //assume it's electron
    
    
	    hgen_ThetaP_ion->Fill(p_gen,theta_gen/DEG/1e-3,crs);
	    hgen_ThetaPhi_ion->Fill(phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);
	    hgen_PhiP_ion->Fill(p_gen,phi_gen/DEG/1e-3,crs);            
	    hgen_Thetat_ion->Fill(-t,theta_gen/DEG/1e-3,crs);      
	    hgen_ThetaPhiP_ion->Fill(p_gen,phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);
	 
	    TVector3 p3(sin(theta_gen/DEG)*cos(phi_gen/DEG),sin(theta_gen/DEG)*sin(phi_gen/DEG),cos(theta_gen/DEG));

	    p3.RotateY(-50e-3);
    

	    theta_x=atan(p3.X()/p3.Z());
	    theta_y=atan(p3.Y()/p3.Z());         
    
	    //	    cout << "test" << "momentum="<< p_gen << " thetaX=" << theta_x << " thetaY=" << theta_y << endl;

	    hgen_ThetaXP_ion->Fill(p_gen,theta_x*DEG,crs);
	    hgen_ThetaYP_ion->Fill(p_gen,theta_y*DEG,crs);
	    hgen_ThetaXY_ion->Fill(theta_x*DEG,theta_y*DEG,crs);  
    
	    hgen_alpha_pRT->Fill(alpha,xinv_pRT,1.);

	}  
  
	flux->GetEntry(i);    

	double px_flux,py_flux,pz_flux,theta_flux,p_flux,phi_flux;    
    
	int acc[2][2][2][8][2];
	for (int i1=0;i1<2;i1++){
	    for (int i2=0;i2<2;i2++){    
		for (int i3=0;i3<2;i3++){    
		    for (int i4=0;i4<8;i4++){    
			for (int i5=0;i5<2;i5++){    
			    acc[i1][i2][i3][i4][i5]=0;	      
			}}}}}
    
	for (int j=0;j<flux_id->size();j++) {
	    int beamline = int(flux_id->at(j))%100000/10000-1;
	    int side     = int(flux_id->at(j))%10000/1000-1;
	    int magnet   = int(flux_id->at(j))%1000/100-1;    
	    int number   = int(flux_id->at(j))%100/10-1;
	    int window   = int(flux_id->at(j))%10/1-1;      
	    acc[beamline][side][magnet][number][window]=1;
	}
    

	for (int j=0;j<flux_id->size();j++) {
	    // # == id name ==============================================================
	    // # digit     beamline                side          magnet      number          window
	    // #        ion    1        upstream    1    dipole    1           n     front      1
	    // #        ele    2       downstream   2  quadrupole  2                 back       2
	    // # =========================================================================
    
	    int beamline = int(flux_id->at(j))%100000/10000-1;
	    int side     = int(flux_id->at(j))%10000/1000-1;
	    int magnet   = int(flux_id->at(j))%1000/100-1;    
	    int number   = int(flux_id->at(j))%100/10-1;
	    int window   = int(flux_id->at(j))%10/1-1;
	    //     cout << flux_id->at(j) << " " << beamline << " " << side << " " << magnet << " " << number << " " << window << endl;

        
	    bool oktofill=false;
	    if (beamline==1-1 && side==2-1 && magnet==1-1 && number==1-1){   
		oktofill=true;
	    }
	    else if (beamline==1-1 && side==2-1 && magnet==2-1 && number==1-1){
		if (acc[1-1][2-1][1-1][1-1][2-1]==1) {oktofill=true; 
		} //quad 1             
	    }
	    else if (beamline==1-1 && side==2-1 && magnet==2-1 && number==2-1){
		if (acc[1-1][2-1][1-1][1-1][2-1]==1 && acc[1-1][2-1][2-1][1-1][2-1]==1) {oktofill=true;
		} //quad 2
	    }       
	    else if (beamline==1-1 && side==2-1 && magnet==2-1 && number==3-1){
		if (acc[1-1][2-1][1-1][1-1][2-1]==1 && acc[1-1][2-1][2-1][1-1][2-1]==1 && acc[1-1][2-1][2-1][2-1][2-1]==1) {oktofill=true;     	
		} //quad 3
	    }
	    else if (beamline==1-1 && side==2-1 && magnet==1-1 && number==2-1){
		if (acc[1-1][2-1][1-1][1-1][2-1]==1 && acc[1-1][2-1][2-1][1-1][2-1]==1 && acc[1-1][2-1][2-1][2-1][2-1]==1 && acc[1-1][2-1][2-1][3-1][2-1]==1) {
		    if(fabs(flux_avg_ly->at(j))<325) oktofill=true;
		    // 	 oktofill=true;
		} //dipole 2
	    }
	    else {
		oktofill=true; //everything else
	    }       
     
	    if (oktofill){
		hhit_alphaXY[beamline][side][magnet][number][window]->Fill(theta_x/1e-3,theta_y/1e-3,crs);	
		hhit_ThetaPhi[beamline][side][magnet][number][window]->Fill(phi_gen,theta_gen,crs); 
	    }    
     
	    //      	hhit_alphaXY[beamline][side][magnet][number][window]->Fill(theta_x/1e-3,theta_y/1e-3);	
	    // 	hhit_ThetaPhi[beamline][side][magnet][number][window]->Fill(phi_gen,theta_gen); 
     
    
	    px_flux=flux_px->at(j)/1e3;  //in MeV, convert to GeV
	    py_flux=flux_py->at(j)/1e3;	//in MeV, convert to GeV
	    pz_flux=flux_pz->at(j)/1e3; //in MeV, convert to GeV
	    p_flux=sqrt(px_flux*px_flux+py_flux*py_flux+pz_flux*pz_flux); 
	    theta_flux=acos(pz_flux/p_flux)*DEG;  //in deg
	    phi_flux=atan2(py_flux,px_flux)*DEG;  //in deg
    


	    if ( acc[0][1][1][0][1]==1 && flux_id->at(j) == 12212){      ///have to pass quad1_back
		hacc_ThetaP_ionquad1back->Fill(p_gen,theta_gen/DEG/1e-3,crs);    
		hacc_ThetaPhi_ionquad1back->Fill(phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);     
		hacc_PhiP_ionquad1back->Fill(p_gen,phi_gen/DEG/1e-3,crs);                 
		hacc_Thetat_ionquad1back->Fill(-t,theta_gen/DEG/1e-3,crs);         
      
		hhit_ThetaP_ionquad1back->Fill(p_flux,theta_flux/DEG/1e-3,crs);    
		hhit_ThetaPhi_ionquad1back->Fill(phi_flux/DEG/1e-3,theta_flux/DEG/1e-3,crs);     
		hhit_PhiP_ionquad1back->Fill(p_flux,phi_flux/DEG/1e-3,crs);                 
		hhit_Thetat_ionquad1back->Fill(-t,theta_flux/DEG/1e-3,crs);               
      
		hacc_ThetaXP_ionquad1back->Fill(p_gen,theta_x/DEG,crs);
		hacc_ThetaYP_ionquad1back->Fill(p_gen,theta_y/DEG,crs);
		hacc_ThetaXY_ionquad1back->Fill(theta_x/DEG,theta_y/DEG,crs);

	    }
    


	    if ( acc[0][1][1][0][1]==1 && flux_id->at(j) == 12222){      ///have to pass quad2_back
		hacc_ThetaP_ionquad2back->Fill(p_gen,theta_gen/DEG/1e-3,crs);    
		hacc_ThetaPhi_ionquad2back->Fill(phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);     
		hacc_PhiP_ionquad2back->Fill(p_gen,phi_gen/DEG/1e-3,crs);                 
		hacc_Thetat_ionquad2back->Fill(-t,theta_gen/DEG/1e-3,crs);         
      
		hhit_ThetaP_ionquad2back->Fill(p_flux,theta_flux/DEG/1e-3,crs);    
		hhit_ThetaPhi_ionquad2back->Fill(phi_flux/DEG/1e-3,theta_flux/DEG/1e-3,crs);     
		hhit_PhiP_ionquad2back->Fill(p_flux,phi_flux/DEG/1e-3,crs);                 
		hhit_Thetat_ionquad2back->Fill(-t,theta_flux/DEG/1e-3,crs);               
      
		hacc_ThetaXP_ionquad2back->Fill(p_gen,theta_x/DEG,crs);
		hacc_ThetaYP_ionquad2back->Fill(p_gen,theta_y/DEG,crs);
		hacc_ThetaXY_ionquad2back->Fill(theta_x/DEG,theta_y/DEG,crs);

	    }
    



	    if ( acc[0][1][1][0][1]==1 && flux_id->at(j) == 12232){      ///have to pass quad3_back
		hacc_ThetaP_ionquad3back->Fill(p_gen,theta_gen/DEG/1e-3,crs);    
		hacc_ThetaPhi_ionquad3back->Fill(phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);     
		hacc_PhiP_ionquad3back->Fill(p_gen,phi_gen/DEG/1e-3,crs);                 
		hacc_Thetat_ionquad3back->Fill(-t,theta_gen/DEG/1e-3,crs);         
      
		hhit_ThetaP_ionquad3back->Fill(p_flux,theta_flux/DEG/1e-3,crs);    
		hhit_ThetaPhi_ionquad3back->Fill(phi_flux/DEG/1e-3,theta_flux/DEG/1e-3,crs);     
		hhit_PhiP_ionquad3back->Fill(p_flux,phi_flux/DEG/1e-3,crs);                 
		hhit_Thetat_ionquad3back->Fill(-t,theta_flux/DEG/1e-3,crs);               
      
		hacc_ThetaXP_ionquad3back->Fill(p_gen,theta_x/DEG,crs);
		hacc_ThetaYP_ionquad3back->Fill(p_gen,theta_y/DEG,crs);
		hacc_ThetaXY_ionquad3back->Fill(theta_x/DEG,theta_y/DEG,crs);
		//       cout << theta << " " << theta_gen/DEG/1e-3 << endl;
	    }
    




	    if ( acc[0][1][1][0][1]==1 && flux_id->at(j) == 12122){ ///have to pass ion dipol2
		hacc_ThetaP_iondipole2back->Fill(p_gen,theta_gen/DEG/1e-3,crs);    
		hacc_ThetaPhi_iondipole2back->Fill(phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);     
		hacc_PhiP_iondipole2back->Fill(p_gen,phi_gen/DEG/1e-3,crs);                 
		hacc_Thetat_iondipole2back->Fill(-t,theta_gen/DEG/1e-3,crs);         
	
		hacc_ThetaPhiP_iondipole2back->Fill(p_gen,phi_gen/DEG/1e-3,theta_gen/DEG/1e-3,crs);      
      
		hhit_ThetaP_iondipole2back->Fill(p_flux,theta_flux/DEG/1e-3,crs);    
		hhit_ThetaPhi_iondipole2back->Fill(phi_flux/DEG/1e-3,theta_flux/DEG/1e-3,crs);     
		hhit_PhiP_iondipole2back->Fill(p_flux,phi_flux/DEG/1e-3,crs);                 
		hhit_Thetat_iondipole2back->Fill(-t,theta_flux/DEG/1e-3,crs);          
      
		hacc_ThetaXP_iondipole2back->Fill(p_gen,theta_x/DEG,crs);
		hacc_ThetaYP_iondipole2back->Fill(p_gen,theta_y/DEG,crs);
		hacc_ThetaXY_iondipole2back->Fill(theta_x/DEG,theta_y/DEG,crs);  
      
		/* if(fabs(flux_avg_ly->at(j))<325) hhit_alpha_pRT_iondipole2back->Fill(alpha,xinv_pRT,crs); */
		/* hhit_alpha_pRT_iondipole2back->Fill(alpha,xinv_pRT,1.); */
		hhit_alpha_pRT_iondipole2back->Fill(alpha,xinv_pRT,crs);

		counter++;
       
	    }    
    
    
	}
    
	// cout << " theta=" <<theta_gen << " phi="<< phi_gen<< endl;
    
    }
    file->Close();

    cout << "counter " << counter <<  endl;

    /* TCanvas *c1 = new TCanvas("c1","c1",1800,800); */
    /* c1->Divide(2,1); */
    /* c1->cd(1); */
    /* hgen_alpha_pRT->Draw("colz"); */
    /* c1->cd(2); */
    /* hhit_alpha_pRT_iondipole2back->Draw("colz"); */
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->Divide(1,1);
    c1->cd(1);
    hhit_alpha_pRT_iondipole2back->Draw("colz");
    TCanvas *c2 = new TCanvas("c2","c2",800,800);
    c2->Divide(1,1);
    c2->cd(1);
    hgen_alpha_pRT->Draw("colz");

    // # == id name ==============================================================
    // # digit     beamline                side          magnet      number          window
    // #        ion    1        upstream    1    dipole    1           n     front      1
    // #        ele    2       downstream   2  quadrupole  2                 back       2
    // # =========================================================================

    /* TCanvas *c_gen_iondipole2back_3D = new TCanvas("gen_iondipole2back_3D","gen_iondipole2back_3D",600,600); */
    /* hgen_ThetaPhiP_ion->Draw("box"); */

    /* TCanvas *c_acc_iondipole2back_3D = new TCanvas("acc_iondipole2back_3D","acc_iondipole2back_3D",600,600); */
    /* hacc_ThetaPhiP_iondipole2back->Divide(hgen_ThetaPhiP_ion); */
    /* hacc_ThetaPhiP_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaPhiP_iondipole2back->Draw("box"); */




    /* TCanvas *c_hit_alphaXY_ion = new TCanvas("hit_alphaXY_ion","hit_alphaXY_ion",900,900); */
    /* hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->SetMarkerColor(kOrange); */
    /* hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->Draw(""); */
    /* hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->SetMarkerColor(kGreen); */
    /* hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->Draw("same"); */
    /* hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->SetMarkerColor(kCyan); */
    /* hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->Draw("same"); */
    /* hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->SetMarkerColor(kViolet); */
    /* hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->Draw("same"); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->SetMarkerColor(kYellow); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->Draw("same"); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->SetMarkerColor(kBlue); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->SetMarkerStyle(kFullSquare); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->SetMarkerSize(0.5); */
    /* cout << hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->Integral() << endl; */
    /* // hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->SetContour(1); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->Draw("same"); */

    /* TCanvas *c_hit_alphaXY_ion_all = new TCanvas("hit_alphaXY_ion_all","hit_alphaXY_ion_all",900,900); */
    /* c_hit_alphaXY_ion_all->Divide(3,2); */
    /* c_hit_alphaXY_ion_all->cd(1); */
    /* hhit_alphaXY[1-1][2-1][1-1][1-1][2-1]->Draw(); */
    /* c_hit_alphaXY_ion_all->cd(2); */
    /* hhit_alphaXY[1-1][2-1][2-1][1-1][2-1]->Draw("same"); */
    /* c_hit_alphaXY_ion_all->cd(3); */
    /* hhit_alphaXY[1-1][2-1][2-1][2-1][2-1]->Draw("same"); */
    /* c_hit_alphaXY_ion_all->cd(4); */
    /* hhit_alphaXY[1-1][2-1][2-1][3-1][2-1]->Draw("same"); */
    /* c_hit_alphaXY_ion_all->cd(5); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][1-1]->Draw("same"); */
    /* c_hit_alphaXY_ion_all->cd(6); */
    /* hhit_alphaXY[1-1][2-1][1-1][2-1][2-1]->Draw("same"); */

    /* TCanvas *c_hit_ThetaPhi_ion_all = new TCanvas("hit_ThetaPhi_ion_all","hit_ThetaPhi_ion_all",900,900); */
    /* c_hit_ThetaPhi_ion_all->Divide(3,2); */
    /* c_hit_ThetaPhi_ion_all->cd(1); */
    /* hhit_ThetaPhi[1-1][2-1][1-1][1-1][2-1]->Draw(); */
    /* c_hit_ThetaPhi_ion_all->cd(2); */
    /* hhit_ThetaPhi[1-1][2-1][2-1][1-1][2-1]->Draw("same"); */
    /* c_hit_ThetaPhi_ion_all->cd(3); */
    /* hhit_ThetaPhi[1-1][2-1][2-1][2-1][2-1]->Draw("same"); */
    /* c_hit_ThetaPhi_ion_all->cd(4); */
    /* hhit_ThetaPhi[1-1][2-1][2-1][3-1][2-1]->Draw("same"); */
    /* c_hit_ThetaPhi_ion_all->cd(5); */
    /* hhit_ThetaPhi[1-1][2-1][1-1][2-1][1-1]->Draw("same"); */
    /* c_hit_ThetaPhi_ion_all->cd(6); */
    /* hhit_ThetaPhi[1-1][2-1][1-1][2-1][2-1]->Draw("same"); */

    /* TCanvas *c_gen = new TCanvas("gen","gen",1500,900); */
    /* c_gen->Divide(2,2); */
    /* c_gen->cd(1); */
    /* hgen_ThetaP->Draw("colz"); */
    /* c_gen->cd(2); */
    /* hgen_ThetaPhi->Draw("colz"); */
    /* c_gen->cd(3); */
    /* hgen_PhiP->Draw("colz"); */

    /* // hhit_ThetaP[5][0]->Add(hhit_ThetaP[4][0],hhit_ThetaP[4][1],1,-1); */
    /* // hhit_ThetaP[5][1]->Add(hhit_ThetaP[4][0],hhit_ThetaP[4][1],1,-1); */

    /* TCanvas *c_gen_ion = new TCanvas("gen_ion","gen_ion",1800,500); */
    /* c_gen_ion->Divide(3,1); */
    /* c_gen_ion->cd(1); */
    /* hgen_ThetaP_ion->Draw("colz"); */
    /* c_gen_ion->cd(2); */
    /* hgen_ThetaPhi_ion->Draw("colz"); */
    /* c_gen_ion->cd(3); */
    /* hgen_PhiP_ion->Draw("colz"); */









    /* TCanvas *c_acc_ionquad1back = new TCanvas("acc_ionquad1back","acc_ionquad1back",1800,500); */
    /* c_acc_ionquad1back->Divide(3,1); */
    /* c_acc_ionquad1back->cd(1); */
    /* hacc_ThetaP_ionquad1back->Divide(hgen_ThetaP_ion); */
    /* hacc_ThetaP_ionquad1back->SetMaximum(1);  */
    /* hacc_ThetaP_ionquad1back->Draw("colz"); */
    /* c_acc_ionquad1back->cd(2); */
    /* hacc_ThetaPhi_ionquad1back->Divide(hgen_ThetaPhi_ion); */
    /* hacc_ThetaPhi_ionquad1back->SetMaximum(1); */
    /* hacc_ThetaPhi_ionquad1back->Draw("colz"); */
    /* c_acc_ionquad1back->cd(3); */
    /* hacc_PhiP_ionquad1back->Divide(hgen_PhiP_ion); */
    /* hacc_PhiP_ionquad1back->SetMaximum(1); */
    /* hacc_PhiP_ionquad1back->Draw("colz"); */

    /* TCanvas *c_hit_ionquad1back = new TCanvas("hit_ionquad1back","hit_ionquad1back",1800,500); */
    /* c_hit_ionquad1back->Divide(3,1); */
    /* c_hit_ionquad1back->cd(1); */
    /* // hhit_ThetaP_ionquad1back->Divide(hgen_ThetaP_ion); */
    /* // hhit_ThetaP_ionquad1back->SetMaximum(1); */
    /* hhit_ThetaP_ionquad1back->Draw("colz"); */
    /* c_hit_ionquad1back->cd(2); */
    /* // hhit_ThetaPhi_ionquad1back->Divide(hgen_ThetaPhi_ion); */
    /* // hhit_ThetaPhi_ionquad1back->SetMaximum(1); */
    /* hhit_ThetaPhi_ionquad1back->Draw("colz"); */
    /* c_hit_ionquad1back->cd(3); */
    /* // hhit_PhiP_ionquad1back->Divide(hgen_PhiP_ion); */
    /* // hhit_PhiP_ionquad1back->SetMaximum(1); */
    /* hhit_PhiP_ionquad1back->Draw("colz"); */












    /* TCanvas *c_acc_ionquad2back = new TCanvas("acc_ionquad2back","acc_ionquad2back",1800,500); */
    /* c_acc_ionquad2back->Divide(3,1); */
    /* c_acc_ionquad2back->cd(1); */
    /* hacc_ThetaP_ionquad2back->Divide(hgen_ThetaP_ion); */
    /* hacc_ThetaP_ionquad2back->SetMaximum(1); */
    /* hacc_ThetaP_ionquad2back->Draw("colz"); */
    /* c_acc_ionquad2back->cd(2); */
    /* hacc_ThetaPhi_ionquad2back->Divide(hgen_ThetaPhi_ion); */
    /* hacc_ThetaPhi_ionquad2back->SetMaximum(1); */
    /* hacc_ThetaPhi_ionquad2back->Draw("colz"); */
    /* c_acc_ionquad2back->cd(3); */
    /* hacc_PhiP_ionquad2back->Divide(hgen_PhiP_ion); */
    /* hacc_PhiP_ionquad2back->SetMaximum(1); */
    /* hacc_PhiP_ionquad2back->Draw("colz"); */

    /* TCanvas *c_hit_ionquad2back = new TCanvas("hit_ionquad2back","hit_ionquad2back",1800,500); */
    /* c_hit_ionquad2back->Divide(3,1); */
    /* c_hit_ionquad2back->cd(1); */
    /* // hhit_ThetaP_ionquad2back->Divide(hgen_ThetaP_ion); */
    /* // hhit_ThetaP_ionquad2back->SetMaximum(1); */
    /* hhit_ThetaP_ionquad2back->Draw("colz"); */
    /* c_hit_ionquad2back->cd(2); */
    /* // hhit_ThetaPhi_ionquad2back->Divide(hgen_ThetaPhi_ion); */
    /* // hhit_ThetaPhi_ionquad2back->SetMaximum(1); */
    /* hhit_ThetaPhi_ionquad2back->Draw("colz"); */
    /* c_hit_ionquad2back->cd(3); */
    /* // hhit_PhiP_ionquad2back->Divide(hgen_PhiP_ion); */
    /* // hhit_PhiP_ionquad2back->SetMaximum(1); */
    /* hhit_PhiP_ionquad2back->Draw("colz"); */


















    /* TCanvas *c_acc_ionquad3back = new TCanvas("acc_ionquad3back","acc_ionquad3back",1800,500); */
    /* c_acc_ionquad3back->Divide(3,1); */
    /* c_acc_ionquad3back->cd(1); */
    /* hacc_ThetaP_ionquad3back->Divide(hgen_ThetaP_ion); */
    /* hacc_ThetaP_ionquad3back->SetMaximum(1); */
    /* hacc_ThetaP_ionquad3back->Draw("colz"); */
    /* c_acc_ionquad3back->cd(2); */
    /* hacc_ThetaPhi_ionquad3back->Divide(hgen_ThetaPhi_ion); */
    /* hacc_ThetaPhi_ionquad3back->SetMaximum(1); */
    /* hacc_ThetaPhi_ionquad3back->Draw("colz"); */
    /* c_acc_ionquad3back->cd(3); */
    /* hacc_PhiP_ionquad3back->Divide(hgen_PhiP_ion); */
    /* hacc_PhiP_ionquad3back->SetMaximum(1); */
    /* hacc_PhiP_ionquad3back->Draw("colz"); */

    /* TCanvas *c_hit_ionquad3back = new TCanvas("hit_ionquad3back","hit_ionquad3back",1800,500); */
    /* c_hit_ionquad3back->Divide(3,1); */
    /* c_hit_ionquad3back->cd(1); */
    /* // hhit_ThetaP_ionquad3back->Divide(hgen_ThetaP_ion); */
    /* // hhit_ThetaP_ionquad3back->SetMaximum(1); */
    /* hhit_ThetaP_ionquad3back->Draw("colz"); */
    /* c_hit_ionquad3back->cd(2); */
    /* // hhit_ThetaPhi_ionquad3back->Divide(hgen_ThetaPhi_ion); */
    /* // hhit_ThetaPhi_ionquad3back->SetMaximum(1); */
    /* hhit_ThetaPhi_ionquad3back->Draw("colz"); */
    /* c_hit_ionquad3back->cd(3); */
    /* // hhit_PhiP_ionquad3back->Divide(hgen_PhiP_ion); */
    /* // hhit_PhiP_ionquad3back->SetMaximum(1); */
    /* hhit_PhiP_ionquad3back->Draw("colz"); */













    /* TCanvas *c_acc_iondipole2back = new TCanvas("acc_iondipole2back","acc_iondipole2back",1800,500); */
    /* c_acc_iondipole2back->Divide(3,1); */
    /* c_acc_iondipole2back->cd(1); */
    /* hacc_ThetaP_iondipole2back->Divide(hgen_ThetaP_ion); */
    /* hacc_ThetaP_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaP_iondipole2back->Draw("colz"); */
    /* c_acc_iondipole2back->cd(2); */
    /* hacc_ThetaPhi_iondipole2back->Divide(hgen_ThetaPhi_ion); */
    /* hacc_ThetaPhi_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaPhi_iondipole2back->Draw("colz"); */
    /* c_acc_iondipole2back->cd(3); */
    /* hacc_PhiP_iondipole2back->Divide(hgen_PhiP_ion); */
    /* hacc_PhiP_iondipole2back->SetMaximum(1); */
    /* hacc_PhiP_iondipole2back->Draw("colz"); */

    /* TCanvas *c_hit_iondipole2back = new TCanvas("hit_iondipole2back","hit_iondipole2back",1800,500); */
    /* c_hit_iondipole2back->Divide(3,1); */
    /* c_hit_iondipole2back->cd(1); */
    /* // hhit_ThetaP_iondipole2back->Divide(hgen_ThetaP_ion); */
    /* // hhit_ThetaP_iondipole2back->SetMaximum(1); */
    /* hhit_ThetaP_iondipole2back->Draw("colz"); */
    /* c_hit_iondipole2back->cd(2); */
    /* // hhit_ThetaPhi_iondipole2back->Divide(hgen_ThetaPhi_ion); */
    /* // hhit_ThetaPhi_iondipole2back->SetMaximum(1); */
    /* hhit_ThetaPhi_iondipole2back->Draw("colz"); */
    /* c_hit_iondipole2back->cd(3); */
    /* // hhit_PhiP_iondipole2back->Divide(hgen_PhiP_ion); */
    /* // hhit_PhiP_iondipole2back->SetMaximum(1); */
    /* hhit_PhiP_iondipole2back->Draw("colz"); */



    /* TCanvas *c_gen_ion_ACC = new TCanvas("gen_ion_ACC","gen_ion_ACC",1800,500); */
    /* c_gen_ion_ACC->Divide(3,1); */
    /* c_gen_ion_ACC->cd(1); */
    /* hgen_ThetaXP_ion->Draw("colz"); */
    /* c_gen_ion_ACC->cd(2); */
    /* hgen_ThetaYP_ion->Draw("colz"); */
    /* c_gen_ion_ACC->cd(3); */
    /* hgen_ThetaXY_ion->Draw("colz"); */

    /* TCanvas *c_acc_ionquad3back_ACC = new TCanvas("acc_ionquad3back_ACC","acc_ionquad3back_ACC",1800,500); */
    /* c_acc_ionquad3back_ACC->Divide(3,1); */
    /* c_acc_ionquad3back_ACC->cd(1); */
    /* hacc_ThetaXP_ionquad3back->Divide(hgen_ThetaXP_ion); */
    /* hacc_ThetaXP_ionquad3back->SetMaximum(1); */
    /* hacc_ThetaXP_ionquad3back->Draw("colz"); */
    /* c_acc_ionquad3back_ACC->cd(2); */
    /* hacc_ThetaYP_ionquad3back->Divide(hgen_ThetaYP_ion); */
    /* hacc_ThetaYP_ionquad3back->SetMaximum(1); */
    /* hacc_ThetaYP_ionquad3back->Draw("colz"); */
    /* c_acc_ionquad3back_ACC->cd(3); */
    /* hacc_ThetaXY_ionquad3back->Divide(hgen_ThetaXY_ion); */
    /* hacc_ThetaXY_ionquad3back->SetMaximum(1); */
    /* hacc_ThetaXY_ionquad3back->Draw("colz"); */

    /* TCanvas *c_acc_iondipole2back_ACC = new TCanvas("acc_iondipole2back_ACC","acc_iondipole2back_ACC",1800,500); */
    /* c_acc_iondipole2back_ACC->Divide(3,1); */
    /* c_acc_iondipole2back_ACC->cd(1); */
    /* hacc_ThetaXP_iondipole2back->Divide(hgen_ThetaXP_ion); */
    /* hacc_ThetaXP_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaXP_iondipole2back->Draw("colz"); */
    /* c_acc_iondipole2back_ACC->cd(2); */
    /* hacc_ThetaYP_iondipole2back->Divide(hgen_ThetaYP_ion); */
    /* hacc_ThetaYP_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaYP_iondipole2back->Draw("colz"); */
    /* c_acc_iondipole2back_ACC->cd(3); */
    /* hacc_ThetaXY_iondipole2back->Divide(hgen_ThetaXY_ion); */
    /* hacc_ThetaXY_iondipole2back->SetMaximum(1); */
    /* hacc_ThetaXY_iondipole2back->Draw("colz"); */

    outputfile->Write();
    outputfile->Flush();
}
