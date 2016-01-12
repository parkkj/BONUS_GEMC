// Code to train to write a Kalman filter
// http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf by Welch and Bishop
//
// Now I try to create fit of a particle in a magnetic field
// Its (x,y,z) coordinates are measuread each Dt
// There is no error
// 
// 
//
// Author: Gabriel Charles
// 2015/09/28
// modification for BONUS12 simulation : K.Park Dec./2015

#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMultiGraph.h"
#include <iostream>
#include <ctime>
#include "BonusHelixFit.h"
#include "TMath.h"
#include "TMatrix.h"




//____________________________________________________________________
void KFilter_KJ1(){


//________________________________________________________________________________________________
// ____________________________________________Variables__________________________________________
//________________________________________________________________________________________________

   const int NHit = 6e4; // 6e5
   double dt = 1e-13; // a measurement is taken each Dt, s 1e-14

   Bool_t iran=kTRUE;      // TRUE==include with momentum spread in x,y
   //Bool_t iran=kFALSE;      // TRUE==include with momentum spread in x,y
  

   const double MeVcToSi = 5.34428595e-22; // constant to change MeV/c in standard units (SI unit: kg*m/s)
   const double MeVc2ToSi = 1.783e-30; // constant to change MeV/c² in standard units
   const double elemcharge = 1.60217662e-19; // elementary charge

   double sigMomloss = 4e-5 * MeVcToSi; // Eloss per step (added by hand nothing physical here: 40eV?? K.Park)
   //double sigMomloss = 1e-2 * MeVcToSi; // Eloss per step 10keV loss from GEANT simulation
   
   double Momloss; // temporary value that will contain the instantaneous Eloss

   double MeasX[NHit];
   double MeasY[NHit];
   double MeasZ[NHit];

 // The same but with random error following a gaussian distribution added
   double MeasXr[NHit];
   double MeasYr[NHit];
   double MeasZr[NHit];

   double MeasPX[NHit];
   double MeasPY[NHit];
   double MeasPZ[NHit];

   double pData[NHit][3];

   double X_gen = 0.; // initial position of the object, in m
   double Y_gen = 0.; // initial position of the object, in m
   double Z_gen = -30e-3; // initial position of the object, in m


   
   double PBeamMC = 120.; // input of spectator proton momentum 
   
   double PX_gen,PY_gen,PZ_gen;
   PX_gen=0;PY_gen=0;PZ_gen=0;
    // initial incoming proton momentum information
   // make sure proton can reach GEM > 72MeV



   // output variable definition : reconstructed Transverse momentum from Helix and KF fits.
   double pTresol_helix=0;
   double pXresol_helix=0;
   double pYresol_helix=0;
   
   double pTresol_kf=0;
   double pXresol_kf=0;
   double pYresol_kf=0;
   
   // output variable definition : reconstructed r, theta, phi from Helix and KF fits.
   double rDiff_helix =0;
   double thetaDiff_helix =0;
   double phiDiff_helix =0;
   double zDiff_helix =0;

   double rDiff_kf =0;
   double thetaDiff_kf =0;
   double phiDiff_kf =0;
   double zDiff_kf =0;

   
   //output file : ASCII file for Helix fit result and Kalman Filter fit result;
   
   ofstream OUT1 ("helix_out.txt", ios::app);
   ofstream OUT2 ("kfilt_out.txt", ios::app);
   
   //output file: ROOT
   TFile *f1 = new TFile("helix.root","RECREATE");
   TNtuple *ntuple1 = new TNtuple("ntuple1","ASCII data from Helix fit","pTresol_helix:pXresol_helix:pYresol_helix:rDiff_helix:thetaDiff_helix:phiDiff_helix:zDiff_helix");

   TFile *f2 = new TFile("kfilt.root","RECREATE");
   TNtuple *ntuple2 = new TNtuple("ntuple2","ASCII data from KF fit","pTresol_kf:pXresol_kf:pYresol_kf:rDiff_kf:thetaDiff_kf:phiDiff_kf:zDiff_kf");

   
   int NEvts = 1000;
     
   for (int iEvt=0; iEvt<NEvts; iEvt++) { // loop for Nevent 

     // random number in each event
   time_t now = time(0);
   TRandom3 *Vrand = new TRandom3();
   Vrand->SetSeed(now);
   //Vrand->SetSeed(1); // fixed random seed for debugging for various covariant error matrix terms with error RMS
   
   srand (time(NULL));
   double rnumber = rand() ;
   TRandom3 ran3;
   // ran3.SetSeed(13579);
   ran3.SetSeed(rnumber);

   
     
     if(iran){
       double PBeamMCx= PBeamMC*ran3.Gaus(0.71,0.001); // mean 0.71 for 70MeV of proton, 0.1% deviation in x
       double PBeamMCy= PBeamMC*ran3.Gaus(0.71,0.001); // mean 0.71 for 70MeV of proton, 0.1% deviation in y
       double PBeamMCz= PBeamMC*ran3.Gaus(0.001,0.001);
     }
     else{
       double PBeamMCx= PBeamMC/TMath::Sqrt(3);
       double PBeamMCy= PBeamMC/TMath::Sqrt(3);
       double PBeamMCz= PBeamMC/TMath::Sqrt(3);    
     }
     
     PX_gen = MeVcToSi*PBeamMCx; 
     PY_gen = MeVcToSi*PBeamMCy; 
     PZ_gen = MeVcToSi*PBeamMCz; 
    

     //     cout << "X=" << PBeamMCx  << " Y=" << PBeamMCy  << " Z=" << PBeamMCz  << endl;
     cout << "--------------------------------------------------------------------------------> " << iEvt+1 << ", rand " << rnumber << endl;
   
   double pT_gen = TMath::Sqrt(PX_gen*PX_gen+PY_gen*PY_gen);
   double phi_gen = TMath::ATan2(PY_gen,PX_gen);
      if(phi_gen<0) phi_gen+=2.*TMath::Pi();
   double theta_gen = TMath::ATan(pT_gen/PZ_gen);
      if(theta_gen<0) theta_gen+=TMath::Pi();
// 0.3 m en 1 ns
   double B = 5.; // Magnetic field in Tesla

   double mymass = 938.*MeVc2ToSi; // mass of the particle 
//   double mymass = 0.511*MeVc2ToSi; // mass of the particle
//cout << mymass << endl;
   double mycharge = 1.*elemcharge;
   double dR = 1e-4; // for now set by hand but should be larger with true or simulated data

   double DistTot=0.;

   double Dt, Vinst;

   int ind_pData=0;
   int ind_kal=0;

   double phi_kal, theta_kal;
   double alpha, beta, gamma; // variables used to simplify the Kalman filter

// Simulationtion parameters (Are these precision of pseudo data??, how we can check ?? K.Park)
   double sigX = 0.3e-3;
   double sigY = 0.3e-3;
   double sigZ = 0.5e-3;
 
// Reconstruction parameters (how we can define ?? K.Park)
   double sigX_Kf = 1e-3;
   double sigY_Kf = 1e-3;
   double sigZ_Kf = 1e-3;
   TMatrix R(3,3);   // measurement noise covariance
      R(0,0) = sigX_Kf*sigX_Kf; R(0,1) = 0.;      	  R(0,2) = 0.;
      R(1,0) = 0.;      	R(1,1) = sigY_Kf*sigY_Kf; R(1,2) = 0.;
      R(2,0) = 0.;      	R(2,1) = 0.;      	  R(2,2) = sigZ_Kf*sigZ_Kf;

   double var = 1e-27; // to give freedom to tune the model as for now there is no process noise 1e-21
   double covar = 1e-27;
   TMatrix Q(6,6);   // process noise covariance
      Q(0,0) = var;    Q(0,1) = 0.;      Q(0,2) = 0.;      Q(0,3) = covar;   Q(0,4) = 0.;   Q(0,5) = 0.;
      Q(1,0) = 0.;     Q(1,1) = var;     Q(1,2) = 0.;      Q(1,3) = 0.;   Q(1,4) = 0.;   Q(1,5) = 0.;
      Q(2,0) = 0.;     Q(2,1) = 0.;      Q(2,2) = var;     Q(2,3) = 0.;   Q(2,4) = 0.;   Q(2,5) = 0.;
      Q(3,0) = -var;   Q(3,1) = -covar;  Q(3,2) = -covar;  Q(3,3) = var;  Q(3,4) = 0.;   Q(3,5) = 0.;
      Q(4,0) = -covar; Q(4,1) = -var;    Q(4,2) = -covar;  Q(4,3) = 0.;   Q(4,4) = var;  Q(4,5) = 0.;
      Q(5,0) = 0.;     Q(5,1) = 0.;      Q(5,2) = -covar;  Q(5,3) = 0.;   Q(5,4) = 0.;   Q(5,5) = var;  

   TMatrix A(6,6); TMatrix AT(6,6);
      A(0,0) = 1.; A(0,1) = 0.; A(0,2) = 0.; A(0,5) = 0.; 
      A(1,0) = 0.; A(1,1) = 1.; A(1,2) = 0.; A(1,5) = 0.;
      A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 1.; A(2,3) = 0.; A(2,4) = 0.;
      A(3,0) = 0.; A(3,1) = 0.; A(3,2) = 0.; A(3,5) = 0.; 
      A(4,0) = 0.; A(4,1) = 0.; A(4,2) = 0.; A(4,5) = 0.; 
      A(5,0) = 0.; A(5,1) = 0.; A(5,2) = 0.; A(5,3) = 0.; A(5,4) = 0.;  A(5,5) = 1.;
      AT.Transpose(A); 

   TMatrix H(3,6); TMatrix HT(6,3);
     H(0,0) = 1; H(0,1) = 0; H(0,2) = 0; H(0,3) = 0; H(0,4) = 0; H(0,5) = 0;
     H(1,0) = 0; H(1,1) = 1; H(1,2) = 0; H(1,3) = 0; H(1,4) = 0; H(1,5) = 0;
     H(2,0) = 0; H(2,1) = 0; H(2,2) = 1; H(2,3) = 0; H(1,4) = 0; H(2,5) = 0;    
     HT.Transpose(H);
   TMatrix X0(6,1);  // Measurement matrix of the first point
   TMatrix P0(6,6);  // error on the position for the first point

   TMatrix Id(6,6);
      Id(0,0) = 1.;  Id(0,1) = 0.;  Id(0,2) = 0.;  Id(0,3) = 0.;  Id(0,4) = 0.;  Id(0,5) = 0.;
      Id(1,0) = 0.;  Id(1,1) = 1.;  Id(1,2) = 0.;  Id(1,3) = 0.;  Id(1,4) = 0.;  Id(1,5) = 0.;
      Id(2,0) = 0.;  Id(2,1) = 0.;  Id(2,2) = 1.;  Id(2,3) = 0.;  Id(2,4) = 0.;  Id(2,5) = 0.;
      Id(3,0) = 0.;  Id(3,1) = 0.;  Id(3,2) = 0.;  Id(3,3) = 1.;  Id(3,4) = 0.;  Id(3,5) = 0.;
      Id(4,0) = 0.;  Id(4,1) = 0.;  Id(4,2) = 0.;  Id(4,3) = 0.;  Id(4,4) = 1.;  Id(4,5) = 0.;
      Id(5,0) = 0.;  Id(5,1) = 0.;  Id(5,2) = 0.;  Id(5,3) = 0.;  Id(5,4) = 0.;  Id(5,5) = 1.; 

   TMatrix Temp(3,1);

   TMatrix Xkm1(6,1); // quantity of interest at step k-1
   TMatrix Pkm1(6,6); // covariance matrix of the error on the above quantity
   TMatrix Xkm(6,1);  // predicted quantity of interest at step k
   TMatrix Pkm(6,6);  // predicted covariance matrix of the error on the above quantity
   TMatrix Kk(6,3);   // Kalman gain at step k
   TMatrix Xk(6,1);   // calculated quantity of interest at step k
   TMatrix Pk(6,6);   // calculated covariance matrix of the error on the above quantity at step k


//________________________________________________________________________________________________
//_____________________________________ Canvas and histograms ____________________________________
//________________________________________________________________________________________________
/*
   TCanvas *c0 = new TCanvas("c0","View of the Kalman filter effect",500,400);
   TCanvas *cX = new TCanvas("cX","View of the Kalman filter for X",500,400);
   TCanvas *cY = new TCanvas("cY","View of the Kalman filter for Y",500,400);
   TCanvas *cZ = new TCanvas("cZ","View of the Kalman filter for Z",500,400);
   TCanvas *cXY = new TCanvas("cxy","View of the Kalman filter for xy",500,400);
   TCanvas *cdR = new TCanvas("cdR","R growth",500,400);

   TCanvas *cPX = new TCanvas("cPX","View of the Kalman filter for PX",500,400);
   TCanvas *cPY = new TCanvas("cPY","View of the Kalman filter for PY",500,400);
   TCanvas *cPZ = new TCanvas("cPZ","View of the Kalman filter for PZ",500,400);
   TCanvas *cPT = new TCanvas("cPT","View of the Kalman filter for PT",500,400);

   TMultiGraph *mg_X = new TMultiGraph();
   TMultiGraph *mg_Y = new TMultiGraph();
   TMultiGraph *mg_Z = new TMultiGraph();
   TMultiGraph *mg_XY = new TMultiGraph();

   TMultiGraph *mg_PX = new TMultiGraph();
   TMultiGraph *mg_PY = new TMultiGraph();
   TMultiGraph *mg_PZ = new TMultiGraph();
   TMultiGraph *mg_PT = new TMultiGraph();
*/
   TGraph2D *gr_P = new TGraph2D();
   TGraph2D *gr_P_Kal = new TGraph2D();
   TGraph *gr_X = new TGraph();
   TGraph *gr_X_Kal = new TGraph();
   TGraph *gr_Y = new TGraph();
   TGraph *gr_Y_Kal = new TGraph();
   TGraph *gr_Z = new TGraph();
   TGraph *gr_Z_Kal = new TGraph();
   TGraph *gr_XY = new TGraph();
   TGraph *gr_XY_Kal = new TGraph();

   TGraph *gr_PX = new TGraph();
   TGraph *gr_PX_Kal = new TGraph();
   TGraph *gr_PY = new TGraph();
   TGraph *gr_PY_Kal = new TGraph();
   TGraph *gr_PZ = new TGraph();
   TGraph *gr_PZ_Kal = new TGraph();
   TGraph *gr_PT = new TGraph();
   TGraph *gr_PT_Kal = new TGraph();

   //   TH1D *h_Rgrowth = new TH1D("h_Rgrowth","R growth",100,1e-9,1e-5);

   
//________________________________________________________________________________________________
// ________________________________________ Event loop ___________________________________________
//________________________________________________________________________________________________

   // Generates values following a gaussian of mean mean and sigma std
      MeasX[0] = X_gen; // Initial X position of the object
         MeasXr[0] = Vrand->Gaus(MeasX[0],sigX);
      MeasY[0] = Y_gen; // Initial Y position of the object
         MeasYr[0] = Vrand->Gaus(MeasY[0],sigY);
      MeasZ[0] = Z_gen; // Initial Z position of the object
         MeasZr[0] = Vrand->Gaus(MeasZ[0],sigZ);

   // Generates values following a gaussian of mean mean and sigma std
      MeasPX[0] = PX_gen; // Initial X momentum of the object
      MeasPY[0] = PY_gen; // Initial Y momentum of the object
      MeasPZ[0] = PZ_gen; // Initial Z momentum of the object
 cout << "Momentum before mom loss " << TMath::Sqrt(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0])/MeVcToSi << endl;
         Momloss = -1.;
         while(Momloss<0){ // energy (in fact momentum) loss
            Momloss = Vrand->Gaus(0.,sigMomloss);
         }
      MeasPX[0] -= Momloss*MeasPX[0]/TMath::Sqrt(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0]); 
      MeasPY[0] -= Momloss*MeasPY[0]/TMath::Sqrt(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0]); 
      MeasPZ[0] -= Momloss*MeasPZ[0]/TMath::Sqrt(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0]); 

 cout << Momloss << endl;
 cout << "Momentum " << TMath::Sqrt(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0])/MeVcToSi << endl;

      gr_P->SetPoint(0,MeasX[0],MeasY[0],MeasZ[0]);
      gr_X->SetPoint(0,0,MeasX[0]);
      gr_Y->SetPoint(0,0,MeasY[0]);
      gr_Z->SetPoint(0,0,MeasZ[0]);

   for(int hit=1; hit<NHit; hit++){

      Vinst = 3.0e8*TMath::Sqrt(1.-1./(1+(MeasPX[hit-1]*MeasPX[hit-1]+MeasPY[hit-1]*MeasPY[hit-1]+MeasPZ[hit-1]*MeasPZ[hit-1])/(mymass*mymass*3.0e8*3.0e8)));
      Dt = dt * TMath::Sqrt(1.-Vinst*Vinst /(3.0e8*3.0e8) ); // because particles are relativistic

      MeasX[hit] = MeasX[hit-1] + MeasPX[hit-1]*Dt/mymass+MeasPY[hit-1]*mycharge*B*Dt*Dt/(2.*mymass*mymass); // X position of the object, no error on the process
         MeasXr[hit] = Vrand->Gaus(MeasX[hit],sigX);
      MeasY[hit] = MeasY[hit-1] + MeasPY[hit-1]*Dt/mymass-MeasPX[hit-1]*mycharge*B*Dt*Dt/(2.*mymass*mymass); // X position of the object, no error on the process
         MeasYr[hit] = Vrand->Gaus(MeasY[hit],sigY);
      MeasZ[hit] = MeasZ[hit-1] + MeasPZ[hit-1]*Dt/mymass;             				             // Z position of the object, no error on the process
         MeasZr[hit] = Vrand->Gaus(MeasZ[hit],sigZ);

	 //      h_Rgrowth->Fill(TMath::Sqrt(MeasX[hit]*MeasX[hit]+MeasY[hit]*MeasY[hit])-TMath::Sqrt(MeasX[hit-1]*MeasX[hit-1]+MeasY[hit-1]*MeasY[hit-1]));

      MeasPX[hit] = MeasPX[hit-1] + MeasPY[hit-1]*mycharge*B*Dt/mymass; // momentum along X
      MeasPY[hit] = MeasPY[hit-1] - MeasPX[hit-1]*mycharge*B*Dt/mymass; // momentum along Y
      MeasPZ[hit] = MeasPZ[hit-1];                                      // momentum along Z
         Momloss = -1.;
         while(Momloss<0){
            Momloss = Vrand->Gaus(0.,sigMomloss);
         }
      MeasPX[hit] -= Momloss*MeasPX[hit]/TMath::Sqrt(MeasPX[hit]*MeasPX[hit]+MeasPY[hit]*MeasPY[hit]+MeasPZ[hit]*MeasPZ[hit]); // Initial X momentum of the object
      MeasPY[hit] -= Momloss*MeasPY[hit]/TMath::Sqrt(MeasPX[hit]*MeasPX[hit]+MeasPY[hit]*MeasPY[hit]+MeasPZ[hit]*MeasPZ[hit]); // Initial Y momentum of the object
      MeasPZ[hit] -= Momloss*MeasPZ[hit]/TMath::Sqrt(MeasPX[hit]*MeasPX[hit]+MeasPY[hit]*MeasPY[hit]+MeasPZ[hit]*MeasPZ[hit]); // Initial Z momentum of the object

      if(hit%3000==0){ // the helix fit can fit 200 points maximum
         pData[ind_pData][0]=MeasXr[hit];
         pData[ind_pData][1]=MeasYr[hit];
         pData[ind_pData][2]=MeasZr[hit];

         gr_P->SetPoint(ind_pData,MeasXr[hit],MeasYr[hit],MeasZr[hit]);
         gr_X->SetPoint(ind_pData,hit,MeasXr[hit]);
         gr_Y->SetPoint(ind_pData,hit,MeasYr[hit]);
         gr_Z->SetPoint(ind_pData,hit,MeasZr[hit]);
         gr_XY->SetPoint(ind_pData,MeasXr[hit],MeasYr[hit]);

         gr_PX->SetPoint(ind_pData,hit,MeasPX[hit]);
         gr_PY->SetPoint(ind_pData,hit,MeasPY[hit]);
         gr_PZ->SetPoint(ind_pData,hit,MeasPZ[hit]);
         gr_PT->SetPoint(ind_pData,hit,TMath::Sqrt(MeasPX[hit]*MeasPX[hit]+MeasPY[hit]*MeasPY[hit]));

         ind_pData++;
      }

      DistTot += TMath::Sqrt((MeasX[hit]-MeasX[hit-1])*(MeasX[hit]-MeasX[hit-1])+(MeasY[hit]-MeasY[hit-1])*(MeasY[hit]-MeasY[hit-1])+(MeasZ[hit]-MeasZ[hit-1])*(MeasZ[hit]-MeasZ[hit-1]));


   } // hit   


// As the Kalman filter needs to be initialized we use the global Helix Fit for this. It will also be used afterwards for comparison

           double tmpR=0.0; double tmpA=0.; double tmpB=0.; double tmpPhi=0.; double tmpTheta=0.; double tmpZ=0.;
           HelixFit(ind_pData,pData,tmpR,tmpA,tmpB,tmpPhi,tmpTheta,tmpZ,1);
	   // K.Park : Call helix_fit from BonusHelixFit.h
           //helix_fit(ind_pData,pData,tmpR,tmpA,tmpB,tmpPhi,tmpTheta,tmpZ,1);

	   
   double pT_gfit = (mycharge/elemcharge) * 300 * B * tmpR * MeVcToSi;
   /*
   cout << "********************** Using the global helix fit ****************************" << endl;
   cout << "Transverse momentum resolution: " << 100.*((mycharge/elemcharge) * 300 * B * tmpR * MeVcToSi - pT_gen)/pT_gen << " %." << endl;
   cout << "Longitudinal momentum resolution: " << 100.*(pT_gfit/TMath::Tan(tmpTheta*TMath::Pi()/180.) - PZ_gen)/PZ_gen << " %." << endl;
   cout << "X momentum resolution: " << 100.*(pT_gfit*TMath::Cos(tmpPhi*TMath::Pi()/180.) - PX_gen)/PX_gen << " %." << endl;
   cout << "Y momentum resolution: " << 100.*(pT_gfit*TMath::Sin(tmpPhi*TMath::Pi()/180.) - PY_gen)/PY_gen << " %." << endl;
   cout << "Delta Phi: " << (tmpPhi*TMath::Pi()/180. - phi_gen)*1000. << " mrad." << endl;
   cout << "Delta Theta: " << (tmpTheta*TMath::Pi()/180. - theta_gen)*1000. << " mrad." << endl;
   cout << "Delta Z: " << tmpZ - Z_gen << " mm." << endl;
   cout << "******************************************************************************" << endl;
   cout << "  " << endl;
   //cout << "_________________" <<  pT_gfit << "________________"<< pT_gen <<endl;
   */

   pTresol_helix = 100.*((mycharge/elemcharge) * 300 * B * tmpR * MeVcToSi - pT_gen)/pT_gen;
   pXresol_helix = 100.*(pT_gfit*TMath::Cos(tmpPhi*TMath::Pi()/180.) - PX_gen)/PX_gen;
   pYresol_helix = 100.*(pT_gfit*TMath::Sin(tmpPhi*TMath::Pi()/180.) - PY_gen)/PY_gen;

   rDiff_helix= TMath::Sqrt(tmpA*tmpA+tmpB*tmpB) - TMath::Sqrt(X_gen*X_gen + Y_gen*Y_gen) ;
   thetaDiff_helix= (tmpTheta*TMath::Pi()/180. - theta_gen)*1000.; // unit is mrad
   phiDiff_helix =(tmpPhi*TMath::Pi()/180. - phi_gen)*1000. ; // unit is mrad
   zDiff_helix = tmpZ - Z_gen; // unit is mm

   //   cout << "tmpR =" << tmpR << ", tmpXY =" << TMath::Sqrt(tmpA*tmpA+tmpB*tmpB) << endl;
   
   //Percentage of transverse, x, y momentum resolutiuon  from Global Helix fit
   OUT1 << setiosflags(ios::fixed) << pTresol_helix  <<", " << pXresol_helix   <<", " << pYresol_helix << ", " << rDiff_helix << ", " << thetaDiff_helix << ", " << phiDiff_helix << ", " << zDiff_helix << endl;


   {
     ntuple1->Fill(pTresol_helix,pXresol_helix,pYresol_helix,rDiff_helix,thetaDiff_helix,phiDiff_helix,zDiff_helix);
   }

   //--Kalman Filter
      // Initialize values (why KF initial only for tmpZ ??)
      //Xkm1(0,0) = 0.1;
      //Xkm1(1,0) = 0.1;

      Xkm1(0,0) = tmpA;
      Xkm1(1,0) = tmpB;      
      Xkm1(2,0) = tmpZ;
      Xkm1(3,0) = -0.5*pT_gfit/(MeasX[NHit]+1e2*tmpR*TMath::Cos(tmpPhi*TMath::Pi()/180.));
      Xkm1(4,0) = -0.5*pT_gfit/(MeasY[NHit]-1e2*tmpR*TMath::Sin(tmpPhi*TMath::Pi()/180.));//-pT_gfit*TMath::Sin(tmpPhi*TMath::Pi()/180.);
      Xkm1(5,0) = pT_gfit/TMath::Tan(tmpTheta*TMath::Pi()/180.);

      double PosErr = 3.;
      double PosErr2 = PosErr*PosErr; 
      double MomErr = 5.0;
      double MomErr2 = MomErr*MomErr;
      double CrossEr = 3.e-6;
      Pkm1(0,0) = PosErr2; Pkm1(0,1) = 0.;      Pkm1(0,2) = 0.;      Pkm1(0,3) = 0.;      Pkm1(0,4) = 0.;      Pkm1(0,5) = 3.;
      Pkm1(1,0) = 0.;  	   Pkm1(1,1) = PosErr2; Pkm1(1,2) = 0.;      Pkm1(1,3) = 0.;      Pkm1(1,4) = 0.;      Pkm1(1,5) = 3.;
      Pkm1(2,0) = 0.;  	   Pkm1(2,1) = 0.;      Pkm1(2,2) = PosErr2; Pkm1(2,3) = 0.;      Pkm1(2,4) = 0.;      Pkm1(2,5) = 0.;
      Pkm1(3,0) = 0.;  	   Pkm1(3,1) = 0.;      Pkm1(3,2) = 0.;      Pkm1(3,3) = MomErr2; Pkm1(3,4) = 0.;      Pkm1(3,5) = 0.;
      Pkm1(4,0) = 0.;      Pkm1(4,1) = 0.;      Pkm1(4,2) = 0.;      Pkm1(4,3) = 0.; 	  Pkm1(4,4) = MomErr2; Pkm1(4,5) = 0.;
      Pkm1(5,0) = 0.;  	   Pkm1(5,1) = 0.;      Pkm1(5,2) = 0.;      Pkm1(5,3) = 3.;      Pkm1(5,4) = 3.;      Pkm1(5,5) = MomErr2;
   int hit = NHit; 
   while(hit>0){

      // First, compute the values for the matrix A at this step
      Vinst = 3.0e8*TMath::Sqrt(1.-1./(1+(Xkm1(3,0)*Xkm1(3,0)+Xkm1(4,0)*Xkm1(4,0)+Xkm1(5,0)*Xkm1(5,0))/(mymass*mymass*3.0e8*3.0e8)));
      Dt = dt * TMath::Sqrt(1.-Vinst*Vinst /(3.0e8*3.0e8)); // because particles are relativistic

      alpha = Dt/mymass;
      beta = mycharge*B*Dt*Dt/(2.*mymass*mymass);
      gamma = mycharge*B*Dt/mymass;

      A(0,3) = -(alpha+beta*gamma)/(1.+gamma*gamma);
      A(0,4) = (alpha*gamma-beta)/(1.+gamma*gamma);
      A(1,3) = (beta-alpha*gamma)/(1.+gamma*gamma);
      A(1,4) = -(alpha+beta*gamma)/(1.+gamma*gamma);
      A(2,5) = -alpha;  
      A(3,3) = 1./(1.+gamma*gamma);
      A(3,4) = -gamma/(1.+gamma*gamma);
      A(4,3) = gamma/(1.+gamma*gamma);
      A(4,4) = 1./(1.+gamma*gamma);

      // Prediction
        Xkm = A*Xkm1; 
        Pkm = A*Pkm1*AT + Q;
//if(hit>NHit-50) cout << hit << "  " << TMath::Sqrt(MeasX[hit]*MeasX[hit]+MeasY[hit]*MeasY[hit]) << "  " << TMath::Sqrt(Xkm(0,0)*Xkm(0,0)+Xkm(1,0)*Xkm(1,0)) << "  " << dR << endl;
      // Update only if next hit has been reconstructed at this radius 
        if(hit%6000==0/*hit<NHit && TMath::Abs(TMath::Sqrt(MeasX[hit]*MeasX[hit]+MeasY[hit]*MeasY[hit])-TMath::Sqrt(Xkm(0,0)*Xkm(0,0)+Xkm(1,0)*Xkm(1,0)))<dR*/){
           Kk = Pkm*HT * (H*Pkm*HT + R).Invert();
           Temp(0,0) = MeasX[hit]; // necessary as operation works only for matrixes
           Temp(1,0) = MeasY[hit];
           Temp(2,0) = MeasZ[hit];
           Xk = Xkm + Kk*(Temp-H*Xkm);
           Pk = (Id-Kk*H) * Pkm;
           hit--;

           gr_P_Kal->SetPoint(ind_kal,Xk(0,0),Xk(1,0),Xk(2,0));
           gr_X_Kal->SetPoint(ind_kal,hit,Xk(0,0));
           gr_Y_Kal->SetPoint(ind_kal,hit,Xk(1,0));
           gr_Z_Kal->SetPoint(ind_kal,hit,Xk(2,0));
           gr_XY_Kal->SetPoint(ind_kal,Xk(0,0),Xk(1,0));

           gr_PX_Kal->SetPoint(ind_kal,hit,Xk(3,0));
           gr_PY_Kal->SetPoint(ind_kal,hit,Xk(4,0));
           gr_PZ_Kal->SetPoint(ind_kal,hit,Xk(5,0));
           gr_PT_Kal->SetPoint(ind_kal,hit,TMath::Sqrt(Xk(3,0)*Xk(3,0)+Xk(4,0)*Xk(4,0)));

           ind_kal++;
      // Save the Xk value and transfer the k values to k-1 
           Xkm1=Xk;
           Pkm1=Pk;
        }
        else if(hit==NHit){
           Kk = Pkm*HT * (H*Pkm*HT + R).Invert();
           Temp(0,0) = MeasX[hit]; // necessary as operation works only for matrixes
           Temp(1,0) = MeasY[hit];
           Temp(2,0) = MeasZ[hit];
           Xk = Xkm + Kk*(Temp-H*Xkm);
           Pk = (Id-Kk*H) * Pkm;
           hit--;

           gr_P_Kal->SetPoint(ind_kal,Xk(0,0),Xk(1,0),Xk(2,0));
           gr_X_Kal->SetPoint(ind_kal,hit,Xk(0,0));
           gr_Y_Kal->SetPoint(ind_kal,hit,Xk(1,0));
           gr_Z_Kal->SetPoint(ind_kal,hit,Xk(2,0));
           gr_XY_Kal->SetPoint(ind_kal,Xk(0,0),Xk(1,0));

           gr_PX_Kal->SetPoint(ind_kal,hit,Xk(3,0));
           gr_PY_Kal->SetPoint(ind_kal,hit,Xk(4,0));
           gr_PZ_Kal->SetPoint(ind_kal,hit,Xk(5,0));
           gr_PT_Kal->SetPoint(ind_kal,hit,TMath::Sqrt(Xk(3,0)*Xk(3,0)+Xk(4,0)*Xk(4,0)));

           ind_kal++;
      // Save the Xk value and transfer the k values to k-1 
           Xkm1=Xk;
           Pkm1=Pk;
        }
        else{ // no udpate of the matrixes
           Xkm1=Xkm;
           Pkm1=Pkm;
           hit--;
        }
   } 
   //--

   phi_kal = TMath::ATan2(Xk(4,0),Xk(3,0));
      if(phi_kal<0) phi_kal+=2.*TMath::Pi();
   theta_kal = TMath::ATan(TMath::Sqrt(Xk(3,0)*Xk(3,0)+Xk(4,0)*Xk(4,0))/Xk(5,0));
      if(theta_kal<0) theta_kal+=TMath::Pi();
      /*
   cout << "************************ Using the Kalman filter *****************************" << endl;
   cout << ind_kal << " points by the Kalman filter." << endl;
   cout << "Transverse momentum resolution: " << 100.*(TMath::Sqrt(Xk(3,0)*Xk(3,0)+Xk(4,0)*Xk(4,0)) - pT_gen)/pT_gen << " %." << endl;
   cout << "Longitudinal momentum resolution: " << 100.*(Xk(5,0) - PZ_gen)/PZ_gen << " %." << endl;
   cout << "X momentum resolution: " << 100.*(Xk(3,0) - PX_gen)/PX_gen << " %." << endl;
   cout << "Y momentum resolution: " << 100.*(Xk(4,0) - PY_gen)/PY_gen << " %." << endl;
   cout << "Delta Phi: " << (phi_kal - phi_gen)*1000. << " mrad." << endl;
   cout << "Delta Theta: " << (theta_kal - theta_gen)*1000. << " mrad." << endl;
   cout << "Delta Z: " << Xk(2,0) - Z_gen << " mm." << endl;
   cout << "******************************************************************************" << endl;
      */

   pTresol_kf = 100.*(TMath::Sqrt(Xk(3,0)*Xk(3,0)+Xk(4,0)*Xk(4,0)) - pT_gen)/pT_gen;
   pXresol_kf = 100.*(Xk(3,0) - PX_gen)/PX_gen;
   pYresol_kf = 100.*(Xk(4,0) - PY_gen)/PY_gen;

   
   rDiff_kf = TMath::Sqrt(Xk(0,0)*Xk(0,0)+Xk(1,0)*Xk(1,0)) - TMath::Sqrt(X_gen*X_gen + Y_gen*Y_gen) ;
   thetaDiff_kf = (theta_kal - theta_gen)*1000. ; // unit is mrad
   phiDiff_kf = (phi_kal - phi_gen)*1000. ; // unit is mrad
   zDiff_kf = Xk(2,0) - Z_gen; // unit is mm

   //Percentage of transverse, x, y momentum resolutiuon  from Global Helix fit
   OUT2 << setiosflags(ios::fixed) << "kf," << pTresol_kf <<", " <<  pXresol_kf <<", " <<  pYresol_kf << rDiff_kf << ", " << thetaDiff_kf << ", " << phiDiff_kf << ", " << zDiff_kf << endl;


   {
     ntuple2->Fill(pTresol_kf,pXresol_kf,pYresol_kf,rDiff_kf,thetaDiff_kf,phiDiff_kf,zDiff_kf);
   }

   
//________________________________________________________________________________________________
// _________________________________________ Displays ____________________________________________
//________________________________________________________________________________________________

cout << "Initial particle speed: " << 3.0e8*TMath::Sqrt(1.-1./(1+(MeasPX[0]*MeasPX[0]+MeasPY[0]*MeasPY[0]+MeasPZ[0]*MeasPZ[0])/(mymass*mymass*3.0e8*3.0e8))) << " m/s." << endl;
cout << "Final particle speed: " << 3.0e8*TMath::Sqrt(1.-1./(1+(MeasPX[NHit-1]*MeasPX[NHit-1]+MeasPY[NHit-1]*MeasPY[NHit-1]+MeasPZ[NHit-1]*MeasPZ[NHit-1])/(mymass*mymass*3.0e8*3.0e8))) << " m/s." << endl;
cout << "Particle path length: " << DistTot*1000 << " mm in " << (double) (NHit*dt) << " s, i.e. a speed of " << DistTot/ (double) (NHit*dt) << " m/s." << endl;


/*  no plotting....
 
   c0->cd();
   gr_P->Draw("surf1");

   cX->cd();
      mg_X->Add(gr_X);
         gr_X_Kal->SetLineColor(kRed);
         gr_X_Kal->SetMarkerColor(kRed);
      mg_X->Add(gr_X_Kal);
      mg_X->Draw("apl");

   cY->cd();
      mg_Y->Add(gr_Y);
         gr_Y_Kal->SetLineColor(kRed);
         gr_Y_Kal->SetMarkerColor(kRed);
      mg_Y->Add(gr_Y_Kal);
      mg_Y->Draw("apl");

   cZ->cd();
      mg_Z->Add(gr_Z);
         gr_Z_Kal->SetLineColor(kRed);
         gr_Z_Kal->SetMarkerColor(kRed);
      mg_Z->Add(gr_Z_Kal);
      mg_Z->Draw("apl");

   cXY->cd();
      mg_XY->Add(gr_XY);
         gr_XY_Kal->SetLineColor(kRed);
         gr_XY_Kal->SetMarkerColor(kRed);
      mg_XY->Add(gr_XY_Kal);
      mg_XY->Draw("apl");


   cPX->cd();
      mg_PX->Add(gr_PX);
         gr_PX_Kal->SetLineColor(kRed);
         gr_PX_Kal->SetMarkerColor(kRed);
      mg_PX->Add(gr_PX_Kal);
      mg_PX->Draw("apl");

   cPY->cd();
      mg_PY->Add(gr_PY);
         gr_PY_Kal->SetLineColor(kRed);
         gr_PY_Kal->SetMarkerColor(kRed);
      mg_PY->Add(gr_PY_Kal);
      mg_PY->Draw("apl");

   cPZ->cd();
      mg_PZ->Add(gr_PZ);
         gr_PZ_Kal->SetLineColor(kRed);
         gr_PZ_Kal->SetMarkerColor(kRed);
      mg_PZ->Add(gr_PZ_Kal);
      mg_PZ->Draw("apl");

   cPT->cd();
      mg_PT->Add(gr_PT);
         gr_PT_Kal->SetLineColor(kRed);
         gr_PT_Kal->SetMarkerColor(kRed);
      mg_PT->Add(gr_PT_Kal);
      mg_PT->Draw("apl");

   cdR->cd();
   h_Rgrowth->Draw();

no plottinng   */

 
}
   OUT1.close();
   OUT2.close();

   f1->Write();
   f2->Write();
   
}

