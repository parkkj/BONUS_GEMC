
#include <stdio.h>
#include <math.h>
#include "BonusHelixFit.h"
#include "rwfthc.cc"
//#define HELIXFIT_DEBUG 1
/*
      SUBROUTINE RWFTHL(NPT,RF,PF,WFI,ZF,WZF,IOPT,
     1                  VV0,EE0,CH2PH,CH2Z,DEL,DELZ)
      implicit none
C
C-----------------------------------------------------------------------
C! Fast helix fit
C
C   A generalization of the TFTHEL routine to allow it to be called
C   from a routine that contains any list of x and y values XF,YF for a
C   set of NPT points to be fitted.
C
C   Input:  NPT    /I     Number of 3-D points to be fit
C           XF     /R     Array of X-values of points to be fit
C           YF     /R     Array of Y-values of points to be fit
C           RF     /R     Array of R-values of points to be fit
C           PF     /R     Array of PHI-values of points to be fit
C           WFI    /R     Array of 1/(sig(rphi))**2 for each point
C           ZF     /R     Array of Z-values of points to be fit
C           WZF    /R     Array of 1/(sig(z))**2 for each point
C           IOPT = 0 -> DISTANCE**2=X**2+Y**2 MINIMISED
C                  1 -> WEIGHTED WITH 1/SIMA(R*PHI)**2
C                  2 -> ERROR MATRIX CALCULATED
C                  3 -> 3-DIMENSIONAL ITERATION
C  OUTPUT:   VV0 = 1/R*CHARGE    POS. IF CLOCKWISE
C                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE)
C                  PHI0         {0,2TMath::Pi()} ANGLE TO X-AXIS at R=D0
C                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,
C                                       +VE IF AXIS ENCIRCLED
C                  Z0           [CM]    Z POS AT R=D0
C          EE0 = INVERSE OF ERROR MATRIX IN TRIANG. FORM
C          CH2PH = CHI SQUARED = SUM (PHI DEVIATIONS/ERRORS)**2
C          CH2Z  = CHI SQUARED = SUM (Z DEVIATIONS/ERRORS)**2
C          DEL = ??
c          DELZ= ??
C  NOTE: DEGREES OF FREEDOM = 2*NPT-5
C----------------------------------------------------------------
C     BASED ON  SUBROUTINE CIRCLE
C     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329
C
C   AUTHORS:  N. CHERNOV, G. OSOSKOV & M. POPPE
C   Modified by:  Fred Weber, 8 Jun 1989
C
C-----------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C" {
#endif 
//file  rwfthl.f
void rwfthl_(int* npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
	     int* iopt, float* vv0, float* ee0, 
	     float* ch2ph, float* ch2z,  /* return values */ 
	     float* del, float* delz);
 
//file  rwfthc.cc
void rwfthc(int npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
	     int iopt, float* vv0, float* ee0, 
	     float* ch2ph, float* ch2z,  /* return values */ 
	     float* del, float* delz);

void rwfthl_c(int npt, float* rf, float* pf, float* wfi, float* zf, float* wzf, 
	      int iopt, float* vv0, float* ee0, 
	      float* ch2ph, float* ch2z,  /* return values */ 
	      float* del, float* delz)
{ 
    rwfthc(  npt,rf,pf,wfi,zf,wzf, iopt,vv0,ee0, ch2ph, ch2z, del,delz);
}


//-----------------------------------------------------------------

#ifdef __cplusplus
}
#endif 

/*------------------------------------------------------------------------\
Function name: void helix_fit(int PointNum,double szPos[][3],
                              double& R, double& X, double& Y));
Calculate the raidus of the trajectory
Input parameters:
	PointNum:       number of x-y points
	szPos[any number>3][3]:  xyz array, in unit of mm
OutPut: R    ->Radius in mm, 
       (x,y) ->origin position (x,y) in mm

\------------------------------------------------------------------------*/


void helix_fit(int PointNum,double szPos[][3], double& R, double& X, double& Y,
	       double& Phi_deg, double& Theta_deg, double& Z0,int fit_track_to_beamline )
{
  const int MAX_HITS = 255;
  
  int jj;
  float my_phi;
  //float xf[MAX_HITS];
  //float yf[MAX_HITS];
  float rf[MAX_HITS];
  float pf[MAX_HITS];
  float wfi[MAX_HITS];
  float zf[MAX_HITS];
  float wzf[MAX_HITS];
  int iopt;
  int npt;
  float vv0[5];
  float ee0[15];
  float ch2ph;
  float ch2z;
  float del[MAX_HITS];
  float delz[MAX_HITS];

  float phi0;
  npt=0;
  if(PointNum>MAX_HITS) PointNum=MAX_HITS-1;
  for (jj=0; jj<PointNum; jj++)
    {// r,phi,z coordinate
#ifdef HELIXFIT_DEBUG
	  printf("point %3d: X=%.2f  Y=%.2f  Z=%.2f\n",jj+1,szPos[jj][0],szPos[jj][1],szPos[jj][2]);
#endif
      rf[jj] = sqrt(pow(szPos[jj][0],2)+pow(szPos[jj][1],2));
      pf[jj] = atan(szPos[jj][1]/szPos[jj][0]); //phi angle
      if(szPos[jj][1]>0 && szPos[jj][0]<0) pf[jj] +=TMath::Pi();
      if(szPos[jj][1]<0 && szPos[jj][0]<0) pf[jj] +=TMath::Pi();
      if(szPos[jj][1]<0 && szPos[jj][0]>0) pf[jj] +=2*TMath::Pi();
      if(pf[jj]>2*TMath::Pi())    pf[jj] -=2*TMath::Pi();
      zf[jj] = szPos[jj][2];
      wfi[jj]= 1.0;
      wzf[jj]= 1.0;
    }
  npt=PointNum;

  if(fit_track_to_beamline)
    {
      rf[npt]= 0.0;
      pf[npt]= 0.0;
      zf[npt]= 0.0; 
      //zf[npt]=Z0; modified by jixie, should not do that because in real data we don't know z0
      wfi[npt]= 1.0;
      wzf[npt]= 0.0; // zero weight for Z on the beamline point
	  //This means that don't calculate the chi square for Z on the beamline point
      npt++;
    }
  iopt= 1; /* tells rwfthl what kind of fit to do */

  rwfthl_c(npt,rf,pf,wfi,zf,wzf,iopt,vv0,ee0, &ch2ph, &ch2z, del,delz);
/*
 OUTPUT:   VV0 = 1/R*CHARGE    POS. IF CLOCKWISE
C                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE)
C                  PHI0         {0,2TMath::Pi()} ANGLE TO X-AXIS at R=D0
C                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,
C                                       +VE IF AXIS ENCIRCLED
C                  Z0           [CM]    Z POS AT R=D0
*/
  //reconstruct the output
  R  = (double)fabs(1.0/vv0[0]); /* minimum distance to z=0 */
  phi0 = vv0[2]; /* in xy plane, direction of track relative to x axis */

  //This is from Nate, I don't understand why they use it like this
  my_phi = phi0+TMath::Pi();
  if (vv0[0]<0.0) my_phi+=TMath::Pi();
  //vv0[0] negtive means curve to anti-CLOCKWISE direction, This is true in BONUS
  if(my_phi>=2.0*TMath::Pi()) my_phi-=2.0*TMath::Pi();

  //center of the circle
  X = (double)(-sin(my_phi)*((-vv0[3])+fabs(1.0/vv0[0])));
  Y = (double)(cos(my_phi)*((-vv0[3])+fabs(1.0/vv0[0])));
  //position of the initial step
  Phi_deg=180.0/3.14159*phi0;
  if(Phi_deg<=-360.0) Phi_deg+=360.;
  else if(Phi_deg>=360.0) Phi_deg-=360.;
  Theta_deg=180.0/3.14159*(3.14159/2.-atan(vv0[1]));
  Z0 =  vv0[4];
///////////////////////////////////I don't need this part
// float dzds,dca,x_close,y_close,z_close,m_s,b_s,chi2;
// dzds = vv0[1];
// dca = fabs(vv0[3]); /* dca = distance of closest approach to beamline */
// x_close = -sin(my_phi)*(-vv0[3]);
// y_close =  cos(my_phi)*(-vv0[3]);
// z_close =  vv0[4]; /* z at r=r_0 */
// m_s = 1/tan(dzds); /* pass it this way for now since other routines expect it */
// b_s = z_close;
// if(npt>5) {chi2 = ch2ph+ch2z/(npt-5);} else {chi2 = 9999.9;}
///////////////////////////////////

#ifdef HELIXFIT_DEBUG
    printf(" +++++ mohammad: BonusHelixFit.cc: fitting %d hits then return (a,b,r)= (%6.4f %6.4f %6.4f)\n ++++++", PointNum,X,Y,R);
#endif
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void HelixFit(int PointNum,double szPos[][3], double& R, double& X, double& Y,
	      double& Phi_deg, double& Theta_deg, double& Z0,int fit_track_to_beamline )
{

   helix_fit(PointNum,szPos, R, X, Y, Phi_deg, Theta_deg, Z0, fit_track_to_beamline);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void HelixFit(int PointNum,double szPos[][3], double& R, double& X, double& Y,
	      float& Phi_deg, float& Theta_deg, float& Z0,int fit_track_to_beamline )
{
  double PPhi_deg,TTheta_deg,ZZ0;
  helix_fit(PointNum,szPos,R,X,Y,PPhi_deg,TTheta_deg,ZZ0,fit_track_to_beamline);
  Phi_deg=(float) PPhi_deg;
  Theta_deg=(float) TTheta_deg;
  Z0=(float) ZZ0;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void HelixFit(int PointNum,double szPos[][3], float& R, float& X, float& Y,
	      float& Phi_deg, float& Theta_deg, float& Z0,int fit_track_to_beamline )
{
  double RR,XX,YY,PPhi_deg,TTheta_deg,ZZ0;
  helix_fit(PointNum,szPos, RR, XX, YY,PPhi_deg, TTheta_deg, ZZ0, fit_track_to_beamline);
  R=(float) RR;
  X=(float) XX;
  Y=(float) YY;
  Phi_deg=(float) PPhi_deg;
  Theta_deg=(float) TTheta_deg;
  Z0=(float) ZZ0;
}






