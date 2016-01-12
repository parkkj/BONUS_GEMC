//Header file for helix fit, c++ code!!!
#ifndef HELIX_FIT
#define HELIX_FIT

#include "BonusHelixFit.c"
///////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------\
Function name: void helix_fit(int PointNum,double szPos[][3],
                              double& R, double& X, double& Y));
Calculate the radius of the trajectory
Input parameters:
	PointNum:       number of x-y points
	szPos[any number>3][3]:  xyz array, in unit of mm
OutPut: R    ->Radius in mm, 
       (x,y) ->origin position (x,y) in mm
       Phi_deg -> The phi angle at the initial step
       Theta_deg -> The theta angle at the initial step
       Z0   ->The vertex position at the initial step
\------------------------------------------------------------------------*/
  void helix_fit(int PointNum,double szPos[][3], double& R, double& X, double& Y, 
		 double& Phi_deg, double& Theta_deg, double& Z0,int fit_track_to_beamline=1);

  void HelixFit(int PointNum,double szPos[][3], double& R, double& X, double& Y, 
		double& Phi_deg, double& Theta_deg, double& Z0,int fit_track_to_beamline=1);
  
  void HelixFit(int PointNum,double szPos[][3], double& R, double& X, double& Y,
		float& Phi_deg, float& Theta_deg, float& Z0,int fit_track_to_beamline=1);
  void HelixFit(int PointNum,double szPos[][3], float& R, float& X, float& Y,
		float& Phi_deg, float& Theta_deg, float& Z0,int fit_track_to_beamline=1);
  

////////////////////////////////////////////////////////////////////////
#endif //#ifdef _HELIX_FIT_
