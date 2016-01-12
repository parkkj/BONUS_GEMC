// gemc headers
// Implementing Electron Transport Formulat from MagBoltz @ 11/19/2015, K.Park
//
#include "bonus_hitprocess.h"
#include "DriftEMagboltz.hh"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

float PAD_W =2.5;
float PAD_L = 4.0;
float PAD_S = 80.0;
float RTPC_L = 400.0;

// GetFunction for drift time and dphi
//get drift time by s
float GetT_s2gem1(float s0_mm,float z0){
  z0+=0;
  float s0=s0_mm/10.;
  float t_us=0.0; 
  float a=-0.1642, b=-0.0947,c=8.8001;
  t_us = a*s0*s0+b*s0+c;
  return t_us*1000.;
}

//get dphi by s
float GetdPhi_s2gem1(float s0_mm,float z0){
  z0+=0;
  float s0=s0_mm/10.;
  float a=0.0287, b=-0.5334, c=2.3475;
  float dphi=a*s0*s0+b*s0+c;
  return dphi;
}

int GetChanId(float z0, float phi_rad)
{
  //need to fill this routine later
  //row shifting in z : 
  //row 0|1|2|3: 0|1|2|3 mm
  float phi_per_pad = PAD_W/PAD_S;
  int row = int(phi_rad/phi_per_pad);
  float z_shift = row%4;
  int col = int((z0-z_shift+RTPC_L/2)/PAD_L);
  int Num_of_Col = int(ceil(RTPC_L/PAD_L));
  return row*Num_of_Col+col;
}

float p_pr=0.1, m_pr=0.9383;
float c=300.0;  //speed of light 3E8 m/s  or 300 mm/ns
float v=p_pr/sqrt(m_pr*m_pr+p_pr*p_pr)*c;  //proton travel speed
//float PAD_S = 80.0;
//int NS_PER_TIC = 100;
float TPC_TZERO = int(PAD_S/v/NS_PER_TIC);
float Kev2ADC =20.0;	 //set to 20.0, so 1 adc means 50 ev


map<string, double> bonus_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
	vector<identifier> identity = aHit->GetId();
	
	// true information
	// for example tInfos.eTot is total energy deposited
	trueInfos tInfos(aHit);
	
	// local variable for each step
	vector<G4ThreeVector> Lpos = aHit->GetLPos();
	
	// energy at each step
	// so tInfos.eTot is the sum of all steps s of Edep[s] 
	vector<double>      Edep = aHit->GetEdep();

	//
	// Get the information x,y,z and Edep at each ionization point
	// 

	double LposX=0.;
	double LposY=0.;
	double LposZ=0.;
	double DiffEdep=0.;
	if(tInfos.eTot > 0)
	  {
	  for(unsigned int s=0; s<tInfos.nsteps; s++)
	    {
	      LposX = Lpos[s].x();
	      LposY = Lpos[s].y();
	      LposZ = Lpos[s].z();
	      DiffEdep = Edep[s];
	      // From this section for electron transport formular
	      // reference from Jixie's G4 original C++ code

	      double r0,phi0_rad;
	      //convert (x0,y0,z0) into (r0,phi0,z0)
	      r0=sqrt(LposX*LposX+LposY*LposY);  //in mm
	      phi0_rad=atan2(LposY,LposX); //return (-Pi, + Pi)
	      if( phi0_rad<0.)  phi0_rad+=2.0*PI;
	      
	      //determine drift_time
	      float t_off;	//time offset
	      float t_s2gem1;   //the driftime(ns) from r0(x0,y0,z0) to the first gem
	      float t_s2pad;    //the driftime(ns) from r0(x0,y0,z0) to the pad
	      
	      //1. calculate toff.
	      //	t_off = GetT_gem2pad(tInfos.lz0); // actually it returns ZERO !
	      t_off = 0;
	      //2.determine drift time
	      t_s2gem1 = GetT_s2gem1(r0,LposZ);
	      t_s2pad = t_s2gem1 + t_off;
	      
	      //GetSByT(float t_s2gem1,float z0)
	      float s_r=-10.0;
	      float t_us=t_s2gem1/1000.;
	      float a=-0.0335, b=-0.3208,c=6.9481; // Gail's MagBoltz number
	      s_r = (a*t_us*t_us+b*t_us+c)*10.;   //in mm
	      	      
	      float dphi_offset;    //all phis are in rad unit
	      float dphi_s2gem1;    //delta phi (phi change) from r0(x0,y0,z0) to the 1st gem
	      float delta_phi;      //delta phi (phi change) from r0(x0,y0,z0) to the pad
	      
	      //these 2 steps add up together is eqal to getdelphi(###)
	      //1. calculate phioff
	      //dphi_offset=GetdPhi_gem2pad(tInfos.lz0); // actually it returns ZERO !
	      dphi_offset= 0;
	      //3.calculate delta phi (phi change) from r0(x0,y0,z0) to the 1st gem
	      dphi_s2gem1=GetdPhi_s2gem1(r0,LposZ);
	      
	      delta_phi = dphi_s2gem1 + dphi_offset;
	      
	      float phi_rad=float(phi0_rad-delta_phi);   //phi at pad pcb board
	      //this is not the center of the pad yet!!!
	      //convert the phi_deg into range [0,2*PI)
	      if( phi_rad<0. )  phi_rad+=2.0*PI;
	      
	      //output////////////////////////////////////////////////////////////////
	      class cellID
	      {
		int phiCellID;
		int zCellID;
		double energyRatio=1.;
		// fraction of total energy assigned to this cell  - e.g. 1 = 100% of total energy
	      };

	      vector<cellID> findPads(G4Threevector xyz)
	      {
		cellID thisCell;
		// calculate cellID
		thisCell.phiCellID =  ...;
		thisCell.zCellID =  ...;
		thisCell.energyRatio =  1.; // no energy sharing with neighbour
		return thisCell;
	      }        



	      //let z_pad=z0, ignore the motoin on z irection
	      int chan=GetChanId(LposZ,phi_rad);
	      //if chan=-1 it means this channel is unreconstrutable
	      int tdc=int(t_s2pad/NS_PER_TIC+TPC_TZERO);
	      int adc=int(Kev2ADC*DiffEdep*1000.);

	      // Calculation from DriftEMagboltz::DriftEl2Pad()
	      //output: chan,adc,tdc (in nano second, not tic unit)
	      //        tdc = (t_s2gem1+toff)/tic+tzero, in tic unit
	      
	      dgtz["hitn"]   = hitn;
	      dgtz["CellID"] = chan;
	      dgtz["ADC"]    = adc;
	      dgtz["TDCL"]   = tdc;
	      dgtz["TDCR"]   = tdc;
	
	      return dgtz;
	    }
	}
}

vector<identifier>  bonus_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	vector<identifier> yid = id;
	G4ThreeVector xyz = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector Lxyz = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()
	  ->GetTopTransform().TransformPoint(xyz);///< Local Coordinates of interaction

	int CellID     = yid[0].id;
	int phiCellID  = yid[1].id;
	int zCellID    = yid[2].id;

	
	vector<cellID> a_hit = bonus.findPads(Lxyz);

	/*  No need to define multi-hit
	vector<cellID> multi_hit = bonus.findPads(Lxyz);
        int n_multi_hits = multi_hit.size(); // this is how many cell are hit in this step
	// closest strip
	yid[2].id = (int) multi_hit[0];
	
	yid[0].id_sharing = multi_hit[1];
	yid[1].id_sharing = multi_hit[1];
	yid[2].id_sharing = multi_hit[1];

	// make a loop for additional pads
        for(unsigned i=0; i<n_multi_hits; i++)
        {        
	  cellID thisCellID;
	  thisCellID.layer = .. ;// dont change
	  thisCellID.phiCell =.. ;// this is the one you change to the shared hit
	  yid.push_back(thisCellID) = .. ;
        
        }
	*/
	
        // make sure identifiers multiply based on vector<cellID>
	yid[0] =  CellID;
        yid[1] =  phiCellID;
        yid[2] =  zCellID;
	yid.push_back(thisCellID);
 	
	//	id[id.size()-1].id_sharing = 1;
	return yid;
}


map< string, vector <int> >  bonus_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}






