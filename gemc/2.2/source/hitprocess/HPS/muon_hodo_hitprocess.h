#ifndef HS_HITPROCESS_H
#define HS_HITPROCESS_H 1

// gemc headers
#include "HitProcess.h"

// Class definition
class muon_hodo_HitProcess : public HitProcess
{
	public:
	
		~muon_hodo_HitProcess(){;}
	
		// - integrateDgt: returns digitized information integrated over the hit
		map<string, double> integrateDgt(MHit*, int);
		
		// - multiDgt: returns multiple digitized information / hit
		map< string, vector <int> > multiDgt(MHit*, int);
		
		// The pure virtual method processID returns a (new) identifier
		// containing hit sharing information
		vector<identifier> processID(vector<identifier>, G4Step*, detector);
	
		// creates the HitProcess
		static HitProcess *createHitClass() {return new muon_hodo_HitProcess;}
};

#endif
