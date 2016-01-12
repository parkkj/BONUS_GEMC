class cellID
{
        int phiCellID;
        int zCellID;
        double energyRatio;  // fraction of total energy assigned to this cell   - e.g. 1 = 100% of total energy
};



vector<cellID> findStrip(G4Threevector xyz)
{

        cellID thisCell;

        // calculate cellID

        thisCell.phiCellID = ….
        thisCell. zCellID = ….
        thisCell. energyRatio = ….
        

        return thisCell;
}        


Add 2 more variables to the identifiers, set to 1
cellPhi manual 1  cellZ manual 1

these will be set by hit process routine


vector<identifier> processID (vector<identifier> id)
{
        vector<identifier> yid = id;   // original geant4 id - this is your sensitive gem foil, some id are empty you fill it

        your id has dimension 3: layer (what you called “which” in the perl script), phiCell, zCell

 vector< cellID > multi_hit = bonus.FindStrip( Lxyz); 
        multi_hit.size() << this is how many cell are hit in this step


        YOU create as many yid as multi_hit.size() 

        make a loop:
        for(unsigned i=0; i<multi_hit.size(); i++)
        {        
                cellID thisCellID;

                thisCellID.layer = // dont change
                
                thisCellID.phiCell // this is the one you change to the shared hit

                yid.push_back(thisCellID);
        
        }


        // make sure identifiers multiply based on vector<cellID>
        yid[0] = “which” in the perl script
        yid[1] =  phiCell
        yid[2] = zCell




         return yid;
}
