//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SILICATrajectory.hh 90240 2015-05-21 09:08:13Z gcosmo $
//
/// \file optical/wls/include/SILICATrajectory.hh
/// \brief Definition of the SILICATrajectory class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SILICATrajectory_h_seen
#define SILICATrajectory_h_seen 1

#include <vector>
#include <stdlib.h>

#include "globals.hh"

#include "G4ios.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VTrajectory.hh"
#include "G4ParticleDefinition.hh"
#include "G4TrajectoryPoint.hh"

typedef std::vector<G4VTrajectoryPoint*> SILICATrajectoryPointContainer;

class SILICATrajectory : public G4VTrajectory
{

//--------
   public: // without description
//--------

// Constructor/Destructor

     SILICATrajectory();
     SILICATrajectory(const G4Track* );
     SILICATrajectory(SILICATrajectory &);
     virtual ~SILICATrajectory();

// Operators

     inline void* operator new(size_t);
     inline void  operator delete(void*);
     inline int operator == (const SILICATrajectory& right) const
     { return (this==&right); }

// Get/Set functions

     inline virtual G4int GetTrackID() const { return fTrackID; }
     inline virtual G4int GetParentID() const { return fParentID; }
     inline virtual G4String GetParticleName() const { return fParticleName; }
     inline virtual G4double GetCharge() const { return fPDGCharge; }
     inline virtual G4int GetPDGEncoding() const { return fPDGEncoding; }
     inline virtual G4ThreeVector GetInitialMomentum() const 
                                                    {return fInitialMomentum;}

// Other member functions

     virtual void ShowTrajectory(std::ostream& os=G4cout) const;
     virtual void AppendStep(const G4Step* aStep);
     virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

     G4ParticleDefinition* GetParticleDefinition();

     virtual int GetPointEntries() const
     { return fpPointsContainer->size(); }
     virtual G4VTrajectoryPoint* GetPoint(G4int i) const
     { return (*fpPointsContainer)[i]; }

    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

// Member data

     SILICATrajectoryPointContainer* fpPointsContainer;

     G4int fTrackID;
     G4int fParentID;
     G4double fPDGCharge;
     G4int    fPDGEncoding;
     G4String fParticleName;
     G4ThreeVector fInitialMomentum;

     G4ParticleDefinition* fParticleDefinition;

};

extern G4ThreadLocal G4Allocator<SILICATrajectory>* SILICATrajectoryAllocator;

inline void* SILICATrajectory::operator new(size_t) {
    if(!SILICATrajectoryAllocator)
      SILICATrajectoryAllocator = new G4Allocator<SILICATrajectory>;
    return (void*) SILICATrajectoryAllocator->MallocSingle();
}

inline void SILICATrajectory::operator delete(void* aTrajectory) {
    SILICATrajectoryAllocator->FreeSingle((SILICATrajectory*)aTrajectory);
}

#endif
