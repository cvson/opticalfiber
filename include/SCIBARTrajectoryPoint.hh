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
// $Id: SCIBARTrajectoryPoint.hh 72065 2013-07-05 09:54:59Z gcosmo $
//
/// \file optical/wls/include/SCIBARTrajectoryPoint.hh
/// \brief Definition of the SCIBARTrajectoryPoint class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SCIBARTrajectoryPoint_h_seen
#define SCIBARTrajectoryPoint_h_seen 1

#include "globals.hh"

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryPoint.hh"

#include "G4StepStatus.hh"

class G4Track;
class G4Step;
class G4VProcess;

class SCIBARTrajectoryPoint : public G4TrajectoryPoint {

//--------
  public: // without description
//--------

// Constructor/Destructor

    SCIBARTrajectoryPoint();
    SCIBARTrajectoryPoint(const G4Track* );
    SCIBARTrajectoryPoint(const G4Step* );
    SCIBARTrajectoryPoint(const SCIBARTrajectoryPoint &right);
    virtual ~SCIBARTrajectoryPoint();

// Operators

    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
    inline int operator==(const SCIBARTrajectoryPoint& right) const
    { return (this==&right); };

// Get/Set functions

    inline G4double GetTime() const { return fTime; };
    inline const G4ThreeVector GetMomentum() const { return fMomentum; };
    inline G4StepStatus GetStepStatus() const { return fStepStatus; };
    inline G4String GetVolumeName() const { return fVolumeName; };

// Get method for HEPRep style attributes

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
  private:
//---------

// Member data

    G4double fTime;
    G4ThreeVector fMomentum;
    G4StepStatus fStepStatus;
    G4String fVolumeName;

};

extern G4ThreadLocal G4Allocator<SCIBARTrajectoryPoint>* SCIBARTrajPointAllocator;

inline void* SCIBARTrajectoryPoint::operator new(size_t)
{
    if(!SCIBARTrajPointAllocator)
      SCIBARTrajPointAllocator = new G4Allocator<SCIBARTrajectoryPoint>;
    return (void *) SCIBARTrajPointAllocator->MallocSingle();
}

inline void SCIBARTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
    SCIBARTrajPointAllocator->FreeSingle(
        (SCIBARTrajectoryPoint *) aTrajectoryPoint);
}

#endif
