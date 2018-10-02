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
// $Id: SILICAPhotonDetHit.hh 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/include/SILICAPhotonDetHit.hh
/// \brief Definition of the SILICAPhotonDetHit class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SILICAPhotonDetHit_h
#define SILICAPhotonDetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

class G4VTouchable;

//--------------------------------------------------
// SILICAPhotonDetHit Class
//--------------------------------------------------

class SILICAPhotonDetHit : public G4VHit
{
  public:

    SILICAPhotonDetHit();
    SILICAPhotonDetHit(G4ThreeVector pExit, G4ThreeVector pArrive, G4double pTime, G4double pEnergyDeposit);
    virtual ~SILICAPhotonDetHit();

    SILICAPhotonDetHit(const SILICAPhotonDetHit &right);
    const SILICAPhotonDetHit& operator=(const SILICAPhotonDetHit& right);

    G4int operator==(const SILICAPhotonDetHit& right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    inline void SetArrivalPos(G4ThreeVector xyz) { fPosArrive = xyz; }
    inline G4ThreeVector GetArrivalPos() { return fPosArrive; }

    inline void SetExitPos(G4ThreeVector xyz) { fPosExit = xyz; }
    inline G4ThreeVector GetExitPos() { return fPosExit; }

    inline void SetArrivalTime(G4double t) { fArrivalTime = t; }
    inline G4double GetArrivalTime() { return fArrivalTime; }
    
    G4double GetEnergyDeposit() {return fenergyDeposit;};
 
  private:

    // the arrival time of the photon
    G4double      fArrivalTime;
    // where the photon hit the detector (detector's coordinate)
    G4ThreeVector fPosArrive;
    // where the photon exited the fiber (world's coordinate)
    G4ThreeVector fPosExit;
    
    G4double fenergyDeposit;

};

//--------------------------------------------------
// Type Definitions
//--------------------------------------------------

typedef G4THitsCollection<SILICAPhotonDetHit> SILICAPhotonDetHitsCollection;

extern G4ThreadLocal G4Allocator<SILICAPhotonDetHit>* SILICAPhotonDetHitAllocator;

//--------------------------------------------------
// Operator Overloads
//--------------------------------------------------

inline void* SILICAPhotonDetHit::operator new(size_t)
{
  if(!SILICAPhotonDetHitAllocator)
      SILICAPhotonDetHitAllocator = new G4Allocator<SILICAPhotonDetHit>;
  return (void *) SILICAPhotonDetHitAllocator->MallocSingle();
}

inline void SILICAPhotonDetHit::operator delete(void *aHit)
{
  SILICAPhotonDetHitAllocator->FreeSingle((SILICAPhotonDetHit*) aHit);
}

#endif
