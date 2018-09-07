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
// $Id: SCIBARPhotonDetHit.hh 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/include/SCIBARPhotonDetHit.hh
/// \brief Definition of the SCIBARPhotonDetHit class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SCIBARPhotonDetHit_h
#define SCIBARPhotonDetHit_h 1

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
// SCIBARPhotonDetHit Class
//--------------------------------------------------

class SCIBARPhotonDetHit : public G4VHit
{
  public:

    SCIBARPhotonDetHit();
    SCIBARPhotonDetHit(G4ThreeVector pExit, G4ThreeVector pArrive, G4double pTime);
    virtual ~SCIBARPhotonDetHit();

    SCIBARPhotonDetHit(const SCIBARPhotonDetHit &right);
    const SCIBARPhotonDetHit& operator=(const SCIBARPhotonDetHit& right);

    G4int operator==(const SCIBARPhotonDetHit& right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    inline void SetArrivalPos(G4ThreeVector xyz) { fPosArrive = xyz; }
    inline G4ThreeVector GetArrivalPos() { return fPosArrive; }

    inline void SetExitPos(G4ThreeVector xyz) { fPosExit = xyz; }
    inline G4ThreeVector GetExitPos() { return fPosExit; }

    inline void SetArrivalTime(G4double t) { fArrivalTime = t; }
    inline G4double GetArrivalTime() { return fArrivalTime; }
 
  private:

    // the arrival time of the photon
    G4double      fArrivalTime;
    // where the photon hit the detector (detector's coordinate)
    G4ThreeVector fPosArrive;
    // where the photon exited the fiber (world's coordinate)
    G4ThreeVector fPosExit;

};

//--------------------------------------------------
// Type Definitions
//--------------------------------------------------

typedef G4THitsCollection<SCIBARPhotonDetHit> SCIBARPhotonDetHitsCollection;

extern G4ThreadLocal G4Allocator<SCIBARPhotonDetHit>* SCIBARPhotonDetHitAllocator;

//--------------------------------------------------
// Operator Overloads
//--------------------------------------------------

inline void* SCIBARPhotonDetHit::operator new(size_t)
{
  if(!SCIBARPhotonDetHitAllocator)
      SCIBARPhotonDetHitAllocator = new G4Allocator<SCIBARPhotonDetHit>;
  return (void *) SCIBARPhotonDetHitAllocator->MallocSingle();
}

inline void SCIBARPhotonDetHit::operator delete(void *aHit)
{
  SCIBARPhotonDetHitAllocator->FreeSingle((SCIBARPhotonDetHit*) aHit);
}

#endif
