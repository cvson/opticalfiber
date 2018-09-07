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
// $Id: SCIBARPrimaryGeneratorAction.cc 88760 2015-03-09 12:21:59Z gcosmo $
//
/// \file optical/wls/src/SCIBARPrimaryGeneratorAction.cc
/// \brief Implementation of the SCIBARPrimaryGeneratorAction class
//
//
#include "G4ios.hh"
#include "G4Event.hh"

#include "G4GeneralParticleSource.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"

#include "Randomize.hh"

#include "SCIBARPrimaryGeneratorAction.hh"

#include "SCIBARDetectorConstruction.hh"
#include "SCIBARPrimaryGeneratorMessenger.hh"

#include "G4SystemOfUnits.hh"

#include "G4AutoLock.hh"

namespace {
  G4Mutex gen_mutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SCIBARPrimaryGeneratorAction::fFirst = false;

SCIBARPrimaryGeneratorAction::
                         SCIBARPrimaryGeneratorAction(SCIBARDetectorConstruction* dc)
{
  fDetector = dc;
  fIntegralTable = NULL;

  fParticleGun = new G4GeneralParticleSource();
 
  fGunMessenger = new SCIBARPrimaryGeneratorMessenger(this);

  // G4String particleName;
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  fTimeConstant = 0.;

//  fParticleGun->SetParticleDefinition(particleTable->
//                               FindParticle(particleName="opticalphoton"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SCIBARPrimaryGeneratorAction::~SCIBARPrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  if (fIntegralTable) {
     fIntegralTable->clearAndDestroy();
     delete fIntegralTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::SetDecayTimeConstant(G4double time)
{
  fTimeConstant = time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::BuildEmissionSpectrum()
{
   if (fIntegralTable) return;

   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

   G4int numOfMaterials = G4Material::GetNumberOfMaterials();

   if(!fIntegralTable)fIntegralTable = new G4PhysicsTable(numOfMaterials);

   for (G4int i=0 ; i < numOfMaterials; i++) {

       G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
                                           new G4PhysicsOrderedFreeVector();

       G4Material* aMaterial = (*theMaterialTable)[i];

       G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                    aMaterial->GetMaterialPropertiesTable();

       if (aMaterialPropertiesTable) {
          G4MaterialPropertyVector* theSCIBARVector =
                      aMaterialPropertiesTable->GetProperty("SCIBARCOMPONENT");

          if (theSCIBARVector) {
             G4double currentIN = (*theSCIBARVector)[0];
             if (currentIN >= 0.0) {
                G4double currentPM = theSCIBARVector->Energy(0);
                G4double currentCII = 0.0;
                aPhysicsOrderedFreeVector->
                                     InsertValues(currentPM , currentCII);
                G4double prevPM  = currentPM;
                G4double prevCII = currentCII;
                G4double prevIN  = currentIN;

                for (size_t j = 1;
                     j < theSCIBARVector->GetVectorLength();
                     j++)
                {
                  currentPM = theSCIBARVector->Energy(j);
                  currentIN = (*theSCIBARVector)[j];
                  currentCII = 0.5 * (prevIN + currentIN);
                  currentCII = prevCII + (currentPM - prevPM) * currentCII;
                  aPhysicsOrderedFreeVector->
                                     InsertValues(currentPM, currentCII);
                  prevPM  = currentPM;
                  prevCII = currentCII;
                  prevIN  = currentIN;
                }
             }
          }
       }
       fIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (!fFirst) {
     fFirst = true;
     BuildEmissionSpectrum();
  }

#ifdef use_sampledEnergy
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4double sampledEnergy = 3*eV;

  for (size_t j=0 ; j<theMaterialTable->size() ; j++) {
      G4Material* fMaterial = (*theMaterialTable)[j];
      if (fMaterial->GetName() == "PMMA" ) {
         G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                      fMaterial->GetMaterialPropertiesTable();
         const G4MaterialPropertyVector* SCIBARIntensity =
                   aMaterialPropertiesTable->GetProperty("SCIBARCOMPONENT");

         if (SCIBARIntensity) {
            G4int MaterialIndex = fMaterial->GetIndex();
            G4PhysicsOrderedFreeVector* SCIBARIntegral =
              (G4PhysicsOrderedFreeVector*)((*fIntegralTable)(MaterialIndex));

            G4double CIImax = SCIBARIntegral->GetMaxValue();
            G4double CIIvalue = G4UniformRand()*CIImax;

            sampledEnergy = SCIBARIntegral->GetEnergy(CIIvalue);
         }
      }
  }

  //fParticleGun->SetParticleEnergy(sampledEnergy);
#endif

  //The code behing this line is not thread safe because polarization
  //and time are randomly selected and GPS properties are global
  G4AutoLock l(&gen_mutex);
  if(fParticleGun->GetParticleDefinition()->GetParticleName()=="opticalphoton"){
    SetOptPhotonPolar();
    SetOptPhotonTime();
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0*deg;
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
  {
     G4cout << "-> warning from SCIBARPrimaryGeneratorAction::SetOptPhotonPolar()"
            << ":  the ParticleGun is not an opticalphoton" << G4endl;
     return;
  }

  G4ThreeVector normal (1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2       = product*product;

  G4ThreeVector e_perpend (0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
  G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

  G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
  fParticleGun->SetParticlePolarization(polar);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SCIBARPrimaryGeneratorAction::SetOptPhotonTime()
{
   G4double time = -std::log(G4UniformRand())*fTimeConstant;
   fParticleGun->SetParticleTime(time);
}
