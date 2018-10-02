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
// $Id: SILICAEventAction.hh 90240 2015-05-21 09:08:13Z gcosmo $
//
/// \file optical/wls/include/SILICAEventAction.hh
/// \brief Definition of the SILICAEventAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef SILICAEventAction_h
#define SILICAEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class SILICARunAction;
class SILICAHistoManager;//newadd
class SILICAEventActionMessenger;

class SILICAEventAction : public G4UserEventAction
{
  public:

    //SILICAEventAction(SILICARunAction*);
    SILICAEventAction(SILICARunAction*, SILICAHistoManager*);
    virtual ~SILICAEventAction();

  public:

    virtual void   BeginOfEventAction(const G4Event*);
    virtual void     EndOfEventAction(const G4Event*);
    void AddDep(G4double de, G4double dl) {fEnergyDep += de; fTrackLDep += dl;};//newadd
    void AddNoPhotonType(G4int nphotonscint, G4int nphotoncheren) {fnphotonscint += nphotonscint; fnphotonscheren += nphotoncheren;};//newadd6
    void AddPhotonEnergy(G4double enscint, G4double encheren){fenscint +=enscint; fencheren +=encheren;};
    void AddNoPhotonPass2end(G4int nphotonpass) {fnphotonpass2end += nphotonpass;};//newadd6
    void AddNoPhotonFail2end(G4int nphotonfail) {fnphotonfail2end += nphotonfail;};//newadd6
    void AddPriTrackx(G4double xval){fprix = xval;};
    void AddPriTracky(G4double yval){fpriy = yval;};
    void AddPriTrackz(G4double zval){fpriz = zval;};
    void AddPriTrackpz(G4double pzval){fpripz = pzval;};
    void AddPriTracktheta(G4double thetaval){fpritheta = thetaval;};
    //void AddPriTrackNophotonpass2end(G4int nphoton2passval){fnphotonpass2end = nphoton2passval;};//newadd2
    //void AddPriTrackNophotonfail2end(G4int nphoton2failval){fnphotonfail2end = nphoton2failval;};//newadd2
    //void AddPriTrackNophotonScint(G4int nphotonscintval){fnphotonscint = nphotonscintval;};//newadd5
    //void AddPriTrackNophotonCheren(G4int nphotoncherenval){fnphotonscheren = nphotoncherenval;};//newadd5
    G4int GetEventNo();
    void SetEventVerbose(G4int);

  private:

    SILICARunAction* fRunAction;
    SILICAHistoManager* fHistoManager;//newadd
    SILICAEventActionMessenger* fEventMessenger;

    G4int fVerboseLevel;
 
    G4int fMPPCCollID;
	 G4double  fEnergyDep;//newadd
   G4double  fTrackLDep;//newadd
    G4double  fEnergyDepHit;//newadd3
    G4double fprix, fpriy, fpriz, fpripz, fpritheta;//newadd
    G4int fnphotonpass2end, fnphotonfail2end;//newadd2
    G4int fnphotonscint, fnphotonscheren;
    G4double fenscint, fencheren;
};

#endif
