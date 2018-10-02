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
/// \file analysis/WLSROOT/src/SILICAHistoManager.cc
/// \brief Implementation of the SILICAHistoManager class
//
//
// $Id: SILICAHistoManager.cc 92443 2015-09-01 13:56:16Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "SILICAHistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAHistoManager::SILICAHistoManager()
 : fFactoryOn(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAHistoManager::~SILICAHistoManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in SILICAHistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  analysisManager->SetNtupleDirectoryName("ntuple");
    
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile("WLSROOT");
  if (! fileOpen) {
    G4cerr << "\n---> SILICAHistoManager::Book(): cannot open " 
           << analysisManager->GetFileName() << G4endl;
    return;
  }
  
  // Create histograms.
  // Histogram ids are generated automatically starting from 0.
  // The start value can be changed by:
  // analysisManager->SetFirstHistoId(1);  
  
  // id = 0
  analysisManager->CreateH1("Edep","Edep in scintillator (MeV)", 200, 0., 50*MeV);
  // id = 1
  analysisManager->CreateH1("Ldep","trackL in scintillator (mm)", 200, 0., 10*cm);
     // id = 2
  analysisManager->CreateH1("dEdx","dE/dx (MeV/cm)", 200, 0., 10*MeV/cm);
     // id = 3
    analysisManager->CreateH1("nhit","Number of hits", 200, 0., 1000);
    // id = 4
    analysisManager->CreateH1("nopticalphoton","Number of Optical Photon", 200, 0., 50000);
    //id =5
    analysisManager->CreateH1("Edephit","Edep in scintillator from hits(MeV)", 200, 0., 10*MeV);
  // Create ntuples.
  // Ntuples ids are generated automatically starting from 0.
  // The start value can be changed by:
  // analysisManager->SetFirstMtupleId(1);  
  
  // Create 1st ntuple (id = 0)
  analysisManager->CreateNtuple("ELdep", "ELdep");
  analysisManager->CreateNtupleDColumn("Edep"); // column Id = 0
  analysisManager->CreateNtupleDColumn("Ldep"); // column Id = 1
  analysisManager->FinishNtuple();

  // Create 2nd ntuple (id = 1)
  //    
  analysisManager->CreateNtuple("PrimTrack", "Primarytrack");
  analysisManager->CreateNtupleDColumn("xpos"); // column Id = 0
  analysisManager->CreateNtupleDColumn("ypos"); // column Id = 1
  analysisManager->CreateNtupleDColumn("zpos"); // column Id = 2
  analysisManager->CreateNtupleDColumn("pz"); // column Id = 3
  analysisManager->CreateNtupleDColumn("theta"); // column Id = 4
    analysisManager->CreateNtupleDColumn("nphotonpass2end"); // column Id = 5//newadd2
    analysisManager->CreateNtupleDColumn("nphotonfail2end"); // column Id = 6//newadd2
    analysisManager->CreateNtupleDColumn("nphotonscintillation"); //column Id=7 newadd5
    analysisManager->CreateNtupleDColumn("nphotoncherenkov"); //column Id=8 newadd5
    analysisManager->CreateNtupleDColumn("totalenscintillation"); //column Id=9 newadd5
    analysisManager->CreateNtupleDColumn("totalencherenkov"); //column Id=10 newadd5
  analysisManager->FinishNtuple();
  
  fFactoryOn = true;       

  G4cout << "\n----> Output file is open in " 
         << analysisManager->GetFileName() << "." 
         << analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::Save()
{
  if (! fFactoryOn) return;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->Write();
  analysisManager->CloseFile(); 
   
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
      
  delete G4AnalysisManager::Instance();
  fFactoryOn = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  analysisManager->FillH1(ih, xbin, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::Normalize(G4int ih, G4double fac)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  G4H1* h1 = analysisManager->GetH1(ih);
  if (h1) h1->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::FillNtuple(G4double energyDep,  G4double trackLDep)
{                
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 0)
  analysisManager->FillNtupleDColumn(0, 0, energyDep);
  analysisManager->FillNtupleDColumn(0, 1, trackLDep);
  analysisManager->AddNtupleRow(0);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//newadd
void SILICAHistoManager::FillPrimTrack(G4double xval, G4double yval, G4double zval, G4double pzval, G4double thetaval, G4int nphotopassval, G4int nphotofailval, G4int nphotoscintval, G4int nphotocherenval, G4double enscintval, G4double encherenval)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Fill 2nd ntuple ( id = 1)
    analysisManager->FillNtupleDColumn(1, 0, xval);//primary track x pos
    analysisManager->FillNtupleDColumn(1, 1, yval);//primary track y pos
    analysisManager->FillNtupleDColumn(1, 2, zval);//primary track z pos
    analysisManager->FillNtupleDColumn(1, 3, pzval);//primary track pz
    analysisManager->FillNtupleDColumn(1, 4, thetaval);//primary track theta
    analysisManager->FillNtupleDColumn(1, 5, nphotopassval);//no of photon go to sensitive detector
    analysisManager->FillNtupleDColumn(1, 6, nphotofailval);//no of photon not go to sensitive detector
    analysisManager->FillNtupleDColumn(1, 7, nphotoscintval);//no of photon from scintillation process
    analysisManager->FillNtupleDColumn(1, 8, nphotocherenval);//no of photon from cherenkov process
    analysisManager->FillNtupleDColumn(1, 9, enscintval);//energy of all scintillation photons
    analysisManager->FillNtupleDColumn(1, 10, encherenval);//energy of all cherenkov photons
    analysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAHistoManager::PrintStatistic()
{
  if (! fFactoryOn) return;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  G4cout << "\n ----> print histograms statistic \n" << G4endl;
  for ( G4int i=0; i<analysisManager->GetNofH1s(); ++i ) {
    G4String name = analysisManager->GetH1Name(i);
    G4H1* h1 = analysisManager->GetH1(i);
    
    G4String unitCategory;
    if (name[0U] == 'E' ) unitCategory = "Energy"; 
    if (name[0U] == 'L' ) unitCategory = "Length";
      if (name[0U] == 'd' ) unitCategory = "Energy/Length";
      if (name[0U] == 'n' ) continue;
         // we use an explicit unsigned int type for operator [] argument
         // to avoid problems with windows compiler

    G4cout << name
           << ": mean = " << G4BestUnit(h1->mean(), unitCategory) 
           << " rms = " << G4BestUnit(h1->rms(), unitCategory ) 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


