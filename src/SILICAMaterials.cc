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
// $Id: SILICAMaterials.cc 82854 2014-07-14 09:08:25Z gcosmo $
//
/// \file optical/wls/src/SILICAMaterials.cc
/// \brief Implementation of the SILICAMaterials class
//
//
#include "SILICAMaterials.hh"

#include "G4SystemOfUnits.hh"

SILICAMaterials* SILICAMaterials::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAMaterials::SILICAMaterials()
{
    fNistMan = G4NistManager::Instance();
    
    fNistMan->SetVerbose(2);
    
    CreateMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAMaterials::~SILICAMaterials()
{
    delete    fPMMA;
    delete    fPethylene;
    delete    fFPethylene;
    delete    fPolystyrene;
    delete    fSilicone;
    delete  fSilica;
    delete  fTECS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SILICAMaterials* SILICAMaterials::GetInstance()
{
    if (fInstance == 0)
    {
        fInstance = new SILICAMaterials();
    }
    return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* SILICAMaterials::GetMaterial(const G4String material)
{
    G4Material* mat =  fNistMan->FindOrBuildMaterial(material);
    
    if (!mat) mat = G4Material::GetMaterial(material);
    if (!mat) {
        std::ostringstream o;
        o << "Material " << material << " not found!";
        G4Exception("SILICAMaterials::GetMaterial","",
                    FatalException,o.str().c_str());
    }
    
    return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SILICAMaterials::CreateMaterials()
{
    G4double density;
    G4int ncomponents;
    G4double fractionmass;
    std::vector<G4int> natoms;
    std::vector<G4double> fractionMass;
    std::vector<G4String> elements;
    
    // Materials Definitions
    // =====================
    
    //--------------------------------------------------
    // Vacuum
    //--------------------------------------------------
    
    fNistMan->FindOrBuildMaterial("G4_Galactic");
    
    //--------------------------------------------------
    // Air
    //--------------------------------------------------
    
    fAir = fNistMan->FindOrBuildMaterial("G4_AIR");
    
    //--------------------------------------------------
    // Silica
    //--------------------------------------------------
    G4Element* elSi  = new G4Element("Silicon"  ,"Si" , 14, 28.085*g/mole);
    G4Element* elO  = new G4Element("Oxygen"  ,"O" , 8, 16.00*g/mole);
    fSilica = new G4Material("Silica",2.203*g/cm3,2);
    fSilica->AddElement(elSi,0.467);
    fSilica->AddElement(elO,1.0-0.467);
    
    //--------------------------------------------------
    // TECS
    //--------------------------------------------------
    G4Element* elC  = new G4Element("Carbon"  ,"C" , 6, 12.011*g/mole);
    G4Element* elF  = new G4Element("Flourine"  ,"F" , 9, 18.998*g/mole);
    fTECS = new G4Material("TECS",1.8*g/cm3,2);
    fTECS->AddElement(elC,0.30);
    fTECS->AddElement(elF,0.7);
    //--------------------------------------------------
    // SILICAfiber PMMA
    //--------------------------------------------------
    
    elements.push_back("C");     natoms.push_back(5);
    elements.push_back("H");     natoms.push_back(8);
    elements.push_back("O");     natoms.push_back(2);
    
    density = 1.190*g/cm3;
    
    fPMMA = fNistMan->
    ConstructNewMaterial("PMMA", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Cladding (polyethylene)
    //--------------------------------------------------
    
    elements.push_back("C");     natoms.push_back(2);
    elements.push_back("H");     natoms.push_back(4);
    
    density = 1.200*g/cm3;
    
    fPethylene = fNistMan->
    ConstructNewMaterial("Pethylene", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Double Cladding (fluorinated polyethylene)
    //--------------------------------------------------
    
    elements.push_back("C");     natoms.push_back(2);
    elements.push_back("H");     natoms.push_back(4);
    
    density = 1.400*g/cm3;
    
    fFPethylene = fNistMan->
    ConstructNewMaterial("FPethylene", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Polystyrene
    //--------------------------------------------------
    
    elements.push_back("C");     natoms.push_back(8);
    elements.push_back("H");     natoms.push_back(8);
    
    density = 1.050*g/cm3;
    
    fPolystyrene = fNistMan->
    ConstructNewMaterial("Polystyrene", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Silicone (Template for Optical Grease)
    //--------------------------------------------------
    
    elements.push_back("C");     natoms.push_back(2);
    elements.push_back("H");     natoms.push_back(6);
    
    density = 1.060*g/cm3;
    
    fSilicone = fNistMan->
    ConstructNewMaterial("Silicone", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Aluminium
    //--------------------------------------------------
    
    fNistMan->FindOrBuildMaterial("G4_Al");
    
    //--------------------------------------------------
    // TiO2
    //--------------------------------------------------
    
    elements.push_back("Ti");     natoms.push_back(1);
    elements.push_back("O");      natoms.push_back(2);
    
    density     = 4.26*g/cm3;
    
    G4Material* TiO2 = fNistMan->
    ConstructNewMaterial("TiO2", elements, natoms, density);
    
    elements.clear();
    natoms.clear();
    
    //--------------------------------------------------
    // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.
    //--------------------------------------------------
    
    density = 1.52*g/cm3;
    
    fCoating =
    new G4Material("Coating", density, ncomponents=2);
    
    fCoating->AddMaterial(TiO2,         fractionmass = 15*perCent);
    fCoating->AddMaterial(fPolystyrene, fractionmass = 85*perCent);
    
    //
    // ------------ Generate & Add Material Properties Table ------------
    //
    
    G4double photonEnergy[] =
    {2.0*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3.0*eV,3.2*eV,3.4*eV,3.6*eV,3.8*eV,4.0*eV,4.2*eV};
    
    const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);
    
    //--------------------------------------------------
    // Air
    //--------------------------------------------------
    const G4int NUMENTRIES = 12;
    G4double refractiveIndex[] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
        1.00, 1.00};
    
    assert(sizeof(refractiveIndex) == sizeof(photonEnergy));
    
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("RINDEX", photonEnergy, refractiveIndex, nEntries);
    
    fAir->SetMaterialPropertiesTable(mpt);
    
    //--------------------------------------------------
    //  PMMA for SILICAfibers
    //--------------------------------------------------
    //Index of refraction from https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
    G4double refractiveIndexSILICAfiber[] =
    { 1.4574,1.4594,1.4615,1.4637,1.4661,1.4687,1.4704,1.4746,1.4779,1.4814,
        1.4851,1.4892};
    
    //Based on NA of fiber of 0.39
    G4double rindexTECS[NUMENTRIES];
    for(int i=0; i<NUMENTRIES; i++) rindexTECS[i] = sqrt(pow(refractiveIndexSILICAfiber[i],2)-0.39*0.39);
    
    assert(sizeof(refractiveIndexSILICAfiber) == sizeof(photonEnergy));
    
    // Transmission curve from https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=123 for 1 cm thick. Transmission scaled by 0.992 to get positive absorption lengths
    G4double absSILICAfiber[] =
    {442.*cm,361.*cm,378.*cm,321.*cm,332.*cm,320.*cm,237.*cm,238.*cm,225.*cm,204.*cm,
        251.*cm,200.*cm};
    
    //Arbitrary
   // G4double absorptionTECS[NUMENTRIES] = {1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm,1.*cm};
    G4double absorptionTECS[NUMENTRIES] = {1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m,1.*m};
    
    assert(sizeof(absSILICAfiber) == sizeof(photonEnergy));
    
    
    
    // Add entries into properties table
    // Add entries into properties table
    G4MaterialPropertiesTable* mptsilica = new G4MaterialPropertiesTable();
    mptsilica->AddProperty("RINDEX",photonEnergy,refractiveIndexSILICAfiber,nEntries);
    mptsilica->AddProperty("ABSLENGTH",photonEnergy,absSILICAfiber,nEntries);
    
    
    //THIS can be WRONG, NEED to check
    G4double silicaFast[] =
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0};
    assert(sizeof(silicaFast) == sizeof(photonEnergy));
    
    // Add entries into properties table
    mptsilica->AddProperty("FASTCOMPONENT",photonEnergy, silicaFast,nEntries);
    mptsilica->AddConstProperty("SCINTILLATIONYIELD",10./keV);
    mptsilica->AddConstProperty("RESOLUTIONSCALE",1.0);
    mptsilica->AddConstProperty("FASTTIMECONSTANT", 2.4*ns);
    fSilica->SetMaterialPropertiesTable(mptsilica);
    // Set the Birks Constant for Silica
    fSilica->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
    
    
    
    // Add entries into properties table
    // Add entries into properties table
    G4MaterialPropertiesTable* mptTECS = new G4MaterialPropertiesTable();
    mptTECS->AddProperty("RINDEX",photonEnergy,rindexTECS,nEntries);
    mptTECS->AddProperty("ABSLENGTH",photonEnergy,absorptionTECS,nEntries);
    
    fTECS->SetMaterialPropertiesTable(mptTECS);
    
    //--------------------------------------------------
    //  Polyethylene
    //--------------------------------------------------
    
    G4double refractiveIndexClad1[] =
    { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
        1.49, 1.49};
    
    assert(sizeof(refractiveIndexClad1) == sizeof(photonEnergy));
    
    G4double absClad[] =
    {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
        20.0*m,20.0*m};
    
    assert(sizeof(absClad) == sizeof(photonEnergy));
    
    // Add entries into properties table
    G4MaterialPropertiesTable* mptClad1 = new G4MaterialPropertiesTable();
    mptClad1->AddProperty("RINDEX",photonEnergy,refractiveIndexClad1,nEntries);
    mptClad1->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);
    
    fPethylene->SetMaterialPropertiesTable(mptClad1);
    
    //--------------------------------------------------
    // Fluorinated Polyethylene
    //--------------------------------------------------
    
    G4double refractiveIndexClad2[] =
    { 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
        1.42, 1.42};
    
    assert(sizeof(refractiveIndexClad2) == sizeof(photonEnergy));
    
    // Add entries into properties table
    G4MaterialPropertiesTable* mptClad2 = new G4MaterialPropertiesTable();
    mptClad2->AddProperty("RINDEX",photonEnergy,refractiveIndexClad2,nEntries);
    mptClad2->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);
    
    fFPethylene->SetMaterialPropertiesTable(mptClad2);
    
    //--------------------------------------------------
    // Silicone
    //--------------------------------------------------
    
    G4double refractiveIndexSilicone[] =
    { 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46, 1.46,
        1.46, 1.46};
    
    assert(sizeof(refractiveIndexSilicone) == sizeof(photonEnergy));
    
    // Add entries into properties table
    G4MaterialPropertiesTable* mptSilicone = new G4MaterialPropertiesTable();
    mptSilicone->
    AddProperty("RINDEX",photonEnergy,refractiveIndexSilicone,nEntries);
    mptSilicone->AddProperty("ABSLENGTH",photonEnergy,absClad,nEntries);
    
    fSilicone->SetMaterialPropertiesTable(mptSilicone);
    
    //--------------------------------------------------
    //  Polystyrene
    //--------------------------------------------------
    
    G4double refractiveIndexPS[] =
    
    { 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58, 1.58,
        1.58, 1.58};
    
    assert(sizeof(refractiveIndexPS) == sizeof(photonEnergy));
    // NO NEED
    G4double absPS[] =
    {250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,250.*cm,
        250.*cm,250.*cm};
    
    assert(sizeof(absPS) == sizeof(photonEnergy));
    
    //NO NEED,can be wrong
    G4double scintilFast[] =
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0};
    
    assert(sizeof(scintilFast) == sizeof(photonEnergy));
    
    // Add entries into properties table
    G4MaterialPropertiesTable* mptPolystyrene = new G4MaterialPropertiesTable();
    mptPolystyrene->AddProperty("RINDEX",photonEnergy,refractiveIndexPS,nEntries);
    mptPolystyrene->AddProperty("ABSLENGTH",photonEnergy,absPS,nEntries);
    mptPolystyrene->
    AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,nEntries);
    mptPolystyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
    mptPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
    //mptPolystyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
    mptPolystyrene->AddConstProperty("FASTTIMECONSTANT", 2.4*ns);
    
    fPolystyrene->SetMaterialPropertiesTable(mptPolystyrene);
    
    // Set the Birks Constant for the Polystyrene scintillator
    
    fPolystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
    
}
