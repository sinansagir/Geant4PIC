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
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"

#include "G4AutoLock.hh"
G4Mutex mutexInPGA = G4MUTEX_INITIALIZER;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fProton = particleTable->FindParticle("proton");
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-20.*cm));
  fParticleGun->SetParticleDefinition(fProton);
    
  // define commands for this class
  //DefineCommands();

  // Reading histogram from file
  // And change particle gun setting
  G4AutoLock lock(&mutexInPGA);
  gRandom->SetSeed(0);
  
  fBeamSourceFile = new TFile("../B5_PIC/Data/Geant4input_QED_L_7um_proton_19.root");
  fBeamSourceReader.SetTree("ana",fBeamSourceFile);
  TTreeReaderValue<Double_t> theKE(fBeamSourceReader, "b_KE");
  //TTreeReaderValue<Double_t> theWe(fBeamSourceReader, "b_We");
  TTreeReaderValue<Int_t>    theNp(fBeamSourceReader, "b_Np");
  TTreeReaderValue<Double_t> thePx(fBeamSourceReader, "b_px");
  TTreeReaderValue<Double_t> thePy(fBeamSourceReader, "b_py");
  TTreeReaderValue<Double_t> thePz(fBeamSourceReader, "b_pz");
  nentries = fBeamSourceReader.GetEntries();
  fBeamSourceReader.SetEntry(1);
   
  fRandomizePrimary = false;
  fParticleGun->SetParticleDefinition(fProton);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
//   G4double envSizeXY = 0;
//   G4double envSizeZ = 0;
// 
//   if (!fEnvelopeBox)
//   {
//     G4LogicalVolume* envLV
//       = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
//     if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
//   }
// 
//   if ( fEnvelopeBox ) {
//     envSizeXY = fEnvelopeBox->GetXHalfLength()*2.;
//     envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
//   }  
//   else  {
//     G4ExceptionDescription msg;
//     msg << "Envelope volume of box shape not found.\n"; 
//     msg << "Perhaps you have changed geometry.\n";
//     msg << "The gun will be place at the center.";
//     G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()",
//      "MyCode0002",JustWarning,msg);
//   }
  
   // Generate kinetic energy and momentum
   fBeamSourceReader.Restart();
   fBeamSourceReader.SetTree("ana",fBeamSourceFile);
   TTreeReaderValue<Double_t> theKE(fBeamSourceReader, "b_KE");
   //TTreeReaderValue<Double_t> theWe(fBeamSourceReader, "b_We");
   TTreeReaderValue<Int_t>    theNp(fBeamSourceReader, "b_Np");
   TTreeReaderValue<Double_t> thePx(fBeamSourceReader, "b_px");
   TTreeReaderValue<Double_t> thePy(fBeamSourceReader, "b_py");
   TTreeReaderValue<Double_t> thePz(fBeamSourceReader, "b_pz");
   TRandom3 rand(0);
   Long64_t coin = rand.Uniform(0,nentries);
   fBeamSourceReader.SetEntry(coin);
   auto ke = (*theKE);
   //auto we = *theWe;
   auto np = *theNp;
   auto pz = *thePx;//Epoch laser beam direction is x, while Geant4 beam axis is z, so swap x and z!
   auto py = *thePy;
   auto px = *thePz;//Epoch laser beam direction is x, while Geant4 beam axis is z, so swap x and z!
   auto parVec = G4ThreeVector(px, py, pz).unit();

   std::cout<<"Coin:"<<coin;
   std::cout<<", ke:"<<ke;
   std::cout<<", px:"<<parVec.x();
   std::cout<<", py:"<<parVec.y();
   std::cout<<", pz:"<<parVec.z()<<std::endl;
   
   fParticleGun->SetNumberOfParticles(np);
   fParticleGun->SetParticleMomentumDirection(parVec);
   fParticleGun->SetParticleEnergy(ke * MeV);
   fParticleGun->GeneratePrimaryVertex(anEvent);
   std::cout<<"Nparticles:"<<fParticleGun->GetNumberOfParticles()<<std::endl;
   std::cout<<"Momentum:"<<fParticleGun->GetParticleMomentum()<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

