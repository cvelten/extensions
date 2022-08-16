// Extra Class for use by TsYScorer
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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id: TsTrackerHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file TsTrackerHit.cc
/// \brief Implementation of the TsTrackerHit class

#include "TsTrackerHit.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<TsTrackerHit>* TsTrackerHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TsTrackerHit::TsTrackerHit()
	: G4VHit(),
	  fTrackID(-1),
	  fEdep(0.),
	  fPos(G4ThreeVector()),
	  fParticleName(""),
	  fParticleFlag(-1),
	  fIncidentEnergy(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TsTrackerHit::~TsTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TsTrackerHit::TsTrackerHit(const TsTrackerHit& right)
	: G4VHit()
{
	fTrackID = right.fTrackID;
	fEdep = right.fEdep;
	fPos = right.fPos;
	fIncidentEnergy = right.fIncidentEnergy;
	fParticleFlag = right.fParticleFlag;
	fParticleName = right.fParticleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TsTrackerHit& TsTrackerHit::operator=(const TsTrackerHit& right)
{
	fTrackID = right.fTrackID;
	fEdep = right.fEdep;
	fPos = right.fPos;
	fIncidentEnergy = right.fIncidentEnergy;
	fParticleFlag = right.fParticleFlag;
	fParticleName = right.fParticleName;

	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsTrackerHit::operator==(const TsTrackerHit& right) const
{
	return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
