// Scorer for ScoreProcesses
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "ScoreProcesses.hh"

#include "TsTrackInformation.hh"

#include "G4VProcess.hh"
#include "G4Track.hh"

ScoreProcesses::ScoreProcesses(TsParameterManager *pM, TsMaterialManager *mM, TsGeometryManager *gM, TsScoringManager *scM, TsExtensionManager *eM,
							   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	fNtuple->RegisterColumnI(&fParentId, "ParentID");
	fNtuple->RegisterColumnD(&fEnergy, "Energy", "MeV");
	fNtuple->RegisterColumnF(&fWeight, "Weight", "");
	fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
	fNtuple->RegisterColumnS(&fOriginProcessName, "Origin Process");
}

G4bool ScoreProcesses::ProcessHits(G4Step *, G4TouchableHistory *)
{
	return false;
}

void ScoreProcesses::UserHookForEndOfTrack(const G4Track *aTrack)
{
	if (fIsActive && aTrack->GetParentID() > 0)
	{
		fParentId = aTrack->GetParentID();
		fEnergy = aTrack->GetVertexKineticEnergy();
		fWeight = aTrack->GetWeight();
		fParticleType = aTrack->GetParticleDefinition()->GetPDGEncoding();
		fOriginProcessName = aTrack->GetCreatorProcess()->GetProcessName();
		fNtuple->Fill();
	}
	else
		fSkippedWhileInactive++;
}
