// Scorer for TsTrackTerminator

#include "TsTrackTerminator.hh"

#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeTable.hh"
#include "G4RunManager.hh"

TsTrackTerminator::TsTrackTerminator(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
									 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
	  fPm(pM), fNbOfMoleculesToScavenge(), fEnergyLossKill(DBL_MAX), fEnergyLossAbort(DBL_MAX), fMaximumTrackLength(DBL_MAX)
{
	ResolveParameters();

	fEnergyLoss = 0.0;
	fTotalTrackLength = 0.0;
}

void TsTrackTerminator::ResolveParameters()
{
	if (fPm->ParameterExists(GetFullParmName("KillPrimaryIfEnergyLossExceeds")))
		fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"), "Energy");
	if (fPm->ParameterExists(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds")))
		fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"), "Energy");
	if (fPm->ParameterExists(GetFullParmName("KillPrimaryBasedOnTrackLength")))
		fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");

	if (fPm->ParameterExists(GetFullParmName("ScavengeTheseMolecules"))) {
		fNbOfMoleculesToScavenge = fPm->GetVectorLength(GetFullParmName("ScavengeTheseMolecules"));

		// Use EITHER ScavengingCapacities OR ScavengerConcentration + ScavengerReactionRate
		if (fPm->ParameterExists(GetFullParmName("ScavengingCapacities")) &&
			(fPm->ParameterExists(GetFullParmName("ScavengerConcentration")) || fPm->ParameterExists(GetFullParmName("ScavengerReactionRate")))) {
			G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
			G4cerr << "You may specify only one of the following set of parameters to describe the scavenging reactions:" << G4endl;
			G4cerr << " (1) " << GetFullParmName("ScavengingCapacities") << G4endl;
			G4cerr << " (2) " << GetFullParmName("ScavengerConcentration") << " and " << GetFullParmName("ScavengerReactionRate") << G4endl;
			fPm->AbortSession(1);
		}
		if ((fPm->ParameterExists(GetFullParmName("ScavengingCapacities")) && fPm->GetVectorLength(GetFullParmName("ScavengingCapacities")) != fNbOfMoleculesToScavenge) ||
			(fPm->ParameterExists(GetFullParmName("ScavengerConcentration")) && fPm->GetVectorLength(GetFullParmName("ScavengerConcentration")) != fNbOfMoleculesToScavenge) ||
			(fPm->ParameterExists(GetFullParmName("ScavengerReactionRate")) && fPm->GetVectorLength(GetFullParmName("ScavengerReactionRate")) != fNbOfMoleculesToScavenge))
		{
			G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
			G4cerr << "The number of elements in " << GetFullParmName("ScavengingCapacities") << G4endl;
			G4cerr << " or " << GetFullParmName("ScavengerConcentration") << " and " << GetFullParmName("ScavengerReactionRate") << G4endl;
			G4cerr << " must match those in " << GetFullParmName("ScavengeTheseMolecules") << G4endl;
			fPm->AbortSession(1);
		}

		G4String* moleculeNames = fPm->GetStringVector(GetFullParmName("ScavengeTheseMolecules"));
		for (auto i = 0; i < fNbOfMoleculesToScavenge; ++i) {
			G4MolecularConfiguration* config = G4MoleculeTable::Instance()->GetConfiguration(moleculeNames[i]);
			fScavengeTheseMolecules.push_back(config->GetMoleculeID());
		}
		delete[] moleculeNames;

		if (fPm->ParameterExists(GetFullParmName("ScavengingCapacities"))) {
			G4double* scavengingCapacities = fPm->GetDoubleVector(GetFullParmName("ScavengingCapacities"), "perTime");
			fScavengingCapacities = std::vector<G4double>(scavengingCapacities, scavengingCapacities + fNbOfMoleculesToScavenge);
			delete[] scavengingCapacities;
		}
		else {
			G4double* scavengerConcentrations = fPm->GetDoubleVector(GetFullParmName("ScavengerConcentrations"), "molar concentration");
			G4double* scavengerReactionRates = fPm->GetDoubleVector(GetFullParmName("ScavengerReactionRates"), "perMolarConcentration perTime");
			for (auto i = 0; i < fNbOfMoleculesToScavenge; ++i) {
				fScavengingCapacities.push_back(scavengerConcentrations[i] * scavengerReactionRates[i]);
			}
			delete[] scavengerReactionRates;
			delete[] scavengerConcentrations;
		}
	}
}

G4bool TsTrackTerminator::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	G4Track* aTrack = aStep->GetTrack();
	if (aTrack->GetTrackID() > -1)
	{
		if (aTrack->GetParentID() == 0)
		{
			fEnergyLoss += aStep->GetPreStepPoint()->GetKineticEnergy() -
						   aStep->GetPostStepPoint()->GetKineticEnergy();
			fTotalTrackLength += aStep->GetStepLength();

			if (fEnergyLoss > fEnergyLossAbort) {
				G4cout << " -- Aborting event " << GetEventID() << G4endl;
				G4RunManager::GetRunManager()->AbortEvent();
				return true;
			}
			if (fEnergyLoss >= fEnergyLossKill) {
				G4cout << " -- Killing primary track of event " << GetEventID() << G4endl;
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				return true;
			}
			if (fTotalTrackLength >= fMaximumTrackLength) {
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				G4cout << " Track killed with track: " << fTotalTrackLength / nm << " nm with energy lost " << fEnergyLoss / keV << " keV " << G4endl;
				return true;
			}
		}
	}
	// else // TrackID < 0
	// {
	// 	if (fNbOfMoleculesToScavenge > 0)
	// 	{
	// 		G4Molecule* molecule = GetMolecule(aTrack);
	// 		G4int moleculeID = molecule->GetMoleculeID();
	// 		G4double t = aTrack->GetGlobalTime();
	// 		for (auto i = 0; i < fNbOfMoleculesToScavenge && fScavengeTheseMolecules[i] == moleculeID; ++i) {
	// 			G4double probability = 1. - std::exp(-fScavengingCapacities[i] * aTrack->GetGlobalTime());
	// 			if (G4UniformRand() < probability) {
	// 				if (fVerbosity > 0) {
	// 					G4cout << "TsTrackTerminator: Scavenged " << molecule->GetName() << " at t=" << aTrack->GetGlobalTime() / ns << " ns" << G4endl;
	// 				}
	// 				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	// 				return true;
	// 			}
	// 		}
	// 	}
	// }

	return false;
}

void TsTrackTerminator::UserHookForChemicalStep(const G4Step* aStep)
{
	auto aTrack = aStep->GetTrack();

	if (aTrack->GetTrackID() > -1) // this is a physical track
		return;

	if (fNbOfMoleculesToScavenge <= 0)
		return;

	if (!GetFilter()->Accept(aStep))
		return;

	G4Molecule* molecule = GetMolecule(aTrack);
	G4int moleculeID = molecule->GetMoleculeID();
	for (auto i = 0; i < fNbOfMoleculesToScavenge && fScavengeTheseMolecules[i] == moleculeID; ++i) {
		G4double probability = 1. - std::exp(-fScavengingCapacities[i] * aTrack->GetLocalTime());
		if (G4UniformRand() < probability) {
			if (fVerbosity > 0 || fTrackingVerbosity > 0) {
				G4cout << "TsTrackTerminator: Scavenged " << molecule->GetName() << " at t=" << aTrack->GetLocalTime() / ns
					   << " ns [global t=" << aTrack->GetGlobalTime() / ns << " ns]"
					   << G4endl;
			}
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return;
		}
	}
}

void TsTrackTerminator::AccumulateEvent()
{
	if (fHaveIncidentParticle) {
		UserHookForEndOfIncidentParticle();
		fHaveIncidentParticle = false;
	}

	UserHookForEndOfEvent();

	fScoredHistories++;
}

void TsTrackTerminator::UserHookForEndOfEvent()
{
	fTotalTrackLength = 0.0;
	fEnergyLoss = 0.0;
}
