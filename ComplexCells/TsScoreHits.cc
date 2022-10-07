// Scorer for Hits
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreHits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4VProcess.hh"

#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

#include "G4Molecule.hh"

TsScoreHits::TsScoreHits(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
	  fNextTimeForTrack(), fHitsMap(), fEnergyDepositedMap()
{
	SetUnit("");

	fTimeCut = 1.0 * us;
	if (fPm->ParameterExists(GetFullParmName("TimeCut")))
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"), "Time");

	if (fPm->ParameterExists(GetFullParmName("TimesToRecord")))
	{
		auto times = fPm->GetDoubleVector(GetFullParmName("TimesToRecord"), "Time");
		for (auto i = 0; i < fPm->GetVectorLength(GetFullParmName("TimesToRecord")); ++i)
		{
			if (times[i] > fTimeCut)
			{
				G4cerr << GetName() << " specified a time to record greater than the time cut (default: 1 us)" << G4endl;
				G4cerr << "at which chemical tracks will be killed" << G4endl;
			}
			fTimesToRecord.push_back(times[i]);
		}

		std::sort(fTimesToRecord.begin(), fTimesToRecord.end());
		// std::reverse(fTimesToRecord.begin(), fTimesToRecord.end());
	}
	else {
		fTimesToRecord.push_back(fTimeCut);
	}

	fIncludeChemistry = false;
	fIncludePhysics = true;
	if (fPm->ParameterExists(GetFullParmName("IncludeChemicalTrack")))
		fIncludeChemistry = fPm->GetBooleanParameter(GetFullParmName("IncludeChemicalTrack"));
	if (fPm->ParameterExists(GetFullParmName("IncludePhysicalTrack")))
		fIncludePhysics = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalTrack"));

	fNtuple->RegisterColumnS(&fParticleName, "Particle Name");
	fNtuple->RegisterColumnS(&fVolumeName, "Volume Name");
	fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume Copy Number");
	fNtuple->RegisterColumnF(&fTime, "Time", "ps");
	fNtuple->RegisterColumnI(&fHits, "Hits");
	fNtuple->RegisterColumnF(&fEnergyDeposited, "Energy Deposited", "keV");
}

G4bool TsScoreHits::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	G4Track* aTrack = aStep->GetTrack();

	if (fIncludeChemistry && aTrack->GetTrackID() < 0) {
		auto globalTime = aStep->GetPreStepPoint()->GetGlobalTime();
		auto itNextTime = fNextTimeForTrack.emplace(aTrack->GetTrackID(), fTimesToRecord.front()).first;

		if (globalTime >= itNextTime->second && itNextTime->second > 0)
		{
			G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

			TsScoreHitsIndex idx;
			idx.IsMolecule = true;
			idx.ParticleName = GetMolecule(aTrack)->GetName();
			idx.VolumeName = touchable->GetVolume()->GetName();
			idx.VolumeCopyNumber = touchable->GetVolume()->GetCopyNo();
			idx.Time = fNextTimeForTrack[aTrack->GetTrackID()];

			++fHitsMap[idx];

			auto itNextTime = std::upper_bound(fTimesToRecord.begin(), fTimesToRecord.end(), globalTime);
			if (itNextTime == fTimesToRecord.end())
				fNextTimeForTrack[aTrack->GetTrackID()] = 0;
			else
				fNextTimeForTrack[aTrack->GetTrackID()] = *itNextTime;
		}
	}
	else if (fIncludePhysics && aStep->GetTotalEnergyDeposit() > 0) {
		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

		TsScoreHitsIndex idx;
		idx.ParticleName = aTrack->GetParticleDefinition()->GetParticleName();
		idx.VolumeName = touchable->GetVolume()->GetName();
		idx.VolumeCopyNumber = touchable->GetVolume()->GetCopyNo();

		++fHitsMap[idx];
		fEnergyDepositedMap[idx] += aStep->GetTotalEnergyDeposit();

		return true;
	}

	return false;
}

void TsScoreHits::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

	TsScoreHits* worker = dynamic_cast<TsScoreHits*>(workerScorer);

	for (auto it = worker->GetHitCountMap().cbegin();
		 it != worker->GetHitCountMap().cend(); ++it)
		fHitsMap[it->first] += it->second;

	for (auto it = worker->GetEnergyDepositedMap().cbegin();
		 it != worker->GetEnergyDepositedMap().cend(); ++it)
		fEnergyDepositedMap[it->first] += it->second;
}

void TsScoreHits::Output()
{
	std::vector<TsScoreHitsIndex> keys;
	for (auto it = fHitsMap.cbegin(); it != fHitsMap.cend(); ++it)
		keys.push_back(it->first);

	for (auto const& key : keys)
	{
		fParticleName = key.ParticleName;
		fVolumeName = key.VolumeName;
		fVolumeCopyNumber = key.VolumeCopyNumber;
		fTime = key.Time;
		fHits = fHitsMap[key];
		fEnergyDeposited = key.IsMolecule ? 0 : fEnergyDepositedMap[key];
		fNtuple->Fill();
	}

	fNtuple->Write();

	// report additional statistics to stdout
	G4cout << G4endl;
	G4cout << "Scorer: " << GetNameWithSplitId() << G4endl;
	if (fNtuple->HasHeaderFile())
		G4cout << "Header   has been written to file: " << fNtuple->GetHeaderFileName() << G4endl;
	G4cout << "Contents has been written to file: " << fNtuple->GetDataFileName() << G4endl;
	if (fPm->ParameterExists(GetFullParmName("Surface")))
		G4cout << "Scored on surface: " << fComponent->GetName() << "/" << GetSurfaceName() << G4endl;
	else
		G4cout << "Scored in component: " << fComponent->GetName() << G4endl;

	UpdateFileNameForUpcomingRun();
}

void TsScoreHits::Clear()
{
	fParticleName = "";
	fVolumeName = "";
	fVolumeCopyNumber = 0;
	// G4float fKineticEnergy;
	fEnergyDeposited = 0;
	fHits = 0;
	fScoredHistories = 0;
	fHitsMap.clear();
	fEnergyDepositedMap.clear();
}
