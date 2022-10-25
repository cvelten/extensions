// Scorer for MoleculeTracker
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

#include "TsScoreMoleculeTracker.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Scheduler.hh"
#include "G4VProcess.hh"

#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

#include "G4Molecule.hh"

TsScoreMoleculeTracker::TsScoreMoleculeTracker(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
	  fReportTimeSteps(true)
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
	fTimesToRecord.push_back(fTimeCut);

	fNtuple->RegisterColumnS(&fParticleName, "Particle Name");
	fNtuple->RegisterColumnS(&fVolumeName, "Volume Name");
	fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume Copy Number");
	fNtuple->RegisterColumnF(&fTime, "Time", "ps");
	fNtuple->RegisterColumnI(&fHits, "Hits");
}

G4bool TsScoreMoleculeTracker::ProcessHits(G4Step*, G4TouchableHistory*)
{
	if (!fIsActive)
		fSkippedWhileInactive++;
	return false;
}

void TsScoreMoleculeTracker::UserHookForPreTimeStepAction()
{
	// if (fReportTimeSteps)
	// {
	// 	G4cout << "TsScoreMoleculeTracker::UserHookForPreTimeStepAction() @ "
	// 		   << G4Scheduler::Instance()->GetTimeStep() / ps << " ps"
	// 		   << G4endl;
	// }
}
void TsScoreMoleculeTracker::UserHookForPostTimeStepAction()
{
	// if (fReportTimeSteps)
	// {
	// 	G4cout << "TsScoreMoleculeTracker::UserHookForPostTimeStepAction() @ "
	// 		   << G4Scheduler::Instance()->GetTimeStep() / ps << " ps"
	// 		   << G4endl;
	// }
}

void TsScoreMoleculeTracker::UserHookForChemicalStep(const G4Step* aStep)
{
	auto aTrack = aStep->GetTrack();

	if (aTrack->GetTrackID() > -1) // this is a physical track
		return;

	auto globalTime = aStep->GetPreStepPoint()->GetGlobalTime();
	auto itNextTime = fNextTimeForTrack.emplace(aTrack->GetTrackID(), fTimesToRecord.front()).first;

	if (globalTime >= itNextTime->second && itNextTime->second > 0)
	{
		if (GetFilter()->Accept(aStep))
		{
			G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

			TsScoreMoleculeTrackerIndex idx;
			idx.IsMolecule = true;
			idx.ParticleName = GetMolecule(aTrack)->GetName();
			idx.VolumeName = touchable->GetVolume()->GetName();
			idx.VolumeCopyNumber = touchable->GetVolume()->GetCopyNo();
			idx.Time = fNextTimeForTrack[aTrack->GetTrackID()];

			++fHitsMap[idx];

			if (fReportTimeSteps)
			{
				G4cout << "TsScoreMoleculeTracker::UserHookForChemicalStep(const G4Step* aStep): "
					   << "TrackID=" << aTrack->GetTrackID() << " @ "
					   << idx.Time / ps << " ps: "
					   << idx.ParticleName << " in " << idx.VolumeName
					   << G4endl;
			}
		}

		auto itNextTime = std::upper_bound(fTimesToRecord.begin(), fTimesToRecord.end(), globalTime);
		if (itNextTime == fTimesToRecord.end())
			fNextTimeForTrack[aTrack->GetTrackID()] = 0; // do not track anymore
		else
			fNextTimeForTrack[aTrack->GetTrackID()] = *itNextTime;
	}

	if (false)
	{
		if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
		{
			G4cout << "Step from " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
				   << " into " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
				   << G4endl;
		}
		else
		{
			G4cout << "Step within " << aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()
				   << G4endl;
		}
		G4cout << "Particle Name: "
			   << (aTrack->GetTrackID() < 0 ? GetMolecule(aTrack)->GetName() : aTrack->GetParticleDefinition()->GetParticleName())
			   << G4endl;
		G4cout << "Time: " << aStep->GetPreStepPoint()->GetGlobalTime() / ps << " ps"
			   << G4endl;
	}
}

void TsScoreMoleculeTracker::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

	TsScoreMoleculeTracker* worker = dynamic_cast<TsScoreMoleculeTracker*>(workerScorer);

	for (auto it = worker->GetHitCountMap().cbegin();
		 it != worker->GetHitCountMap().cend(); ++it)
		fHitsMap[it->first] += it->second;
}

void TsScoreMoleculeTracker::Output()
{
	std::vector<TsScoreMoleculeTrackerIndex> keys;
	for (auto it = fHitsMap.cbegin(); it != fHitsMap.cend(); ++it)
		keys.push_back(it->first);

	for (auto const& key : keys)
	{
		fParticleName = key.ParticleName;
		fVolumeName = key.VolumeName;
		fVolumeCopyNumber = key.VolumeCopyNumber;
		fTime = key.Time;
		fHits = fHitsMap[key];
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

void TsScoreMoleculeTracker::Clear()
{
	fParticleName = "";
	fVolumeName = "";
	fVolumeCopyNumber = 0;
	fHits = 0;
	fScoredHistories = 0;
	fHitsMap.clear();
}