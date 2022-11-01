// Scorer for MoleculeTracker

#include "TsScoreMoleculeTracker.hh"

#include "TsGeometryManager.hh"

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
	  fTimeCut(1 * us), fTimesToRecord(),
	  fKillMoleculeInVolumes(), fKillMoleculeNotInVolumes(),
	  fNextTimeForTrack(), fHitsMap()
{
	SetUnit("");

	ResolveParameters();

	fNtuple->RegisterColumnS(&fParticleName, "Particle Name");
	fNtuple->RegisterColumnS(&fVolumeName, "Volume Name");
	fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume Copy Number");
	fNtuple->RegisterColumnF(&fTime, "Time", "ps");
	fNtuple->RegisterColumnI(&fHits, "Hits");
}

void TsScoreMoleculeTracker::ResolveParameters()
{
	if (fPm->ParameterExists(GetFullParmName("TimeCut")))
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"), "Time");

	if (fPm->ParameterExists(GetFullParmName("TimesToRecord")))
	{
		auto times = fPm->GetDoubleVector(GetFullParmName("TimesToRecord"), "Time");
		for (auto i = 0; i < fPm->GetVectorLength(GetFullParmName("TimesToRecord")); ++i)
		{
			if (times[i] > fTimeCut)
			{
				G4cout << GetName() << " specified a time to record greater than the time cut (default: 1 us)" << G4endl;
				G4cout << "at which chemical tracks will be killed" << G4endl;
			}
			fTimesToRecord.push_back(times[i]);
		}

		std::sort(fTimesToRecord.begin(), fTimesToRecord.end());
		// std::reverse(fTimesToRecord.begin(), fTimesToRecord.end());
	}
	fTimesToRecord.push_back(fTimeCut);

	if (fPm->ParameterExists(GetFullParmName("KillMoleculeIfInComponents")))
	{
		auto n = fPm->GetVectorLength(GetFullParmName("KillMoleculeIfInComponents"));
		auto killInVolumes = fPm->GetStringVector(GetFullParmName("KillMoleculeIfInComponents"));
		for (auto i = 0; i < n; ++i)
		{
			TsVGeometryComponent* component = fGm->GetComponent(killInVolumes[i]);
			if (component == nullptr) {
				G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
				G4cerr << GetName() << "/KillMoleculeIfInComponents = " << killInVolumes[i] << " refers to an unknown Component." << G4endl;
				fPm->AbortSession(1);
			}
			auto physicalVolumes = component->GetAllPhysicalVolumes();
			fKillMoleculeInVolumes.insert(fKillMoleculeInVolumes.end(), physicalVolumes.begin(), physicalVolumes.end());
		}
		std::sort(fKillMoleculeInVolumes.begin(), fKillMoleculeInVolumes.end());
	}

	if (fPm->ParameterExists(GetFullParmName("KillMoleculeIfNotInComponents")))
	{
		auto n = fPm->GetVectorLength(GetFullParmName("KillMoleculeIfNotInComponents"));
		auto killNotInVolumes = fPm->GetStringVector(GetFullParmName("KillMoleculeIfNotInComponents"));
		for (auto i = 0; i < n; ++i)
		{
			TsVGeometryComponent* component = fGm->GetComponent(killNotInVolumes[i]);
			if (component == nullptr) {
				G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
				G4cerr << GetName() << "/KillMoleculeIfNotInComponents = " << killNotInVolumes[i] << " refers to an unknown Component." << G4endl;
				fPm->AbortSession(1);
			}
			auto physicalVolumes = component->GetAllPhysicalVolumes();
			fKillMoleculeNotInVolumes.insert(fKillMoleculeNotInVolumes.end(), physicalVolumes.begin(), physicalVolumes.end());
		}
		std::sort(fKillMoleculeNotInVolumes.begin(), fKillMoleculeNotInVolumes.end());
	}
}

G4bool TsScoreMoleculeTracker::ProcessHits(G4Step*, G4TouchableHistory*)
{
	if (!fIsActive)
		fSkippedWhileInactive++;
	return false;
}

void TsScoreMoleculeTracker::UserHookForPreTimeStepAction()
{
	if (fVerbosity > 1 && fTrackingVerbosity > 1)
	{
		G4cout << "TsScoreMoleculeTracker::UserHookForPreTimeStepAction() @ "
			   << G4Scheduler::Instance()->GetTimeStep() / ps << " ps "
			   << "with " << G4Scheduler::Instance()->GetNTracks() << " tracks in the stack."
			   << G4endl;
	}
}
void TsScoreMoleculeTracker::UserHookForPostTimeStepAction()
{
	if (fVerbosity > 1 && fTrackingVerbosity > 1)
	{
		G4cout << "TsScoreMoleculeTracker::UserHookForPostTimeStepAction() @ "
			   << G4Scheduler::Instance()->GetTimeStep() / ps << " ps "
			   << "with " << G4Scheduler::Instance()->GetNTracks() << " tracks in the stack."
			   << G4endl;
	}
}

void TsScoreMoleculeTracker::UserHookForBeginOfChemicalTrack(const G4Track* aTrack)
{
	if ((fKillMoleculeInVolumes.size() > 0 && TrackWithinVolumes(aTrack, fKillMoleculeInVolumes)) ||
		(fKillMoleculeNotInVolumes.size() > 0 && !TrackWithinVolumes(aTrack, fKillMoleculeNotInVolumes)))
	{
		if (fVerbosity > 2 || fTrackingVerbosity > 2)
		{
			G4cout << "TrackID = " << aTrack->GetTrackID() << G4endl;
			if (fKillMoleculeNotInVolumes.size() > 0 && !TrackWithinVolumes(aTrack, fKillMoleculeNotInVolumes)) {
				G4cout << "Not present in ";
				for (auto v : fKillMoleculeNotInVolumes) {
					G4cout << v->GetName() << " ";
				}
				G4cout << G4endl;
			}
			if (fKillMoleculeInVolumes.size() > 0 && TrackWithinVolumes(aTrack, fKillMoleculeInVolumes)) {
				G4cout << "Present in ";
				for (auto v : fKillMoleculeNotInVolumes) {
					G4cout << v->GetName() << " ";
				}
				G4cout << G4endl;
			}
		}

		if (fVerbosity > 0 || fTrackingVerbosity > 0)
			G4cout << "TsScoreMoleculeTracker::UserHookForBeginOfChemicalTrack(const G4Track*): "
				   << "Killed track " << aTrack->GetTrackID()
				   << G4endl;
		const_cast<G4Track*>(aTrack)->SetTrackStatus(fStopAndKill);
	}
}

void TsScoreMoleculeTracker::UserHookForChemicalStep(const G4Step* aStep)
{
	auto aTrack = aStep->GetTrack();

	if (aTrack->GetTrackID() > -1) // this is a physical track
		return;

	auto globalTime = aStep->GetPreStepPoint()->GetGlobalTime();
	auto itNextTime = fNextTimeForTrack.emplace(aTrack->GetTrackID(), fTimesToRecord.front()).first->second;
	if (itNextTime > 0) {
		if (globalTime >= itNextTime)
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

				if (fTrackingVerbosity > 2)
				{
					G4cout << "TsScoreMoleculeTracker::UserHookForChemicalStep(const G4Step*): "
						   << "TrackID=" << aTrack->GetTrackID() << " @ "
						   << idx.Time / ps << " ps: " << idx.ParticleName << " in " << idx.VolumeName
						   << G4endl;
				}
			}

			// Set next tracking time for track
			auto itNextTrackingTime = std::upper_bound(fTimesToRecord.begin(), fTimesToRecord.end(), globalTime);
			if (itNextTrackingTime == fTimesToRecord.end())
				fNextTimeForTrack[aTrack->GetTrackID()] = 0; // do not track anymore
			else
				fNextTimeForTrack[aTrack->GetTrackID()] = *itNextTrackingTime;
		}
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
	if (fVerbosity > 0) {
		G4cout << G4endl;
		G4cout << "Scorer: " << GetNameWithSplitId() << G4endl;
		if (fNtuple->HasHeaderFile())
			G4cout << "Header   has been written to file: " << fNtuple->GetHeaderFileName() << G4endl;
		G4cout << "Contents has been written to file: " << fNtuple->GetDataFileName() << G4endl;
		if (fPm->ParameterExists(GetFullParmName("Surface")))
			G4cout << "Scored on surface: " << fComponent->GetName() << "/" << GetSurfaceName() << G4endl;
		else
			G4cout << "Scored in component: " << fComponent->GetName() << G4endl;
	}

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

bool TsScoreMoleculeTracker::TrackWithinVolumes(const G4Track* aTrack, const std::vector<G4VPhysicalVolume*>& volumes) const
{
	if (volumes.size() == 0) return false;

	G4VPhysicalVolume* pv = aTrack->GetTouchable()->GetVolume();
	if (pv == nullptr)
		return false;

	return std::binary_search(volumes.cbegin(), volumes.cend(), pv);
}
