#ifndef TsScoreMoleculeTracker_hh
#define TsScoreMoleculeTracker_hh

#include "TsVNtupleScorer.hh"

struct TsScoreMoleculeTrackerIndex;

class TsScoreMoleculeTracker : public TsVNtupleScorer
{
public:
	TsScoreMoleculeTracker(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreMoleculeTracker() override = default;

	void ResolveParameters() override;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
	void AbsorbResultsFromWorkerScorer(TsVScorer*) override;

	void UserHookForPreTimeStepAction() override;
	void UserHookForPostTimeStepAction() override;
	void UserHookForBeginOfChemicalTrack(const G4Track*) override;
	void UserHookForChemicalStep(const G4Step*) override;

protected:
	void Output() override;
	void Clear() override;

	bool TrackWithinVolumes(const G4Track*, const std::vector<G4VPhysicalVolume*>&) const;

private:
	G4String fParticleName;
	G4String fVolumeName;
	G4int fVolumeCopyNumber;
	G4int fHits;
	G4float fTime;

	G4double fTimeCut;

	std::vector<G4double> fTimesToRecord;

	std::vector<G4VPhysicalVolume*> fKillMoleculeInVolumes;
	std::vector<G4VPhysicalVolume*> fKillMoleculeNotInVolumes;

	std::map<G4int, G4double> fNextTimeForTrack;
	std::map<TsScoreMoleculeTrackerIndex, G4int> fHitsMap;

public:
	const std::map<TsScoreMoleculeTrackerIndex, G4int>& GetHitCountMap() const { return fHitsMap; }
};

struct TsScoreMoleculeTrackerIndex
{
	G4String ParticleName;
	G4String VolumeName;
	G4int VolumeCopyNumber;
	G4bool IsMolecule;
	G4double Time;

	TsScoreMoleculeTrackerIndex() : IsMolecule(false), Time(0) {}

	bool operator<(TsScoreMoleculeTrackerIndex const& other) const
	{
		return ParticleName < other.ParticleName ||
			   (ParticleName == other.ParticleName && VolumeName < other.VolumeName) ||
			   (ParticleName == other.ParticleName && VolumeName == other.VolumeName && VolumeCopyNumber < other.VolumeCopyNumber) ||
			   (ParticleName == other.ParticleName && VolumeName == other.VolumeName && VolumeCopyNumber == other.VolumeCopyNumber && Time < other.Time);
	}
};

#endif
