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

#ifndef TsMoleculeTracker_hh
#define TsMoleculeTracker_hh

#include "TsVNtupleScorer.hh"

struct TsMoleculeTrackerIndex;

class TsMoleculeTracker : public TsVNtupleScorer
{
public:
	TsMoleculeTracker(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsMoleculeTracker() = default;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
	void AbsorbResultsFromWorkerScorer(TsVScorer*) override;

	void UserHookForChemicalStep(const G4Step*) override;

protected:
	void
	Output() override;
	void Clear() override;

private:
	G4String fParticleName;
	G4String fVolumeName;
	G4int fVolumeCopyNumber;
	G4int fHits;
	G4float fTime;
	// G4float fGlobalTime;
	// G4int fTrackID;

	G4double fTimeCut;

private:
	std::vector<G4double> fTimesToRecord;
	std::map<G4int, G4double> fNextTimeForTrack;
	std::map<TsMoleculeTrackerIndex, G4int> fHitsMap;

public:
	const std::map<TsMoleculeTrackerIndex, G4int>& GetHitCountMap() const { return fHitsMap; }
};

struct TsMoleculeTrackerIndex
{
	G4String ParticleName;
	G4String VolumeName;
	G4int VolumeCopyNumber;
	G4bool IsMolecule;
	G4double Time;

	TsMoleculeTrackerIndex() : IsMolecule(false), Time(0) {}

	bool operator<(TsMoleculeTrackerIndex const& other) const
	{
		return ParticleName < other.ParticleName ||
			   VolumeName < other.VolumeName ||
			   VolumeCopyNumber < other.VolumeCopyNumber ||
			   Time < other.Time;
	}
};

#endif
