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

#ifndef TsScoreHits_hh
#define TsScoreHits_hh

#include "TsVNtupleScorer.hh"

struct TsScoreHitsIndex;

class TsScoreHits : public TsVNtupleScorer
{
public:
	TsScoreHits(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
				G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsScoreHits() = default;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
	void AbsorbResultsFromWorkerScorer(TsVScorer*) override;

protected:
	void Output() override;
	void Clear() override;

private:
	G4String fParticleName;
	G4String fVolumeName;
	G4int fVolumeCopyNumber;
	// G4float fKineticEnergy;
	G4float fEnergyDeposited;
	G4int fHits;
	G4float fTime;
	// G4float fGlobalTime;
	// G4int fTrackID;

	G4double fTimeCut;

private:
	std::vector<G4double> fTimesToRecord;
	std::map<G4int, G4double> fNextTimeForTrack;

	std::map<TsScoreHitsIndex, G4int> fHitsMap;
	std::map<TsScoreHitsIndex, G4double> fEnergyDepositedMap;

	G4bool fIncludeChemistry;
	G4bool fIncludePhysics;

public:
	const std::map<TsScoreHitsIndex, G4int>& GetHitCountMap() const { return fHitsMap; }
	const std::map<TsScoreHitsIndex, G4double>& GetEnergyDepositedMap() const { return fEnergyDepositedMap; }
};

struct TsScoreHitsIndex
{
	G4String ParticleName;
	G4String VolumeName;
	G4int VolumeCopyNumber;
	G4bool IsMolecule;
	G4double Time;

	TsScoreHitsIndex() : IsMolecule(false), Time(0) {}

	bool operator<(TsScoreHitsIndex const& other) const
	{
		return ParticleName < other.ParticleName ||
			   VolumeName < other.VolumeName ||
			   VolumeCopyNumber < other.VolumeCopyNumber ||
			   Time < other.Time;
	}
};

#endif
