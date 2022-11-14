#ifndef TsTrackTerminator_hh
#define TsTrackTerminator_hh

#include "TsVScorer.hh"

class TsTrackTerminator : public TsVScorer
{
public:
	TsTrackTerminator(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsTrackTerminator() override = default;

	void ResolveParameters() override;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

	void UserHookForChemicalStep(const G4Step*) override;
	void UserHookForEndOfEvent() override;

	void RestoreResultsFromFile() override {}
	void AccumulateEvent() override;
	void AbsorbResultsFromWorkerScorer(TsVScorer*) override {}

protected:
	void Output() override {}
	void Clear() override {}

private:
	TsParameterManager* fPm;

	G4int fNbOfMoleculesToScavenge;
	std::vector<G4int> fScavengeTheseMolecules;
	std::vector<G4double> fScavengingCapacities;

	G4double fEnergyLossKill;
	G4double fEnergyLossAbort;
	G4double fEnergyLoss;
	G4double fMaximumTrackLength;
	G4double fTotalTrackLength;
};

#endif
