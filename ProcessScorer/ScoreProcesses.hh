#ifndef ScoreProcesses_hh
#define ScoreProcesses_hh

#include "TsVNtupleScorer.hh"

class ScoreProcesses : public TsVNtupleScorer
{
public:
	ScoreProcesses(TsParameterManager *pM, TsMaterialManager *mM, TsGeometryManager *gM, TsScoringManager *scM, TsExtensionManager *eM,
				   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;
	void UserHookForEndOfTrack(const G4Track *) override;

protected:
	// Output variables
	G4int fParentId;
	G4double fEnergy;
	G4float fWeight;
	G4int fParticleType;
	G4String fOriginProcessName;
};
#endif
