#ifndef TsScoreTupleExtended_hh
#define TsScoreTupleExtended_hh

#include "TsVNtupleScorer.hh"

class TsScoreTupleExtended : public TsVNtupleScorer
{
public:
	TsScoreTupleExtended(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsScoreTupleExtended() = default;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
	G4int fEvt;
	G4String fParticleName;
	G4String fProcessName;
	G4String fVolumeName;
	G4int fVolumeCopyNumber;
	G4float fPosX;
	G4float fPosY;
	G4float fPosZ;
	G4float fKineticEnergy;
	G4float fEnergyDeposited;
	G4float fTime;
	G4float fGlobalTime;
	G4int fTrackID;
	G4int fParentAID;
	G4int fParentBID;
	G4int fMoleculeID;
	G4int fStepNumber;
	G4float fVertexPositionX;
	G4float fVertexPositionY;
	G4float fVertexPositionZ;

	G4double fTimeCut;

private:
	std::vector<G4double> fTimesToRecord;

	std::map<G4int, G4double> fNextTimeForTrack;

	G4bool fIncludeChemistry;
	G4bool fIncludePhysics;
	G4bool fIncludeKineticEnergy;
	G4bool fIncludeEventID;
	G4bool fIncludeTrackID;
	G4bool fIncludeParentID;
	G4bool fIncludeStepNumber;
	G4bool fIncludeParticleName;
	G4bool fIncludeProcessName;
	G4bool fIncludeVolumeName;
	G4bool fIncludeVolumeCopyNumber;
	G4bool fIncludeGlobalTime;
	G4bool fIncludeEnergyDeposited;
	G4bool fIncludeVertex;
};
#endif
