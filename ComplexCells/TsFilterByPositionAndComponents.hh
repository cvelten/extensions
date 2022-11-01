#ifndef TsFilterByPositionAndComponents_hh
#define TsFilterByPositionAndComponents_hh

#include "TsVFilterIRT.hh"

class G4Navigator;
class G4VPhysicalVolume;
class TsVScorer;

class TsFilterByPositionAndComponents : public TsVFilterIRT
{
public:
	TsFilterByPositionAndComponents(G4String name, TsParameterManager* pM, TsGeometryManager* gM, TsVScorer* scorer, TsVFilterIRT* parentFilter = nullptr);
	// TsFilterByPositionAndComponents(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
	// 								TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
	virtual ~TsFilterByPositionAndComponents();

	G4bool Accept(const G4ThreeVector&) const override;
	G4bool Accept(const G4Step*) const override;
	// G4bool AcceptTrack(const G4Track*) const override;

	void ResolveParameters();
	void CacheGeometryPointers();

private:
	TsVFilterIRT* fParentFilter;

	G4Navigator* fNavigator;
	G4bool fNamesAreVolumes;

	std::vector<G4String> fNames;
	std::vector<G4VPhysicalVolume*> fVolumes;
};

#endif