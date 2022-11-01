#ifndef TsFilterByComponentAndChildren_hh
#define TsFilterByComponentAndChildren_hh

#include "TsVFilter.hh"

class G4VPhysicalVolume;

class TsFilterByComponentAndChildren : public TsVFilter
{
public:
	TsFilterByComponentAndChildren(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
								   TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
	virtual ~TsFilterByComponentAndChildren() = default;

	G4bool Accept(const G4Step*) const override;
	G4bool AcceptTrack(const G4Track*) const override;

	void ResolveParameters() override;
	void CacheGeometryPointers() override;

private:
	G4bool fIncludeChildren;
	G4bool fNamesAreVolumes;

	std::vector<G4String> fNames;
	std::vector<G4VPhysicalVolume*> fVolumes;
};

#endif