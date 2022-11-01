#ifndef TsVFilterIRT_hh
#define TsVFilterIRT_hh

#include "TsVFilter.hh"

class TsParameterManager;
class TsGeometryManager;
class TsVScorer;

class TsVFilterIRT //: public TsVFilter
{
public:
	TsVFilterIRT(G4String, TsParameterManager* pM, TsGeometryManager* gM, TsVScorer*);
	// TsVFilterIRT(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
	// 			 TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
	virtual ~TsVFilterIRT() = default;

	virtual G4bool Accept(const G4ThreeVector&) const = 0;
	virtual G4bool Accept(const G4Step*) const = 0;
	// virtual G4bool AcceptTrack(const G4Track*) const = 0;

public:
	inline G4String GetName() const { return fName; }
	inline G4String GetFullName() const { return fFullName; }
	G4String GetFullParmName(const char*) const;

	inline const TsVScorer* GetScorer() const { return fScorer; }

protected:
	TsParameterManager* fPm;
	TsGeometryManager* fGm;

private:
	G4String fName, fFullName;
	TsVScorer* fScorer;
};

#endif