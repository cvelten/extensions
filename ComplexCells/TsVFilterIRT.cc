// Extra Class for TsEmDNAChemistry

#include "TsVFilterIRT.hh"

#include "TsGeometryManager.hh"
#include "TsParameterManager.hh"
#include "TsVScorer.hh"

/*TsVFilterIRT::TsVFilterIRT(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
						   TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter)
	: TsVFilter(name, pM, mM, gM, fM, generator, scorer, parentFilter)
{}*/

TsVFilterIRT::TsVFilterIRT(G4String name, TsParameterManager* pM, TsGeometryManager* gM, TsVScorer* scorer)
	: fPm(pM), fGm(gM), fName(name), fScorer(scorer)
{
	fFullName = scorer->GetName() + "/" + name;
}

G4String TsVFilterIRT::GetFullParmName(const char* parm) const
{
	return fScorer->GetFullParmName(parm);
}