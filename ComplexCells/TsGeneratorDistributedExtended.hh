//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsGeneratorDistributed_hh
#define TsGeneratorDistributed_hh

#include "TsVGenerator.hh"

class TsGeneratorDistributedExtended : public TsVGenerator
{
public:
	TsGeneratorDistributedExtended(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName);
	~TsGeneratorDistributedExtended() = default;

	void GeneratePrimaries(G4Event*);
};

#endif
