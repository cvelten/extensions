//
// ********************************************************************
// *                                                                  *
// * This  code implementation is an extension to developments by the *
// * TOPAS collaboration.                                             *
// * This extension  is an freely  available in  accordance  with the *
// * freeBSD license:                                                 *
// * Copyright (c) <2015>, <Harald Paganetti>                         *
// * All rights reserved.                                             *
// * Redistribution    and   use in   source and   binary    forms,   *
// * with or without modification, are permitted provided that the    *
// * following conditions are met:                                    *
// *                                                                  *
// *                                                                  *
// * 1. Redistributions of source code must retain the above          *
// * copyright notice, this                                           *
// * list of conditions and the following disclaimer.                 *
// * 2. Redistributions in binary form must reproduce the above       *
// * copyright notice, this list of conditions and the following      *
// * disclaimer in the documentation and/or other materials provided  *
// * with the distribution.                                           *
// *                                                                  *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
// * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
// * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
// * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
// * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR             *
// * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT *
// * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF *
// * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED  *
// * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      *
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING   *
// * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF   *
// * THE POSSIBILITY OF SUCH DAMAGE.                                  *
// *                                                                  *
// * The views and conclusions contained in the software and          *
// * documentation are those of the authors and should not be         *
// * interpreted as representing official policies, either expressed  *
// * or implied, of the FreeBSD Project.                              *
// *                                                                  *
// * Contacts: Jan Schuemann, jschuemann@mgh.harvard.edu              *
// *           Harald Paganetti, hpaganetti@mgh.harvard.edu           *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreDoseSqrtBeta_Tabulated_hh
#define TsScoreDoseSqrtBeta_Tabulated_hh

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "TsVScoreBiologicalEffect.hh"

class TsScoreDoseSqrtBeta_Tabulated : public TsVScoreBiologicalEffect
{
public:
	TsScoreDoseSqrtBeta_Tabulated(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
								  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer = false);
	virtual ~TsScoreDoseSqrtBeta_Tabulated();

	G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
	TsVModelBiologicalEffect* ConstructModel(G4String cellLine);
};

class TsModelBeta_Tabulated : public TsVModelBiologicalEffect
{
public:
	TsModelBeta_Tabulated(const G4String& cellLine, const G4String& modelName, TsParameterManager* pM);

	G4double InterpolateBeta(G4int particleZ, G4double kineticEnergyPerNucleon);

private:
	G4bool fUseReferenceBeta;
	G4double fBetax;
	G4int fNumberOfEnergyBins;
	G4int fNumberOfParticleNames;
	std::map<G4int, G4double*> fBeta;	// key: particle Z, value: beta binned by energy
	G4double* fKineticEnergyPerNucleon; // Energy bins for alpha and beta
};

#endif
