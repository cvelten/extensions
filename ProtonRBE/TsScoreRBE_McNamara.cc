// Scorer for RBE_McNamara
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

#include "TsScoreRBE_McNamara.hh"

#include "TsParameterManager.hh"

TsScoreRBE_McNamara::TsScoreRBE_McNamara(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsScoreRBE_MinMax(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	;
}

TsScoreRBE_McNamara::~TsScoreRBE_McNamara()
{
	;
}

TsVModelBiologicalEffect* TsScoreRBE_McNamara::ConstructModel(G4String cellLine)
{
	return new TsModelRBE_McNamara(cellLine, fPm, fOutputQuantity);
}

TsModelRBE_McNamara::TsModelRBE_McNamara(const G4String& cellLine, TsParameterManager* pM, const G4String& outputQuantity)
	: TsModelRBE_MinMax()
{
	fRBEmaxOffset = 0.99064;
	fRBEmaxSlope = 0.35605 * gray / (keV / um);
	fRBEminOffset = 1.1012;
	fRBEminSlope = -0.0038703 * 1 / (keV / um);

	G4String name = "Sc/" + cellLine + "/AlphaBetaRatiox";
	fAlphaBetax = pM->GetDoubleParameter(name, "Dose");

	if (outputQuantity == "alpha" || outputQuantity == "survivalfraction") {
		name = "Sc/" + cellLine + "/Alphax";
		fAlphax = pM->GetDoubleParameter(name, "perDose");
	}

	if (outputQuantity == "beta" || outputQuantity == "survivalfraction") {
		name = "Sc/" + cellLine + "/Betax";
		fBetax = pM->GetDoubleParameter(name, "perDoseSquare");
	}
}
