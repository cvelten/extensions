// Scorer for TupleExtended

#include "TsScoreTupleExtended.hh"

#include "G4Molecule.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"

TsScoreTupleExtended::TsScoreTupleExtended(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
	  fTimesToRecord(), fNextTimeForTrack(),
	  fIncludeKineticEnergy(false),
	  fIncludeEventID(false), fIncludeTrackID(false), fIncludeParentID(false), fIncludeStepNumber(false),
	  fIncludeParticleName(false), fIncludeProcessName(false), fIncludeVolumeName(false), fIncludeVolumeCopyNumber(false),
	  fIncludeGlobalTime(false), fIncludeEnergyDeposited(false), fIncludeVertex(false)
{
	SetUnit("");

	fNtuple->RegisterColumnI(&fMoleculeID, "MoleculeID or ParticlePDG");
	fNtuple->RegisterColumnF(&fPosX, "Position X", "um");
	fNtuple->RegisterColumnF(&fPosY, "Position Y", "um");
	fNtuple->RegisterColumnF(&fPosZ, "Position Z", "um");
	fNtuple->RegisterColumnF(&fTime, "Time", "ps");

	fTimeCut = 1.0 * us;
	if (fPm->ParameterExists(GetFullParmName("TimeCut")))
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"), "Time");

	if (fPm->ParameterExists(GetFullParmName("TimeToRecord")))
	{
		auto times = fPm->GetDoubleVector(GetFullParmName("TimeToRecord"), "Time");
		for (auto i = 0; i < fPm->GetVectorLength(GetFullParmName("TimeToRecord")); ++i)
		{
			if (times[i] > fTimeCut)
			{
				G4cerr << GetName() << " specified a Time to record greater than the Time cut (default: 1 us)" << G4endl;
				G4cerr << "at which chemical tracks will be killed" << G4endl;
			}
			fTimesToRecord.push_back(times[i]);
		}

		std::sort(fTimesToRecord.begin(), fTimesToRecord.end());
		// std::reverse(fTimesToRecord.begin(), fTimesToRecord.end());
	}
	else {
		fTimesToRecord.push_back(fTimeCut);
	}

	fIncludeChemistry = true;
	fIncludePhysics = true;
	if (fPm->ParameterExists(GetFullParmName("IncludeChemicalTrack")))
		fIncludeChemistry = fPm->GetBooleanParameter(GetFullParmName("IncludeChemicalTrack"));

	if (fPm->ParameterExists(GetFullParmName("IncludePhysicalTrack")))
		fIncludePhysics = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalTrack"));

	if (fPm->ParameterExists(GetFullParmName("IncludeKineticEnergy")))
		fIncludeKineticEnergy = fPm->GetBooleanParameter(GetFullParmName("IncludeKineticEnergy"));

	if (fPm->ParameterExists(GetFullParmName("IncludeEventID")))
		fIncludeEventID = fPm->GetBooleanParameter(GetFullParmName("IncludeEventID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeTrackID")))
		fIncludeTrackID = fPm->GetBooleanParameter(GetFullParmName("IncludeTrackID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeParentID")))
		fIncludeParentID = fPm->GetBooleanParameter(GetFullParmName("IncludeParentID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeStepNumber")))
		fIncludeStepNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeStepNumber"));

	if (fPm->ParameterExists(GetFullParmName("IncludeParticleName")))
		fIncludeParticleName = fPm->GetBooleanParameter(GetFullParmName("IncludeParticleName"));

	if (fPm->ParameterExists(GetFullParmName("IncludePhysicalProcessName")))
		fIncludeProcessName = fPm->GetBooleanParameter(GetFullParmName("IncludePhysicalProcessName"));

	if (fPm->ParameterExists(GetFullParmName("IncludeVolumeName")))
		fIncludeVolumeName = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeName"));

	if (fPm->ParameterExists(GetFullParmName("IncludeVolumeCopyNumber")))
		fIncludeVolumeCopyNumber = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeCopyNumber"));

	if (fPm->ParameterExists(GetFullParmName("IncludeGlobalTime")))
		fIncludeGlobalTime = fPm->GetBooleanParameter(GetFullParmName("IncludeGlobalTime"));

	if (fPm->ParameterExists(GetFullParmName("IncludeEnergyDeposited")))
		fIncludeEnergyDeposited = fPm->GetBooleanParameter(GetFullParmName("IncludeEnergyDeposited"));

	if (fPm->ParameterExists(GetFullParmName("IncludeVertexPosition")))
		fIncludeVertex = fPm->GetBooleanParameter(GetFullParmName("IncludeVertexPosition"));

	if (fIncludeEventID)
		fNtuple->RegisterColumnI(&fEvt, "EventID");

	if (fIncludeTrackID)
		fNtuple->RegisterColumnI(&fTrackID, "TrackID");

	if (fIncludeStepNumber)
		fNtuple->RegisterColumnI(&fStepNumber, "Step Number");

	if (fIncludeParticleName)
		fNtuple->RegisterColumnS(&fParticleName, "Particle Name");

	if (fIncludeProcessName)
		fNtuple->RegisterColumnS(&fProcessName, "Process Name");

	if (fIncludeVolumeName)
		fNtuple->RegisterColumnS(&fVolumeName, "Volume Name");

	if (fIncludeVolumeCopyNumber)
		fNtuple->RegisterColumnI(&fVolumeCopyNumber, "Volume copy Number");

	if (fIncludeParentID) {
		fNtuple->RegisterColumnI(&fParentAID, "ParentA ID");
		fNtuple->RegisterColumnI(&fParentBID, "ParentB ID");
	}

	if (fIncludeVertex) {
		fNtuple->RegisterColumnF(&fVertexPositionX, "Vertex Position X", "um");
		fNtuple->RegisterColumnF(&fVertexPositionY, "Vertex Position Y", "um");
		fNtuple->RegisterColumnF(&fVertexPositionZ, "Vertex Position Z", "um");
	}

	if (fIncludeGlobalTime)
		fNtuple->RegisterColumnF(&fGlobalTime, "Global Time", "ps");

	if (fIncludeEnergyDeposited)
		fNtuple->RegisterColumnF(&fEnergyDeposited, "Energy Deposited", "keV");

	if (fIncludeKineticEnergy)
		fNtuple->RegisterColumnF(&fKineticEnergy, "Kinetic Energy", "keV");
}

G4bool TsScoreTupleExtended::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	G4Track* aTrack = aStep->GetTrack();

	if (fIncludeChemistry && aTrack->GetTrackID() < 0) {
		fGlobalTime = aStep->GetPreStepPoint()->GetGlobalTime();
		fTime = fNextTimeForTrack.emplace(aTrack->GetTrackID(), fTimesToRecord.front()).first->second;

		if (fTime > 0) {
			if (fGlobalTime >= fTime)
			{
				G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

				fVolumeName = touchable->GetVolume()->GetName();
				fVolumeCopyNumber = touchable->GetVolume()->GetCopyNo();

				fEvt = GetEventID();
				fKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
				G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
				fPosX = pos.x();
				fPosY = pos.y();
				fPosZ = pos.z();
				if (fIncludeVertex) {
					G4ThreeVector vpos = aTrack->GetVertexPosition();
					fVertexPositionX = vpos.x();
					fVertexPositionY = vpos.y();
					fVertexPositionZ = vpos.z();
				}

				fParentAID = -1;
				fParentBID = -1;
				fTrackID = aTrack->GetTrackID();

				fParticleName = GetMolecule(aTrack)->GetName();
				fMoleculeID = GetMolecule(aTrack)->GetMoleculeID();

				GetMolecule(aTrack)->GetParentID(fParentAID, fParentBID);
				fEnergyDeposited = 0.0;
				fProcessName = "none";
				fStepNumber = aTrack->GetCurrentStepNumber();

				fNtuple->Fill();

				return true;
			}

			// Set next tracking time for track
			auto itNextTrackingTime = std::upper_bound(fTimesToRecord.cbegin(), fTimesToRecord.cend(), fGlobalTime);
			if (itNextTrackingTime == fTimesToRecord.cend())
				fNextTimeForTrack[aTrack->GetTrackID()] = 0; // do not track anymore
			else
				fNextTimeForTrack[aTrack->GetTrackID()] = *itNextTrackingTime;
		}
	}
	else if (fIncludePhysics && aStep->GetTotalEnergyDeposit() > 0) {
		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		fVolumeName = touchable->GetVolume()->GetName();
		fVolumeCopyNumber = touchable->GetVolume()->GetCopyNo();

		fEvt = GetEventID();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		fPosX = pos.x();
		fPosY = pos.y();
		fPosZ = pos.z();
		if (fIncludeVertex) {
			G4ThreeVector vpos = aTrack->GetVertexPosition();
			fVertexPositionX = vpos.x();
			fVertexPositionY = vpos.y();
			fVertexPositionZ = vpos.z();
		}
		fParentAID = aTrack->GetParentID();
		fParentBID = -1;
		fTrackID = aTrack->GetTrackID();
		fEnergyDeposited = aStep->GetTotalEnergyDeposit();
		fStepNumber = aTrack->GetCurrentStepNumber();
		fKineticEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
		fMoleculeID = aTrack->GetParticleDefinition()->GetPDGEncoding();
		fParticleName = aTrack->GetParticleDefinition()->GetParticleName();
		fProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

		fNtuple->Fill();

		return true;
	}

	return false;
}
