// Extra Class for OnlyIncludeIfWithinComponentOrItsChildren

#include "TsFilterByPositionAndComponents.hh"

#include "TsGeometryManager.hh"
#include "TsVGeometryComponent.hh"

#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

// TsFilterByPositionAndComponents::TsFilterByPositionAndComponents(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
// 																 TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter)
//	: TsVFilterIRT(name, pM, mM, gM, fM, generator, scorer, parentFilter),
TsFilterByPositionAndComponents::TsFilterByPositionAndComponents(G4String name, TsParameterManager* pM, TsGeometryManager* gM, TsVScorer* scorer, TsVFilterIRT* parentFilter)
	: TsVFilterIRT(name, pM, gM, scorer), fParentFilter(parentFilter), fNamesAreVolumes(false)
{
	ResolveParameters();

	pM->SetNeedsSteppingAction();

	fNavigator = new G4Navigator();
	fNavigator->SetWorldVolume(fGm->GetPhysicalVolume("World"));
	// if (scorer->GetComponent() == nullptr || scorer->GetComponent()->GetEnvelopePhys() == nullptr) {
	// 	G4cerr << "Navigator set to WORLD" << G4endl;
	// 	fNavigator->SetWorldVolume(fGm->GetPhysicalVolume("World"));
	// }
	// else {
	// 	G4cerr << "Navigator set to " << scorer->GetComponent()->GetEnvelopePhys()->GetName() << G4endl;
	// 	fNavigator->SetWorldVolume(scorer->GetComponent()->GetEnvelopePhys());
	// }
}

TsFilterByPositionAndComponents::~TsFilterByPositionAndComponents()
{
	delete fNavigator;
}

void TsFilterByPositionAndComponents::ResolveParameters()
{
	G4String parName = GetName() + "/NamesAreVolumes";
	G4cerr << "Expecting: " << GetFullParmName(parName) << G4endl;
	if (fPm->ParameterExists(GetFullParmName(parName)))
		fNamesAreVolumes = fPm->GetBooleanParameter(GetFullParmName(parName));

	G4cerr << "Expecting: " << GetFullParmName(GetName()) << G4endl;
	auto n = fPm->GetVectorLength(GetFullParmName(GetName()));
	G4String* names = fPm->GetStringVector(GetFullParmName(GetName()));
	fNames = std::vector<G4String>(names, names + n);
	delete[] names;

	CacheGeometryPointers();
}

void TsFilterByPositionAndComponents::CacheGeometryPointers()
{
	fVolumes.clear();

	if (fNamesAreVolumes)
	{
		for (auto name : fNames)
		{
			G4VPhysicalVolume* volume = fGm->GetPhysicalVolume(name);
			if (volume == nullptr) {
				G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
				G4cerr << GetName() << " = " << name << " refers to an unknown Volume." << G4endl;
				fPm->AbortSession(1);
			}
			fVolumes.push_back(volume);
		}
	}
	else
	{
		for (auto name : fNames)
		{
			TsVGeometryComponent* component = fGm->GetComponent(name);
			if (component == nullptr) {
				G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
				G4cerr << GetName() << " = " << name << " refers to an unknown Component." << G4endl;
				fPm->AbortSession(1);
			}
			auto physicalVolumes = component->GetAllPhysicalVolumes();
			fVolumes.insert(fVolumes.end(), physicalVolumes.begin(), physicalVolumes.end());
		}
	}

	std::sort(fVolumes.begin(), fVolumes.end());

	G4cerr << "FilterVolumes:" << G4endl;
	for (auto volume : fVolumes)
		G4cerr << volume->GetName() << G4endl;
	G4cerr << G4endl;
}

G4bool TsFilterByPositionAndComponents::Accept(const G4ThreeVector& pos) const
{
	if (fParentFilter != nullptr && !fParentFilter->Accept(pos)) return false;

	auto pv = fNavigator->LocateGlobalPointAndSetup(pos);
	if (pv == nullptr)
	{
		G4cerr << pos << " is not associated with a PV" << G4endl;
		return false;
	}
	if (std::binary_search(fVolumes.cbegin(), fVolumes.cend(), pv))
	{
		G4cerr << pos << " found in volumes" << G4endl;
		return true;
	}
	G4cerr << pos << " NOT found in volumes" << G4endl;
	return false;
}

G4bool TsFilterByPositionAndComponents::Accept(const G4Step* aStep) const
{
	if (fParentFilter != nullptr && !fParentFilter->Accept(aStep)) return false;

	auto touchable = aStep->GetPreStepPoint()->GetTouchable();

	G4VPhysicalVolume* pv;
	for (auto depth = 0; depth - 1 < touchable->GetHistoryDepth(); ++depth)
	{
		pv = touchable->GetVolume(depth);
		if (pv == nullptr)
			continue;
		if (std::binary_search(fVolumes.begin(), fVolumes.end(), pv))
			return true;
	}
	return false;
}
