// Filter for OnlyIncludeIfWithinComponentOrItsChildren

#include "TsFilterByComponentAndChildren.hh"

#include "TsVGeometryComponent.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

TsFilterByComponentAndChildren::TsFilterByComponentAndChildren(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
															   TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter)
	: TsVFilter(name, pM, mM, gM, fM, generator, scorer, parentFilter),
	  fIncludeChildren(), fNamesAreVolumes()
{
	ResolveParameters();

	pM->SetNeedsSteppingAction();
}

void TsFilterByComponentAndChildren::ResolveParameters()
{
	G4String parName;

	parName = GetName() + "/IncludeChildren";
	if (fPm->ParameterExists(GetFullParmName(parName)))
		fIncludeChildren = fPm->GetBooleanParameter(GetFullParmName(parName));
	else
		fIncludeChildren = true;

	parName = GetName() + "/NamesAreVolumes";
	if (fPm->ParameterExists(GetFullParmName(parName)))
		fNamesAreVolumes = fPm->GetBooleanParameter(GetFullParmName(parName));

	auto n = fPm->GetVectorLength(GetFullParmName(GetName()));
	G4String* names = fPm->GetStringVector(GetFullParmName(GetName()));
	fNames = std::vector<G4String>(names, names + n);
	delete[] names;

	CacheGeometryPointers();
}

void TsFilterByComponentAndChildren::CacheGeometryPointers()
{
	fVolumes.clear();

	if (fNamesAreVolumes)
	{
		for (auto name : fNames)
		{
			G4VPhysicalVolume* volume = GetPhysicalVolume(name);
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
			TsVGeometryComponent* component = GetComponent(name);
			if (component == nullptr) {
				G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
				G4cerr << GetName() << " = " << name << " refers to an unknown Component." << G4endl;
				fPm->AbortSession(1);
			}
			auto physicalVolumes = component->GetAllPhysicalVolumes(false);
			fVolumes.insert(fVolumes.end(), physicalVolumes.begin(), physicalVolumes.end());
		}
	}

	std::sort(fVolumes.begin(), fVolumes.end());
}

G4bool TsFilterByComponentAndChildren::Accept(const G4Step* aStep) const
{
	if (fParentFilter && !fParentFilter->Accept(aStep)) return false;

	const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();

	G4VPhysicalVolume* pv;
	size_t maxDepth = fIncludeChildren ? touchable->GetHistoryDepth() + 1 : 1;
	for (size_t depth = 0; depth < maxDepth; ++depth)
	{
		pv = touchable->GetVolume(depth);
		if (pv == nullptr)
			continue;
		if (std::binary_search(fVolumes.begin(), fVolumes.end(), pv))
			return true;
	}
	return false;
}

G4bool TsFilterByComponentAndChildren::AcceptTrack(const G4Track*) const
{
	G4cerr << "Topas is exiting due to a serious error in source setup." << G4endl;
	G4cerr << "Sources cannot be filtered by " << GetName() << G4endl;
	fPm->AbortSession(1);
	return false;
}
