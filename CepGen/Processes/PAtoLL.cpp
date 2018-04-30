#include "CepGen/Processes/PAtoLL.h"
#include "CepGen/Physics/PDG.h"

using namespace CepGen::Process;

PAtoLL::PAtoLL() : GenericKTProcess( "patoll", "ɣɣ → l⁺l¯", { { PDG::Photon, PDG::Photon } }, { PDG::Muon, PDG::Muon } )
{
}

void
PAtoLL::preparePhaseSpace()
{
}

double
PAtoLL::computeKTFactorisedMatrixElement()
{
}

void
PAtoLL::fillCentralParticlesKinematics()
{
}
