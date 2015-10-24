#include "Particle.h"

Particle::Particle() :
  id(-1), pdgId((ParticleCode)0), charge(999.), name(""), role(-1),
  helicity(0.),
  status(0), fMass(-1.),
  fIsPrimary(true)
{
  for (int i=0; i<4; i++) fP4[i] = 0.;
}

Particle::Particle(int role_, ParticleCode pdgId_) :
  id(-1), pdgId(pdgId_), charge(999.), name(""), role(role_),
  status(0), fMass(-1.),
  fIsPrimary(true)
{
  for (int i=0; i<4; i++) fP4[i] = 0.;
  if (pdgId!=0) M(-1.);
}

Particle::~Particle()
{
  //DebugInsideLoop("Destructor called");
}

Particle&
Particle::operator=(const Particle &part_)
{
  this->pdgId = part_.pdgId;
  this->role = part_.role;
  if (this->id==-1) {
    this->id = part_.id;
  }
  this->P(part_.Px(), part_.Py(), part_.Pz(), part_.E());
  this->M(part_.fMass);

  return *this;
}

Particle&
Particle::operator+=(const Particle &part_)
{
  pdgId = (part_.pdgId==pdgId) ? pdgId : (ParticleCode)(-1);
  role = (part_.role==role) ? role : -1;
  for (int i=0; i<4; i++) fP4[i]+= part_.fP4[i];
  M(-1.);

  return *this;
}

Particle&
Particle::operator-(const Particle &part_)
{
  pdgId = (part_.pdgId==pdgId) ? pdgId : (ParticleCode)(-1);
  role = (part_.role==role) ? role : -1;
  for (int i=0; i<4; i++) fP4[i]-= part_.fP4[i];
  M(-1.);

  return *this;
}

bool
Particle::Valid()
{
  if (pdgId==0) return false;
  if (P()==0. and M()==0.) return false;
  return true;
}

std::string
Particle::GetLHEline(bool revert_)
{
  std::stringstream line;

  if (revert_) {
    fP4[2] = -fP4[2];
  }

  line << pdgId << "\t";
  line << "1 1 2 0 0" << "\t";
  line << Px() << "\t";
  line << Py() << "\t";
  line << Pz() << "\t";
  line << E() << "\t";
  line << M() << "\t";
  line << "0." << "\t";
  line << "0."; //FIXME iz!!!
  return line.str();
}

bool
Particle::M(double m_)
{
  double mass;
  if (m_>=0.) fMass = m_;
  else if (pdgId!=0) {
    mass = GetMassFromPDGId(pdgId);
    if (mass<0.) return false;
    if (fP4[3]<0.) {
      fMass = mass;
      fP4[3] = std::pow(P(), 2)+M2();
      return true;
    }
    if (std::pow(E(), 2)-std::pow(P(), 2)!=std::pow(mass, 2)) mass = std::sqrt(std::pow(E(), 2)-std::pow(P(), 2));
    fMass = mass;
    return true;
  }
  else if (E()>=0. and P()>=0.) {
    fMass = std::sqrt(std::pow(E(), 2)-std::pow(P(), 2));
    return true;
  }
  return false;
}

void
Particle::SetMother(Particle* part_)
{
  fMothers.insert(part_->id);
  fIsPrimary = false;
  
  Debug(Form("Particle %2d (pdgId=%4d) is the new mother of %2d (pdgId=%4d)",
                  part_->id+1, part_->pdgId, id+1, pdgId));
  
  part_->AddDaughter(this);
};

bool
Particle::AddDaughter(Particle* part_)
{
  std::pair<ParticlesIds::iterator,bool> ret;
  ret = fDaughters.insert(part_->id);

  if (Logger::GetInstance()->Level>=Logger::Debug) {
    std::ostringstream os;
    ParticlesIds::iterator it;
    for (it=fDaughters.begin(); it!=fDaughters.end(); it++) os << Form("\n\t * %d", *it);
    Debug(Form("Particle %2d (pdgId=%4d) has now %2d daughter(s):"
                    "%s", role, pdgId, NumDaughters(), os.str().c_str()));
  }
  
  if (ret.second) {
    Debug(Form("Particle %2d (pdgId=%4d) is a new daughter of %2d (pdgId=%4d)",
                    part_->role, part_->pdgId, role, pdgId));
    
    if (!part_->Primary() && part_->GetMothersIds().size()<1) {
      part_->SetMother(this);
    }
  }

  return ret.second;
}

std::vector<int>
Particle::GetDaughters()
{
  std::vector<int> out;
  ParticlesIds::iterator it;
  
  if (fDaughters.empty()) return out;
  
  out.reserve(fDaughters.size());
  
  DebugInsideLoop(Form("Reserved %d slot(s) for the daughter particle(s)", fDaughters.size()));
  
  for (it=fDaughters.begin(); it!=fDaughters.end(); it++) {
    if (*it==-1) continue;
    out.push_back(*it);
  }
  std::sort(out.begin(), out.end());
  
  DebugInsideLoop(Form("Returning a vector containing %d particle(s)", out.size()))
  
  return out;
}
 
void
Particle::Dump()
{  
  if (!Valid()) throw Exception(__PRETTY_FUNCTION__, Form("Particle with role \"%d\" is invalid", role), Fatal);
  
  std::vector<int> daugh;
  std::ostringstream osm, osd;
  if (!Primary()) {
    ParticlesIds::iterator m;
    osm << "\n\t\tMothers = ";
    for (m=fMothers.begin(); m!=fMothers.end(); m++) osm << (*m) << " ";
  }
  daugh = GetDaughters();
  for (unsigned int i=0; i<NumDaughters(); i++) osd << "\n\t\t* Id = " << daugh[i];
  Info(Form("Id:\t%4d\t"
            "role:\t%4d\n\t"
            "Status:\t%4d\t"
            "PDG Id:\t%4d\n\t"
            "(E,P) = (%4.2f, %4.2f, %4.2f, %4.2f) GeV\t"
            "(|P| = p = %4.2f GeV)\n\t"
            " Pt = %4.2f GeV\teta = %4.3f\n\t"
            " M = %4.2f GeV\n\t"
            "Primary? %d%s\n\t"
            "Daughters (%d)",
            id, role, status, pdgId,
            E(), Px(), Py(), Pz(), P(), Pt(), Eta(), M(),
            Primary(), osm.str().c_str(), NumDaughters(), osd.str().c_str()));
}

//double*
void
Particle::LorentzBoost(double m_, double p_[4])
{
  double pf4, fn;

  if (p_[3]!=m_) {
    pf4 = 0.;
    for (int i=0; i<4; i++) {
      pf4 += fP4[i]*p_[i];
    }
    pf4 /= m_;
    fn = (pf4+E())/(p_[3]+m_);
    /*for (int i=0; i<3; i++) {
      __tmp3[i] = fP4[i] + fn*p_[i];
    }*/
    for (int i=0; i<3; i++) {
      fP4[i] += fn*p_[i];
    }
  }
  //return __tmp3;
}

double*
Particle::LorentzBoost(double p_[3])
{
  double p2, gamma, bp, gamma2;

  p2 = std::pow(p_[0], 2)+std::pow(p_[1], 2)+std::pow(p_[2], 2);
  gamma = 1./std::sqrt(1.-p2);
  bp = 0.;
  for (int i=0; i<3; i++) bp+= p_[i]*fP4[i];
  //bp = p_[0]*fP4[0]+p_[1]*fP4[1]+p_[2]*fP4[2];

  if (p2>0.) gamma2 = (gamma-1.)/p2;
  else gamma2 = 0.;

  for (int i=0; i<3; i++) {
    __tmp3[i] = fP4[i] + gamma2*bp*p_[i]+gamma*p_[i]*E();
  }
  //E(gamma*E()+bp);
  return __tmp3;
}

void
Particle::RotateThetaPhi(double theta_, double phi_)
{
  double rotmtx[3][3], mom[3]; //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
  rotmtx[0][0] = -sin(phi_); rotmtx[0][1] = -cos(theta_)*cos(phi_); rotmtx[0][2] =  sin(theta_)*cos(phi_);
  rotmtx[1][0] =  cos(phi_); rotmtx[1][1] = -cos(theta_)*sin(phi_); rotmtx[1][2] =  sin(theta_)*sin(phi_);
  rotmtx[2][0] =  0.;        rotmtx[2][1] =  sin(theta_);           rotmtx[2][2] =  cos(theta_);

  for (int i=0; i<3; i++) {
    mom[i] = 0.;
    for (int j=0; j<3; j++) {
      mom[i] += rotmtx[i][j]*fP4[j];
    }
  }

  std::copy(mom, mom+3, fP4);
  //fP4[0] *= sin(theta_)*cos(phi_);
  //fP4[1] *= sin(theta_)*sin(phi_);
  //fP4[2] *= cos(theta_);
}

double
Particle::GetMassFromPDGId(Particle::ParticleCode pdgId_)
{
  switch (abs(pdgId_)) {
  case dQuark:      return 0.33;           // mass from PYTHIA6.4
  case uQuark:      return 0.33;           // mass from PYTHIA6.4
  case Electron:    return 0.510998928e-3;
  case Muon:        return 0.1056583715;
  case Tau:         return 1.77682;
  case Gluon:       return 0.;
  case Photon:      return 0.;
  case PiPlus:      return 0.13957018;
  case PiZero:      return 0.1349766;
  case JPsi:        return 20.;            // J/psi //FIXME FIXME FIXME
  case ud0Diquark:  return 0.57933;
  case ud1Diquark:  return 0.77133;
  case uu1Diquark:  return 0.77133;
  case Proton:      return 0.938272046;
  case Neutron:     return 0.939565346;
  default:          return -1.;
  }
}

double
Particle::GetWidthFromPDGId(Particle::ParticleCode pdgId_)
{
  switch (abs(pdgId_)) {
  case JPsi:  return 5.; //FIXME
  default:     return -1.;
  }
}

std::ostream&
operator<<(std::ostream& os, const Particle::ParticleCode& pc)
{
  switch (pc) {
  case Particle::dQuark:     os << "d quark"; break;
  case Particle::uQuark:     os << "u quark"; break;
  case Particle::Electron:   os << "Electron"; break;
  case Particle::Muon:       os << "Muon"; break;
  case Particle::Tau:        os << "Tau"; break;
  case Particle::Gluon:      os << "Gluon"; break;
  case Particle::Photon:     os << "Photon"; break;
  case Particle::PiPlus:     os << "Pi+"; break;
  case Particle::PiZero:     os << "Pi0"; break;
  case Particle::Rho770_0:   os << "Rho(770)0"; break;
  case Particle::Omega782:   os << "Omega(782)"; break;
  case Particle::JPsi:       os << "J/Psi"; break;
  case Particle::Phi1680:    os << "Phi(1680)"; break;
  case Particle::Upsilon1S:  os << "Upsilon(1S)"; break;
  case Particle::Upsilon2S:  os << "Upsilon(2S)"; break;
  case Particle::Upsilon3S:  os << "Upsilon(3S)"; break;
  case Particle::ud0Diquark: os << "(ud)0 di-quark"; break;
  case Particle::ud1Diquark: os << "(ud)1 di-quark"; break;
  case Particle::uu1Diquark: os << "(uu)1 di-quark"; break;
  case Particle::Proton:     os << "Proton"; break;
  case Particle::Neutron:    os << "Neutron"; break;
  case Particle::Pomeron:    os << "Pomeron"; break;
  case Particle::Reggeon:    os << "Reggeon"; break;
  }
  return os;
}
