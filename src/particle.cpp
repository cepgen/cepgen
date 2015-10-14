#include "particle.h"

Particle::Particle() :
  id(-1), pdgId((ParticleCode)0), charge(999.), name(""), role(-1),
  helicity(0.),
  status(0), _m(-1.),
  _isPrimary(true)
{
  for (int i=0; i<4; i++) _p4[i] = 0.;
}

Particle::Particle(int role_, ParticleCode pdgId_) :
  id(-1), pdgId(pdgId_), charge(999.), name(""), role(role_),
  status(0), _m(-1.),
  _isPrimary(true)
{
  for (int i=0; i<4; i++) _p4[i] = 0.;
  if (this->pdgId!=0) {
    this->M(-1.);
  }
}

Particle::~Particle()
{
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Destructor called" << std::endl;
#endif
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
  this->M(part_._m);

  return *this;
}

Particle&
Particle::operator+=(const Particle &part_)
{
  pdgId = (part_.pdgId==pdgId) ? pdgId : (ParticleCode)(-1);
  role = (part_.role==role) ? role : -1;
  for (int i=0; i<4; i++) _p4[i]+= part_._p4[i];
  M(-1.);

  return *this;
}

Particle&
Particle::operator-(const Particle &part_)
{
  pdgId = (part_.pdgId==pdgId) ? pdgId : (ParticleCode)(-1);
  role = (part_.role==role) ? role : -1;
  for (int i=0; i<4; i++) _p4[i]-= part_._p4[i];
  M(-1.);

  return *this;
}

bool
Particle::Valid()
{
  if (this->pdgId==0) return false;
  if (this->P()==0. and this->M()==0.) return false;
  return true;
}

std::string
Particle::GetLHEline(bool revert_)
{
  std::stringstream line;

  if (revert_) {
    _p4[2] = -_p4[2];
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
  if (m_>=0.) this->_m = m_;
  else if (this->pdgId!=0) {
    mass = GetMassFromPDGId(this->pdgId);
    if (mass<0.) return false;
    if (this->_p4[3]<0.) {
      this->_m = mass;
      this->_p4[3] = std::pow(this->P(), 2)+this->M2();
      return true;
    }
    if (std::pow(this->E(), 2)-std::pow(this->P(), 2)!=std::pow(mass, 2)) mass = std::sqrt(std::pow(this->E(), 2)-std::pow(this->P(), 2));
    this->_m = mass;
    return true;
  }
  else if (this->E()>=0. and this->P()>=0.) {
    this->_m = std::sqrt(std::pow(this->E(), 2)-std::pow(this->P(), 2));
    return true;
  }
  return false;
}

void
Particle::SetMother(Particle* part_)
{
  this->_moth.insert(part_->id);
  this->_isPrimary = false;
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Particle "
	    << part_->id+1 << " (pdgId=" << part_->pdgId << ") is the new mother of "
	    << this->id+1 << " (pdgId=" << this->pdgId << ")" << std::endl;
#endif
  part_->AddDaughter(this);
};

bool
Particle::AddDaughter(Particle* part_)
{
  std::pair<ParticlesIds::iterator,bool> ret;
  ret = this->_daugh.insert(part_->id);
#ifdef DEBUG
  ParticlesIds::iterator it;
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Particle "
	    << this->role << " (pdgId=" << this->pdgId << ") has now "
	    << this->NumDaughters() << " daughter(s) : " << std::endl;
  for (it=this->_daugh.begin(); it!=this->_daugh.end(); it++) {
    std::cout << " * " << *it << std::endl;
  }
#endif

  if (ret.second) {
#ifdef DEBUG
    std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Particle " 
	      << part_->role << " (pdgId=" << part_->pdgId << ") is a new daughter of "
	      << this->role << " (pdgId=" << this->pdgId << ")" << std::endl;
#endif
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
  
  if (this->_daugh.empty()) return out;
  
  out.reserve(this->_daugh.size());
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Reserved " << this->_daugh.size() << " slot(s) for the daughter particle(s)" << std::endl;
#endif
  
  for (it=this->_daugh.begin(); it!=this->_daugh.end(); it++) {
    if (*it==-1) continue;
#ifdef DEBUG
    std::cout << " * " << *it << std::endl;
#endif
    out.push_back(*it);
  }
  std::sort(out.begin(), out.end());
#ifdef DEBUG
  std::cout << __PRETTY_FUNCTION__ << " [DEBUG] Returning a vector containing " << out.size() << " particle(s)" << std::endl;
#endif
  return out;
}
 
void
Particle::Dump()
{
  std::vector<int> daugh;

  if (this->Valid()) {
    std::cout << __PRETTY_FUNCTION__ << std::endl
	            << "  Id = " << this->id << std::endl
	            << "  Role = " << this->role << std::endl
	            << "  Status = " << this->status << std::endl
	            << "  PDG id = " << this->pdgId << std::endl
	            << "  P = (" << this->Px() << ", "
                           << this->Py() << ", "
                           << this->Pz() << ") GeV" << std::endl
	            << "  |P| = " << this->P() << " GeV" << std::endl
	            << "  Pt = " << this->Pt() << " GeV" << std::endl
	            << "  E = " << this->E() << " GeV" << std::endl
	            << "  M = " << this->M() << " GeV" << std::endl
	            << "  eta = " << this->Eta() << std::endl
	            << "  Is valid ? " << this->Valid() << std::endl
	            << "  Is primary ? " << this->_isPrimary << std::endl;
    if (!this->Primary()) {
      ParticlesIds::iterator m;
      std::cout << "  Mothers = ";
      for (m=_moth.begin(); m!=_moth.end(); m++) std::cout << (*m) << " ";
      std::cout << std::endl;
    }

    daugh = this->GetDaughters();

    std::cout << "  Daughters (" << this->NumDaughters() << ")" << std::endl;
    for (unsigned int i=0; i<this->NumDaughters(); i++) {
      std::cout << "   * Id = " << daugh[i] << std::endl; 
    }
  }
  else {
    std::cout << __PRETTY_FUNCTION__ << " ERROR: Particle with role \"" << this->role << "\" is invalid" << std::endl;
  }
}

//double*
void
Particle::LorentzBoost(double m_, double p_[4])
{
  double pf4, fn;

  if (p_[3]!=m_) {
    pf4 = 0.;
    for (int i=0; i<4; i++) {
      pf4 += this->_p4[i]*p_[i];
    }
    pf4 /= m_;
    fn = (pf4+this->E())/(p_[3]+m_);
    /*for (int i=0; i<3; i++) {
      __tmp3[i] = this->_p4[i] + fn*p_[i];
    }*/
    for (int i=0; i<3; i++) {
      this->_p4[i] += fn*p_[i];
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
  for (int i=0; i<3; i++) bp+= p_[i]*_p4[i];
  //bp = p_[0]*_p4[0]+p_[1]*_p4[1]+p_[2]*_p4[2];

  if (p2>0.) gamma2 = (gamma-1.)/p2;
  else gamma2 = 0.;

  for (int i=0; i<3; i++) {
    __tmp3[i] = this->_p4[i] + gamma2*bp*p_[i]+gamma*p_[i]*E();
  }
  //this->E(gamma*E()+bp);
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
      mom[i] += rotmtx[i][j]*this->_p4[j];
    }
  }

  std::copy(mom, mom+3, this->_p4);
  //this->_p4[0] *= sin(theta_)*cos(phi_);
  //this->_p4[1] *= sin(theta_)*sin(phi_);
  //this->_p4[2] *= cos(theta_);
}

double
Particle::GetMassFromPDGId(Particle::ParticleCode pdgId_)
{
  switch (abs(pdgId_)) {
  case QUARK_D:     return 0.33;           // mass from PYTHIA6.4
  case QUARK_U:     return 0.33;           // mass from PYTHIA6.4
  case ELECTRON:    return 0.510998928e-3;
  case MUON:        return 0.1056583715;
  case TAU:         return 1.77682;
  case GLUON:       return 0.;
  case PHOTON:      return 0.;
  case PI_PLUS:     return 0.13957018;
  case PI_0:        return 0.1349766;
  case J_PSI:       return 20.;            // J/psi //FIXME FIXME FIXME
  case DIQUARK_UD0: return 0.57933;
  case DIQUARK_UD1: return 0.77133;
  case DIQUARK_UU1: return 0.77133;
  case PROTON:      return 0.938272046;
  case NEUTRON:     return 0.939565346;
  default:          return -1.;
  }
}

double
Particle::GetWidthFromPDGId(Particle::ParticleCode pdgId_)
{
  switch (abs(pdgId_)) {
  case J_PSI:  return 5.; //FIXME
  default:     return -1.;
  }
}
