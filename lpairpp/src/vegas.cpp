#include "vegas.h"

//#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])

Vegas::Vegas(const int dim_, double f_(double*,size_t,void*), Parameters* inParam_) :
  _ndim(dim_), _ndo(50),
  _nTreatCalls(0), _mbin(3),
  _ffmax(0.), _correc(0.), _corre2(0.), _fmax2(0.), _fmdiff(0.), _fmold(0.),
  _j(0),
  _force_correction(false)
{
  /* x content :
      0 = t1 mapping
      1 = t2 mapping
      2 = s2 mapping
      3 = yy4 definition
      4 = w4 mapping
      5 = xx6 definition
      6 = phicm6 definition
    ( 7 = xq, wx mappings   ) <- single- and double-dissociative only
  */

  _xl = new double[dim_];
  _xu = new double[dim_];
  for (int i=0; i<MAX_ND; i++) {
    _xi[i] = new double[dim_];
    _d[i] = new double[dim_];
    _di[i] = new double[dim_];
  }
  // ...
  _n = new int[dim_];
  _nm = new int[7000];
  _fmax = new double[7000];
  // ...
  
  for (int i=0; i<dim_; i++) {
    _xl[i] = 0.;
    _xu[i] = 1.;
  }
  _nIter = inParam_->itvg;
  _ncalls = inParam_->ncvg;

#ifdef DEBUG
  std::cout << "[Vegas::Vegas] [DEBUG]"
            << "\n  Number of integration dimensions : " << dim_
            << "\n  Number of iterations : " << _nIter
            << "\n  Number of function calls : " << _ncalls
            << std::endl;
#endif

  _ip = inParam_;
  _f = f_;
}

Vegas::~Vegas()
{
#ifdef DEBUG
  std::cout << "[Vegas::~Vegas] [DEBUG] Destructor called" << std::endl;
#endif
  delete[] _xl;
  delete[] _xu;
  delete[] _n;
  delete[] _nm;
  delete[] _fmax;
  for (int i=0; i<MAX_ND; i++) {
    delete[] _xi[i];
    delete[] _d[i];
    delete[] _di[i];
  }
}

int Vegas::Integrate(double *result_, double *abserr_)
{
  unsigned int ndm, nd, npg;
  int ng, now;
  unsigned int i, j;
  int k;
  const double one = 1.;
  int it, iaj, iaj1;
  double avgi, sd, tsi, ti;
  double chi2a, rel, ti2, si, si2, swgt, schi, scalls;
  double calls, dxg, dv2g, xnd, xo, xn, xjac, rc, dr;
  double f, f2, wgt;
  double fb, f2b;
  double x[_ndim], dx[_ndim], kg[_ndim], dt[_ndim], qran[_ndim], xin[MAX_ND], r[MAX_ND];
  int ia[_ndim];

  for (j=0; j<MAX_ND; j++) {
    xin[j] = r[j] = 0.;
    for (i=0; i<_ndim; i++) {
      _d[j][i] = _di[j][i] = _xi[j][i] = 0.;
    }
  }
  for (i=0; i<_ndim; i++) {
    x[i] = dx[i] = kg[i] = dt[i] = qran[i] = 0.;
  }

  // Parameters
  int mds = 1;
  double acc = 1.e-4;
  double alph = 1.5;

  if (_ip->itvg<0) {
    std::cerr << "[Vegas::MyIntegrate] [ERROR] Vegas called with a negative number of maximum iterations. No execution." << std::endl;
    return -1;
  }
  _ndo = 1;
  for (j=0; j<_ndim; j++) {
    _xi[0][j] = one;
  }

  // Entry VEGAS1
  it = 0;
  si = si2 = swgt = schi = scalls = 0.;

  // Entry VEGAS2
  nd = MAX_ND;
  ng = 1;
  if (mds!=0) {
    ng = std::pow(_ip->ncvg/2., 1./_ndim);
    mds = 1;
    if (2*ng-MAX_ND>=0) {
      mds = -1;
      npg = ng/MAX_ND+1;
      nd = ng/npg;
      ng = npg*nd;
    }
  }

  k = std::pow(ng, _ndim);
  npg = _ip->ncvg/k;
  if (npg<2) npg = 2;
  calls = npg*k;
  dxg = one/ng;
  dv2g = std::pow(dxg, 2*_ndim)/npg/npg/(npg-one);
  xnd = nd;
  ndm = nd-1;
  dxg *= xnd;
  xjac = one;
  //std::cout << k << "\t" << ng << "\t" << calls << "\t" << dxg << "\t" << dv2g << "\t" << xnd << "\t" << ndm << std::endl;
  //exit(0);
  for (i=0; i<_ndim; i++) {
    dx[i] = _xu[i]-_xl[i];
    xjac *= dx[i];
    //std::cout << i << "\t" << xjac << "\t" << dx[i] << std::endl;
  }

  // Rebin preserving bin density
  if (nd!=_ndo) {
    rc = _ndo/xnd;
    for (j=0; j<_ndim; j++) {
      k = 0;
      xo = dr = xn = 0.;
      i = -1;
      do {
	dr += one;
	xn = _xi[k][j];
	k++;
      line5:
	now += 1;
      } while (rc>dr);
      i++;
      dr -= rc;
      xin[i] = xn-(xn-xo)*dr;
      //std::cout << i << "\t" << j << "\t" << k << "\t" << dr << "\t" << xin[i] << std::endl;
      if (i<ndm-1) goto line5; //FIXME need to remove these gotos
      for (i=0; i<ndm; i++) {
	_xi[i][j] = xin[i];
	//std::cout << i << "\t" << j << "\t" << _xi[i][j] << std::endl;
      }
      _xi[nd-1][j] = one;
    }
    _ndo = nd;
  }

  // Entry VEGAS3
  // Main integration loop
  do {
    it++;
    tsi = ti = 0.;
    for (j=0; j<_ndim; j++) {
      kg[j] = 1;
      for (i=1; i<nd; i++) {
	_di[i][j] = _d[i][j] = ti;
      }
    }
  line11:
    f2b = fb = 0.;
    for (k=0; k<(int)npg; k++) {
      for (j=0; j<_ndim; j++) {
	qran[j] = (double)rand()/RAND_MAX;
	//if (j%2==0) qran[j] = 0.7;
	//else qran[j] = 0.78;
	//qran[j] = 0.7;
      }
      wgt = xjac;
      for (j=0; j<_ndim; j++) {
	//iaj = ia[j] = xn = (kg[j]-qran[j])*dxg+one;
	xn = (kg[j]-qran[j])*dxg;
	iaj = ia[j] = int(xn);
	iaj1 = iaj-1;
	if (iaj<=1) {
	  xo = _xi[iaj][j];
	  rc = (xn-iaj)*xo;
	}
	else {
	  xo = _xi[iaj][j]-_xi[iaj1][j];
	  rc = _xi[iaj1][j]+(xn-iaj)*xo;
	}
	x[j] = _xl[j]+rc*dx[j];
	if (x[j]>1. or x[j]<0.)
	  std::cout << "-------> j=" << j 
		    << "\tx[j]=" << x[j] 
		    << "\txo=" << xo 
		    << "\trc=" << rc 
		    << "\tiaj=" << iaj 
		    << "\t(xn-iaj)=" << (xn-iaj) 
		    << "\txi[iaj1][j]=" << _xi[iaj1][j]
		    << "\txi[iaj][j]=" << _xi[iaj][j]
		    << std::endl;
	wgt *= (xo*xnd);
      }
      f = this->F(x)*wgt;
      
      f2 = std::pow(f, 2);
      fb += f;
      f2b += f2;
      
      for (j=0; j<_ndim; j++) {
	iaj = ia[j];
	_di[iaj][j] += f/calls;
	if (mds>=0) _d[iaj][j] += f2;
      }
    }

    f2b *= npg;
    f2b = sqrt(f2b);
    f2b = fabs((f2b-fb)*(f2b+fb));
    ti += fb;
    tsi += f2b;
    if (mds<0) {
      for (j=0; j<_ndim; j++) {
	iaj = ia[j];
	_d[iaj][j] += f2b;
      }
    }
    for (k=_ndim-1; k>=0; k--) {
      kg[k] = fmod(kg[k], ng)+1;
      if (kg[k]!=1) goto line11; //FIXME need to remove these goto
    }
    /*for (j=0; j<_ndim; j++) {
      for (i=0; i<ndm; i++) {
	std::cout << "->" << i << "\t" << j << "\t" << _d[i][j] << std::endl;
      }
    }
    exit(0);*/
    
    // Final results for this iteration
    ti /= calls;
    tsi *= dv2g;
    ti2 = std::pow(ti, 2);
    
    if (tsi==0) wgt = 0.;
    else wgt = ti2/tsi;
    
    si += ti*wgt;
    si2 += ti2;
    swgt += wgt;
    schi += ti2*wgt;
    
    if (swgt==0) avgi = ti;
    else avgi = si/swgt;
    
    if (si2==0) sd = tsi;
    else sd = swgt*it/si2;
    
    scalls += calls;
    chi2a = 0.;
    
    if (it>1) chi2a = sd*(schi/swgt-pow(avgi, 2))/(it-1);
    
    if (sd!=0.) sd = sqrt(one/sd);
    else sd = tsi;
    
    std::cout << "--> iteration " 
	      << std::setfill(' ') << std::setw(3) << it << " : "
	      << "average = " << std::setprecision(5) << std::setw(14) << avgi 
	      << "sigma = " << std::setprecision(5) << std::setw(14) << sd 
	      << "chi2 = " << chi2a << std::endl;

    // Refine grid
    if (sd!=0.) rel = fabs(sd/avgi);
    else rel = 0.;
    
    if (rel<=fabs(acc) or it>=_ip->itvg) now = 2;
    
    for (j=0; j<_ndim; j++) {
      xo = _d[0][j];
      xn = _d[1][j];
      _d[0][j] = (xo+xn)/2.;
      dt[j] = _d[0][j];
      for (i=1; i<ndm; i++) {
	_d[i][j] = xo+xn;
	xo = xn;
	xn = _d[i+1][j];
	_d[i][j] = (_d[i][j]+xn)/3.;
	dt[j] += _d[i][j];
      }
      _d[nd][j] = (xn+xo)/2.;
      dt[j] += _d[nd][j];
    }
    
    for (j=0; j<_ndim; j++) {
      rc = 0.;
      for (i=0; i<nd; i++) {
	r[i] = 0.;
	if (_d[i][j]>0.) {
	  xo = dt[j]/_d[i][j];
	  r[i] = std::pow((xo-one)/xo/log(xo), alph);
	}
	rc += r[i];
      }
      rc /= xnd;
      dr = xn = 0.;
      k = 0;
      i = 0;
      do {
	dr += r[k];
	xo = xn;
	xn = _xi[k][j];
	k++;
      line26:
	now += 1;
      } while (rc>dr);
      //for (i=0; i<ndm; i++) {
      dr -= rc;
      if (dr==0.) xin[i] = xn;
      else xin[i] = xn-(xn-xo)*dr/r[k-1];
      i++;
      if (i<ndm) goto line26; //FIXME need to remove these goto

      for (i=0; i<ndm; i++) {
	_xi[i][j] = xin[i];
      }
      _xi[nd-1][j] = one;
    }
  } while (it<_ip->itvg and fabs(acc)<rel);

  *result_ = avgi;
  *abserr_ = sd;
  return 0;
}

void
Vegas::Generate()
{
  //count_ = 1;
  std::ofstream of;
  std::string fn;
  int i;
  
  //fn = "test";
  //of.open(fn.c_str());
  //this->SetGen(_ip->file);
  this->SetGen();
  std::cout << "[Vegas::Generate] [DEBUG] " << _ip->maxgen << " events will be generated" << std::endl;
  i = 0;
  while (i<_ip->maxgen) {
    if (this->GenerateOneEvent()) i++;
  }
  std::cout << "[Vegas::Generate] [DEBUG] " << i << " events generated" << std::endl;
  //of.close();
}

bool
Vegas::GenerateOneEvent()
{
  // Inherited from GMUGNA
  double ami, max;
  double y;
  int jj, jjj;
  double x[_ndim];
  
  y = -1.;
  ami = 1./_mbin;
  max = pow(_mbin, _ndim);

  // Correction cycles are started
  if (_j!=0) {
  line4:
#ifdef DEBUG
    std::cout << "[Vegas::GenerateOneEvent] [DEBUG] Correction cycles are started."
	      << "\n\tj = " << _j
	      << "\n\tcorrec = " << _correc
	      << "\n\tcorre2 = " << _corre2
	      << std::endl;
#endif
    if (_correc<1.) {
      if ((double)rand()/RAND_MAX>=_correc) {
	goto line7; //FIXME need to remove these goto
      }
      _correc = -1.;
    }
    else {
      _correc -= 1.;
    }
    // Select x values in Vegas bin
    for (unsigned int k=0; k<_ndim; k++) {
      x[k] = ((double)rand()/RAND_MAX+_n[k])*ami;
    }
    // Compute weight for x value
    if (_ip->ntreat>0) _weight = Treat(x);
    else _weight = this->F(x);
    // Parameter for correction of correction
    if (_weight>_fmax[_j]) {
      if (_weight>_fmax2) _fmax2 = _weight;
      _corre2 -= 1.;
      _correc += 1.;
    }
    // Accept event
    if (_weight>=_fmdiff*(double)rand()/RAND_MAX+_fmold) { // FIXME!!!!
      return this->StoreEvent(x);
    }
    goto line4;
  line7:
    // Correction if too big weight is found while correction
    // (All your bases are belong to us...)
    if (_fmax2>_fmax[_j]) {
      _fmold = _fmax[_j];
      _fmax[_j] = _fmax2;
      _fmdiff = _fmax2-_fmold;
      if (_fmax2<_ffmax) {
	_correc = (_nm[_j]-1.)*_fmdiff/_ffmax-_corre2;
      }
      else {
	_ffmax = _fmax2;
	_correc = (_nm[_j]-1.)*_fmdiff/_ffmax*_fmax2/_ffmax-_corre2;
      }
      _corre2 = 0.;
      _fmax2 = 0.;
      goto line4;
      //return this->GenerateOneEvent(); //GOTO 4
    }
  }

  // Normal generation cycle
  // Select a Vegas bin and reject if fmax is too little
 line1:
  //double* Vegas::SelectBin() //FIXME need to implement it this way instead of these bloody goto...!
  do {
    // ...
    _j = (double)rand()/RAND_MAX*max;
    y = (double)rand()/RAND_MAX*_ffmax;
    _nm[_j] += 1;
  } while (y>_fmax[_j]);
  // Select x values in this Vegas bin
  jj = _j;
  for (unsigned int i=0; i<_ndim; i++) {
    jjj = jj/_mbin;
    _n[i] = jj-jjj*_mbin;
    x[i] = ((double)rand()/RAND_MAX+_n[i])*ami;
    jj = jjj;
  }

  // Get weight for selected x value
  if (_ip->ntreat>0) _weight = this->Treat(x);
  else _weight = this->F(x);

  // Eject if weight is too low
  if (y>_weight) {
    //std::cout << "ERROR : y>weight => " << y << ">" << _weight << ", " << _j << std::endl;
    //_force_correction = false;
    //_force_returnto1 = true;
    //return this->GenerateOneEvent();
    goto line1;
  }

  if (_weight<=_fmax[_j]) _j = 0;
  // Init correction cycle if weight is higher than fmax or ffmax
  else if (_weight<=_ffmax) {
    _fmold = _fmax[_j];
    _fmax[_j] = _weight;
    _fmdiff = _weight-_fmold;
    _correc = (_nm[_j]-1.)*_fmdiff/_ffmax-1.;
  }
  else {
    _fmold = _fmax[_j];
    _fmax[_j] = _weight;
    _fmdiff = _weight-_fmold;
    _ffmax = _weight;
    _correc = (_nm[_j]-1.)*_fmdiff/_ffmax*_weight/_ffmax-1.;
  }
#ifdef DEBUG
  std::cout << "[Vegas::GenerateOneEvent] [DEBUG] correc = " << _correc << ", j = " << _j << std::endl;
#endif
  // Return with an accepted event
  return this->StoreEvent(x);
}

bool
Vegas::StoreEvent(double *x_)
{
  if (_weight<=0.) {
#ifdef DEBUG
    std::cout << "[Vegas::StoreEvent] [DEBUG] Tried to store event while the weight is <= 0 : " << _weight << std::endl;
#endif
    return false;
  }
  _ip->store = true;
  //if (_ip->ntreat>0) _weight = Treat(x_,true);
  if (_ip->ntreat>0) _weight = Treat(x_);
  else _weight = this->F(x_);
  _ip->ngen += 1;
  _ip->store = false;
  if (_ip->ngen%1000==0) {
    std::cout << "[Vegas::StoreEvent] Generated events : " << _ip->ngen << std::endl;
  }
  return true;
}

void
Vegas::SetGen()
{
  int max;
  int jj, jjj;
  double sum, sum2, sum2p;
  int n[10];
  int npoin = _ip->npoints;
  double fsum, fsum2;
  double z;
  double x[_ndim];
  double sig2;
  double av, av2;
#ifdef DEBUG
  double eff, eff1, eff2;
  double sig, sigp;
#endif

  _ip->ngen = 0;
#ifdef DEBUG
  std::cout << "[Vegas::SetGen] [DEBUG] maxgen = " << _ip->maxgen << std::endl;
  _ip->Dump();
#endif
  // ...
  sum = 0.;
  sum2 = 0.;
  sum2p = 0.;
  max = pow(_mbin, _ndim);

  for (int i=0; i<max; i++) {
    _nm[i] = 0;
    _fmax[i] = 0.;
  }

  for (int i=1; i<=max; i++) {
    jj = i-1;
    for (unsigned int j=1; j<=_ndim; j++) {
      jjj = jj/_mbin;
      n[j-1] = jj-jjj*_mbin;
      jj = jjj;
    }
    fsum = fsum2 = 0.;
    for (int j=1; j<=npoin; j++) {
      for (unsigned int k=1; k<=_ndim; k++) {
        x[k-1] = ((double)rand()/RAND_MAX+n[k-1])/_mbin;
      }
      if (_ip->ntreat>0) z = this->Treat(x);
      else z = this->F(x);
      if (z>_fmax[i-1]) _fmax[i-1] = z;
      fsum += z;
      fsum2 += std::pow(z, 2);
    }
    av = fsum/npoin;
    av2 = fsum2/npoin;
    sig2 = av2-pow(av, 2);
#ifdef DEBUG
    sig = sqrt(sig2);
#endif
    sum += av;
    sum2 += av2;
    sum2p += sig2;
    if (_fmax[i-1]>_ffmax) _ffmax = _fmax[i-1];
#ifdef DEBUG
    eff = 1.e4;
    if (_fmax[i-1]!=0.) eff = _fmax[i-1]/av;
    std::cout << "[Vegas::SetGen] [DEBUG] in iteration #" << i << " :"
	      << "\n\tav   = " << av
	      << "\n\tsig  = " << sig
	      << "\n\tfmax = " << _fmax[i]
	      << "\n\teff  = " << eff
	      << "\n\tn = (";
    for (unsigned int j=0; j<_ndim; j++) {
      std::cout << n[j];
      if (j!=_ndim-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
#endif
  }

  sum = sum/max;
  sum2 = sum2/max;
  sum2p = sum2p/max;

#ifdef DEBUG
  sig = sqrt(sum2-pow(sum, 2));
  sigp = sqrt(sum2p);
  eff1 = 0.;
  for (int i=0; i<max; i++) {
    eff1 += _fmax[i];
  }
  eff1 = eff1/(max*sum);
  eff2 = _ffmax/sum;
  std::cout << "[Vegas::SetGen] [DEBUG]"
            << "\n\tAverage function value     =  sum   = " << sum
            << "\n\tAverage function value**2  =  sum2  = " << sum2
            << "\n\tOverall standard deviation =  sig   = " << sig
	    << "\n\tAverage standard deviation =  sigp  = " << sigp
	    << "\n\tMaximum function value     = ffmax  = " << _ffmax
	    << "\n\tAverage inefficiency       =  eff1  = " << eff1 
	    << "\n\tOverall inefficiency       =  eff2  = " << eff2 
            << "\n\teff = " << eff 
	    << std::endl;
#endif
}

void
Vegas::DumpGrid()
{
  unsigned int i,j;
  // DUMP THE GRID
  for (i=0; i<_ndim; i++) {
    for (j=0; j<MAX_ND; j++) {
      std::cout << i << "\t" << j << "\t" << _xi[j][i] << std::endl;
    }
  }
}

double
Vegas::Treat(double *x_, Parameters* ip_, bool storedbg_)
{
  double w, xx, y, dd, f;
  int j;
  double z[_ndim];

  if (_nTreatCalls==0) {
    _nTreatCalls = 1;
    _rTreat = std::pow(_ndo, _ndim);
    if (storedbg_ && remove("test_vegas")!=0) {
      std::cerr << "Error while trying to delete test_vegas" << std::endl;
    }
    //this->DumpGrid();
  }
  w = _rTreat;
  //std::cout << "rtreat = " << _rTreat << std::endl;
  unsigned int i;
  for (i=0; i<_ndim; i++) {
    xx = x_[i]*_ndo-1;
    j = xx;
    y = xx-j;
    if (j<=0) {
      dd = _xi[0][i];
    }
    else {
      dd = _xi[j+1][i]-_xi[j][i];
    }
    z[i] = _xi[j+1][i]-dd*(1.-y);
    w = w*dd;
  }

  f = this->F(z, ip_);

  if (storedbg_) {
    std::ofstream df;
    df.open("test_vegas", std::ios::app);
    df << w 
       << "\t" << w*f;
    for (unsigned int i=0; i<_ndim; i++) {
      df << "\t" << z[i];
    }
    for (unsigned int i=0; i<_ndim; i++) {
      df << "\t" << x_[i];
    }
    df << std::endl;
    df.close();
  }
#ifdef DEBUG
  std::cout << "[Vegas::Treat] [DEBUG] w = " << w << ", dd = " << dd << ", ndo = " << _ndo << ", r = " << _rTreat << std::endl;
#endif
  return w*f;
}
