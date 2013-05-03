#include "../include/gnuplot.h"

Gnuplot::Gnuplot(std::string outFile_)
{
  // with -persist option you will see the windows as your program ends
  _pipe = popen("gnuplot -persist","w");
  //without that option you will not see the window
  // because I choose the terminal to output files so I don't want to see the window
  //_pipe=popen("gnuplot","w");
  if (!_pipe) {
    std::cerr << "Gnuplot not found !" << std::endl;
  }
  if (outFile_!="") {
    this->SetOutputFile(outFile_);
  }
}

Gnuplot::~Gnuplot()
{
  fprintf(_pipe,"exit\n");
  pclose(_pipe);
}

void Gnuplot::SetOutputFile(std::string outFile_)
{
  //fprintf(_pipe,"set term epslatex\n");
  fprintf(_pipe,"set term pngcairo\n");
  fprintf(_pipe,"set key left bottom\n");
  fprintf(_pipe,"set output '%s'\n",outFile_.c_str());
  fflush(_pipe);
}

void Gnuplot::SetTitle(std::string title_)
{
  fprintf(_pipe,"set title '%s'\n",title_.c_str());
  fflush(_pipe);
}

void Gnuplot::SetXAxisTitle(std::string title_)
{
  fprintf(_pipe,"set xlabel '%s'\n",title_.c_str());
  fflush(_pipe);
}

void Gnuplot::SetYAxisTitle(std::string title_)
{
  fprintf(_pipe,"set ylabel '%s'\n",title_.c_str());
  fflush(_pipe);
}

void Gnuplot::SetLogy(bool logy_)
{
  if (logy_) {
    fprintf(_pipe,"set logscale y\n");
  }
  else {
    fprintf(_pipe,"unset logscale y\n");  
  }
  fflush(_pipe);
}

void Gnuplot::SetHistogram(int num_, double low_, double high_)
{
  _isHist = true;
  fprintf(_pipe,"set style data histograms\n");
  fprintf(_pipe,"set style histogram gap 0.\n");
  fprintf(_pipe,"set style fill solid 1.0 border -1\n");
  _histBounds = new double[num_+2];
  _histValues = new int[num_];
  _histUnderflow = 0.;
  _histOverflow = 0.;
  _histNum = num_;
  _histLow = low_;
  _histHigh = high_;
  for (int i=0; i<=num_+1; i++) {
    _histBounds[i] = _histLow+i*(_histHigh-_histLow)/_histNum;
  }
}

int Gnuplot::Fill(double value_)
{
  if (value_<_histLow) {
    _histUnderflow += 1;
    return 2;
  }
  if (value_>_histHigh) {
    _histOverflow += 1;
    return 3;
  }
  for (int i=0; i<_histNum; i++) {
    if (value_>=_histBounds[i] && value_<_histBounds[i+1]) {
      _histValues[i] += 1;
      //std::cout << "--> value " << value_ << " in bin " << i << std::endl;
      return 1;
    }
  }
  return 0;
}

int Gnuplot::DrawHistogram()
{
  if (!_isHist) {
    return -1;
  }
  /*for (int i=0; i<_histNum; i++) {
    std::cout << "[" << i << "] = " << _histValues[i] << std::endl;
  }
  std::cout << "underflow : " << _histUnderflow << std::endl;
  std::cout << "overflow : " << _histOverflow << std::endl;*/
  
  std::ofstream tmp;
  tmp.open("/tmp/gp.tmp");
  for (int i=0; i<_histNum; i++) {
    tmp << _histBounds[i] << "\t" << _histValues[i] << std::endl;
  }
  //std::string ef = "everyn(col) = (int(column(col)) %5 == 0) ? stringcolumn(1) : \"\"";
  std::string ef = "everyNth(lab,N) =((int(column(0)) % N == 0) ? stringcolumn(lab) : \"\"); ";
  fprintf(_pipe,"set xtics auto\n");
  fprintf(_pipe,"%s\n",ef.c_str());
  //fprintf(_pipe,"plot '/tmp/gp.tmp' using 2:(everyn(0))\n");
  //fprintf(_pipe,"plot '/tmp/gp.tmp' using 2:xtic(1) t col\n");
  fprintf(_pipe,"plot '/tmp/gp.tmp' using 2:xtic(everyNth(1, 10))\n");
  //fprintf(_pipe,"set xtics -1.,0.5\n");
  //fprintf(_pipe,"replot\n");
  
  fflush(_pipe);
  delete[] _histBounds;
  return 0;
}

void Gnuplot::operator<<(const std::string& command)
{
  fprintf(_pipe,"%s\n",command.c_str());
  fflush(_pipe);
  // flush is necessary, nothing gets plotted else
};

