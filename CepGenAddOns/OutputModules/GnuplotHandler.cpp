#include "CepGen/IO/ExportHandler.h"
#include "CepGen/Utils/String.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>

namespace cepgen
{
  namespace io
  {
    Gnuplot::Gnuplot(std::string outFile_) :
      _type(GP_CLASSIC),
      _isPlottable(false),
      _outputFile("")
    {
      std::stringstream of;

      // with -persist option you will see the windows as your program ends
      _pipe = popen("gnuplot -persist","w");
      //without that option you will not see the window
      // because I choose the terminal to output files so I don't want to see the window
      //_pipe=popen("gnuplot","w");
      if (!_pipe)
        std::cerr << "Gnuplot not found !" << std::endl;
      if (outFile_!="")
        this->SetOutputFile(outFile_);
      of.str(""); of << "/tmp/" << utils::randomString( 5 ) << ".tmp";
      _tmpFile = of.str();
    }

    Gnuplot::~Gnuplot()
    {
      fprintf(_pipe,"exit\n");
      pclose(_pipe);
      remove(_tmpFile.c_str());
    }

    void
    Gnuplot::SetOutputFile(std::string outFile_)
    {
      _outputFile = outFile_;
      //fprintf(_pipe,"set term epslatex\n");
      fprintf(_pipe,"set term pngcairo transparent enhanced font 'arial,10' fontscale 1.0 size 800, 600\n");
      //fprintf(_pipe,"set term pngcairo\n");
      fprintf(_pipe,"set key right top\n");
      //fprintf(_pipe,"set key left bottom\n");
      fprintf(_pipe,"set output '%s'\n",outFile_.c_str());
      fflush(_pipe);
    }

    void Gnuplot::SetTitle(std::string title_)
    {
      _title = title_;
      fprintf(_pipe,"set title '%s'\n",title_.c_str());
      fflush(_pipe);
    }

    void
    Gnuplot::SetName(std::string name_)
    {
      _name = name_;
    }

    void
    Gnuplot::SetXAxisTitle(std::string title_)
    {
      fprintf(_pipe,"set xlabel '%s'\n",title_.c_str());
      fflush(_pipe);
    }

    void
    Gnuplot::SetYAxisTitle(std::string title_)
    {
      fprintf(_pipe,"set ylabel '%s'\n",title_.c_str());
      fflush(_pipe);
    }

    void
    Gnuplot::SetLogy(bool logy_)
    {
      if (logy_)
        fprintf(_pipe,"set logscale y\n");
      else
        fprintf(_pipe,"unset logscale y\n");
      fflush(_pipe);
    }

    void
    Gnuplot::SetGrid(bool grid_)
    {
      if (grid_)
        fprintf(_pipe,"set grid x y mx my\n");
      else
        fprintf(_pipe,"unset grid x y mx my\n");
      fflush(_pipe);
    }

    void
    Gnuplot::SetHistogram(int num_, double low_, double high_, std::string name_)
    {
      _type = GP_HISTOGRAM;
      if (name_!="")
        this->SetName(name_);
      fprintf(_pipe,"set style data histograms\n");
      fprintf(_pipe,"set style histogram gap 0.\n");
      //fprintf(_pipe,"set style fill solid 1.0 border -1\n");
      fprintf(_pipe,"set style fill transparent pattern 2 bo\n");
      _histBounds = new double[num_+2];
      _histValues = new double[num_];
      _histUnderflow = 0.;
      _histOverflow = 0.;
      _histNum = num_;
      _histLow = low_;
      _histHigh = high_;
      for (int i=0; i<=num_+1; i++)
        _histBounds[i] = _histLow+i*(_histHigh-_histLow)/_histNum;
    }

    int
    Gnuplot::Fill(double value_, double weight_)
    {
      _isPlottable = true;
      if (value_<_histLow) {
        _histUnderflow += weight_;
        CG_WARNING("GnuplotHandler") << " value in underflow bin (" << value_ << ").";
        return 2;
      }
      if (value_>_histHigh) {
        _histOverflow += weight_;
        CG_WARNING("GnuplotHandler") << " value in overflow bin (" << value_ << ").";
        return 3;
      }
      for (int i=0; i<_histNum; i++) {
        if (value_>=_histBounds[i] && value_<_histBounds[i+1]) {
          _histValues[i] += weight_;
          CG_DEBUG("GnuplotHandler") << " value in good range (" << value_ << "), bin " << i << ".";
          return 1;
        }
      }
      return 0;
    }

    int
    Gnuplot::DrawHistogram()
    {
      if (_type!=GP_HISTOGRAM || !_isPlottable)
        return -1;

      std::ofstream tmp;
      std::stringstream of;
      std::string ef;

      if (_outputFile=="" && _name!="") {
        of << _name << ".png";
        this->SetOutputFile(of.str());
        std::cout << __PRETTY_FUNCTION__ << " [DEBUG] output name = " << of.str() << std::endl;
      }

      tmp.open(_tmpFile.c_str());
      for (int i=0; i<_histNum; i++)
        tmp << _histBounds[i] << "\t" << _histValues[i] << std::endl;
      ef = "everyNth(lab,N) =((int(column(0)) % N == 0) ? stringcolumn(lab) : \"\"); ";
      const std::string title = ( !_title.empty() ) ? _title : _name;
      fprintf(_pipe,"set xtics auto\n");
      fprintf(_pipe,"%s\n",ef.c_str());
      fprintf(_pipe,"plot '%s' using 2:xtic(everyNth(1, %d)) w histeps t \"%s\"\n", _tmpFile.c_str(), (int)(_histNum/10), title.c_str());

      fflush(_pipe);
      delete[] _histBounds;
      return 0;
    }

    void
    Gnuplot::operator<<(const std::string& command)
    {
      fprintf(_pipe,"%s\n",command.c_str());
      fflush(_pipe);
      _isPlottable = true; //FIXME need to think about that...
      // flush is necessary, nothing gets plotted else
    }
  }
}

