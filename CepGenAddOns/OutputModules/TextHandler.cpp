#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Parameters.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Plotter.h"
#include "CepGen/Version.h"

#include <fstream>

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHandler : public ExportModule {
    public:
      explicit TextHandler(const ParametersList&);
      ~TextHandler();
      static std::string description() { return "Text-based histogramming tool"; }

      void initialise(const Parameters&) override;
      void setCrossSection(double cross_section, double) override { cross_section_ = cross_section; }
      void operator<<(const Event&) override;

    private:
      std::ofstream file_, hist_file_;
      std::string hist_filename_;
      //--- variables definition
      const std::vector<std::string> variables_;
      const bool save_banner_, save_variables_;
      const bool show_hists_, save_hists_;
      const std::string separator_;

      const utils::EventBrowser browser_;

      std::ostringstream oss_vars_;

      double cross_section_;

      //--- kinematic variables
      double sqrts_;
      unsigned long num_evts_;
      struct Hist1DInfo {
        std::string var;
        utils::Hist1D hist;
      };
      std::vector<Hist1DInfo> hists_;
      struct Hist2DInfo {
        std::string var1, var2;
        utils::Hist2D hist;
      };
      std::vector<Hist2DInfo> hists2d_;
    };

    TextHandler::TextHandler(const ParametersList& params)
        : ExportModule(params),
          file_(params.get<std::string>("filename", "output.txt")),
          hist_filename_(params.get<std::string>("histFilename", "output.hists.txt")),
          variables_(params.get<std::vector<std::string>>("variables")),
          save_banner_(params.get<bool>("saveBanner", true)),
          save_variables_(params.get<bool>("saveVariables", true)),
          show_hists_(params.get<bool>("showHistograms", true)),
          save_hists_(params.get<bool>("saveHistograms", false)),
          separator_(params.get<std::string>("separator", "\t")),
          cross_section_(1.),
          sqrts_(0.),
          num_evts_(0ul) {
      //--- first extract list of variables to store in output file
      oss_vars_.clear();
      std::string sep;
      for (const auto& var : variables_)
        oss_vars_ << sep << var, sep = separator_;
      //--- then extract list of variables to be plotted in histogram
      const auto& hist_vars = params.get<ParametersList>("histVariables");
      for (const auto& key : hist_vars.keys()) {
        const auto& vars = utils::split(key, ':');
        if (vars.size() < 1 || vars.size() > 2)
          throw CG_FATAL("TextHandler") << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvar = hist_vars.get<ParametersList>(key);
        if (vars.size() == 1) {  // 1D histogram
          if (hvar.has<std::vector<double>>("xbins"))
            hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(hvar.get<std::vector<double>>("xbins"))});
          else if (hvar.has<Limits>("xrange"))
            hists_.emplace_back(Hist1DInfo{vars.at(0),
                                           utils::Hist1D(hvar.get<int>("nbinsX", hvar.get<int>("nbins", 25)),
                                                         hvar.get<Limits>("xrange", Limits(0., 1.)))});
          else {
            CG_WARNING("TextHandler") << "Neither xrange nor xbins found in parameters for 1D plot of variable \""
                                      << vars.at(0) << "\".";
            continue;
          }
          auto& hist = hists_.rbegin()->hist;
          hist.setLog(hvar.get<bool>("log", false));
          hist.setName(key);
          hist.setXlabel(vars.at(0));
          hist.setYlabel("d(sig)/d" + vars.at(0) + " (pb/bin)");
        } else if (vars.size() == 2) {  // 2D histogram
          if (hvar.has<std::vector<double>>("xbins") && hvar.has<std::vector<double>>("ybins"))
            hists2d_.emplace_back(Hist2DInfo{
                vars.at(0),
                vars.at(1),
                utils::Hist2D(hvar.get<std::vector<double>>("xbins"), hvar.get<std::vector<double>>("ybins"))});
          else if (hvar.has<Limits>("xrange"))
            hists2d_.emplace_back(Hist2DInfo{vars.at(0),
                                             vars.at(1),
                                             utils::Hist2D(hvar.get<int>("nbinsX", hvar.get<int>("nbins", 25)),
                                                           hvar.get<Limits>("xrange", Limits(0., 1.)),
                                                           hvar.get<int>("nbinsY", 50),
                                                           hvar.get<Limits>("yrange", Limits(0., 1.)))});
          else {
            CG_WARNING("TextHandler")
                << "Neither (x/y)range nor (x/y)bins found in parameters for 1D plot of variables \"" << vars << "\".";
            continue;
          }
          auto& hist = hists2d_.rbegin()->hist;
          hist.setName(key);
          hist.setXlabel(vars.at(0));
          hist.setYlabel(vars.at(1));
          hist.setName("d^2(sig)/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)");
          hist.setLog(hvar.get<bool>("log", false));
        }
      }
      if (save_hists_ && !hists_.empty())
        hist_file_.open(hist_filename_);
    }

    TextHandler::~TextHandler() {
      //--- finalisation of the output file
      file_.close();
      //--- histograms printout
      if (!show_hists_ && !save_hists_)
        return;
      for (auto& h_var : hists_) {
        h_var.hist.scale(cross_section_ / (num_evts_ + 1));
        std::ostringstream os;
        h_var.hist.draw(os);
        if (show_hists_)
          CG_INFO("TextHandler") << os.str();
        if (save_hists_)
          hist_file_ << "\n" << os.str() << "\n";
      }
      for (const auto& h_var : hists2d_) {
        std::ostringstream os;
        h_var.hist.draw(os);
        if (show_hists_)
          CG_INFO("TextHandler") << os.str();
        if (save_hists_)
          hist_file_ << "\n" << os.str() << "\n";
      }
      if (save_hists_)
        CG_INFO("TextHandler") << "Saved " << utils::s("histogram", hists_.size()) << " into \"" << hist_filename_
                               << "\".";
    }

    void TextHandler::initialise(const Parameters& params) {
      sqrts_ = params.kinematics.sqrtS();
      num_evts_ = 0ul;
      if (save_banner_)
        file_ << banner(params, "#") << "\n";
      if (save_variables_)
        file_ << "# " << oss_vars_.str() << "\n";
      if (save_hists_ && !hists_.empty())
        hist_file_ << banner(params, "#") << "\n";
    }

    void TextHandler::operator<<(const Event& ev) {
      //--- write down the variables list in the file
      if (!variables_.empty()) {
        std::string sep;
        for (const auto& var : variables_)
          file_ << sep << browser_.get(ev, var), sep = separator_;
        file_ << "\n";
      }
      //--- increment the corresponding histograms
      for (auto& h_var : hists_)
        h_var.hist.fill(browser_.get(ev, h_var.var));
      for (auto& h_var : hists2d_)
        h_var.hist.fill(browser_.get(ev, h_var.var1), browser_.get(ev, h_var.var2));
      ++num_evts_;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("text", TextHandler)
