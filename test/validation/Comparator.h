/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Validation_Comparator_h
#define CepGen_Validation_Comparator_h

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

namespace cepgen::validation {
  class Comparator : public SteeredObject<Comparator> {
  public:
    explicit Comparator(Generator& gen, const ParametersList& params)
        : SteeredObject(params),
          gen_(gen),
          top_label_(steer<std::string>("topLabel")),
          path_tmpl_(steer<std::string>("pathTemplate")),
          num_events_(steer<int>("numEvents")) {}
    inline ~Comparator() {
      try {
        finalise();
      } catch (const Exception& err) {
        CG_ERROR("Comparator") << "Caught exception while finalising the comparison:\n" << err.what();
      }
    }
    virtual void initialise() = 0;
    Comparator& book(const std::string& name, const std::string& var, const std::string& unit, utils::Hist1D hist) {
      hist.xAxis().setLabel(var + (unit.empty() ? "" : " (" + unit + ")"));
      hist.yAxis().setLabel("d$\\sigma$/d" + var + " (pb" + (unit.empty() ? "" : "/" + unit) + ")");
      m_hist1d_tmpl_.insert(std::make_pair(name, std::move(hist)));
      m_draw_modes_[name] = cepgen::utils::Drawer::Mode::nostack | cepgen::utils::Drawer::Mode::grid;
      return *this;
    }
    inline void loop(const std::string& sample_name) {
      if (!initialised_) {
        initialise();
        initialised_ = true;
      }
      addSample(sample_name);
      weight_ = (double)gen_.computeXsection() / num_events_;
      gen_.generate(num_events_, [&](const Event& evt, size_t) { process(evt); });
    }
    Comparator& fill(const std::string& plot_name, double value) {
      m_hist1ds_.at(plot_name).at(this_sample_).fill(value, weight_);
      return *this;
    }

  protected:
    virtual void process(const Event&) = 0;
    Comparator& addSample(const std::string& sample_name) {
      this_sample_ = sample_name;
      if (std::find(samples_.begin(), samples_.end(), this_sample_) == samples_.end()) {
        samples_.emplace_back(sample_name);
        for (const auto& h_tmpl : m_hist1d_tmpl_)
          m_hist1ds_[h_tmpl.first].insert(std::make_pair(this_sample_, h_tmpl.second));
      }
      if (ref_sample_.empty())
        setReferenceSample(this_sample_);
      return *this;
    }
    Comparator& setReferenceSample(const std::string& sample_name) {
      ref_sample_ = sample_name;
      return *this;
    }
    utils::Drawer::Mode& drawMode(const std::string& plot_name) { return m_draw_modes_[plot_name]; }
    void finalise() {
      auto plotter = steer<ParametersList>("plotter");
      if (plotter.empty())
        return;
      auto plt = cepgen::DrawerFactory::get().build(plotter.set<std::string>("format", "png,pdf"));
      for (auto& plot : m_hist1ds_) {
        cepgen::utils::DrawableColl coll;
        for (auto& gr : plot.second) {
          std::string chi2_info;
          if (gr.first != ref_sample_) {  // do not compute chi^2 test for reference sample
            size_t ndf;
            const auto chi2 = gr.second.chi2test(plot.second.at(ref_sample_), ndf);
            chi2_info = cepgen::utils::format(", $\\chi^{2}$/ndf = %.2g/%zu", chi2, ndf);
          }
          gr.second.setTitle(std::string(gr.first + chi2_info));
          coll.emplace_back(&gr.second);
        }
        plt->draw(coll, path_tmpl_ + plot.first, top_label_, m_draw_modes_[plot.first]);
      }
    }

  private:
    std::map<std::string, utils::Hist1D> m_hist1d_tmpl_;
    std::map<std::string, std::map<std::string, utils::Hist1D> > m_hist1ds_;
    std::map<std::string, utils::Drawer::Mode> m_draw_modes_;
    std::vector<std::string> samples_;
    Generator& gen_;
    bool initialised_{false};
    const std::string top_label_, path_tmpl_;
    const int num_events_;
    std::string ref_sample_, this_sample_;
    double weight_{0.};
  };
}  // namespace cepgen::validation

#endif
