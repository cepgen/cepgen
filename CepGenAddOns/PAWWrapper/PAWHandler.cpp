#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGenAddOns/PAWWrapper/PAWTreeInfo.h"

extern "C" {
#include <cfortran/cfortran.h>
#include <cfortran/hbook.h>
//#include <cfortran/packlib.h>
}

#define PAWC_SIZE 50000
#define PAWC_ALIGNMENT 64

typedef struct {
  float PAW[PAWC_SIZE] __attribute__((aligned(PAWC_ALIGNMENT)));
} PAWC_DEF;
#define PAWC COMMON_BLOCK(PAWC, pawc)
COMMON_BLOCK_DEF(PAWC_DEF, PAWC);
//PAWC_DEF PAWC;

typedef struct {
  int iquest[100] __attribute__((aligned(PAWC_ALIGNMENT)));
} QUEST_DEF;
#define QUEST COMMON_BLOCK(QUEST, quest)
COMMON_BLOCK_DEF(QUEST_DEF, QUEST);
//QUEST_DEF QUEST;

namespace cepgen {
  /**
     * Handler for the storage of events in a PAW/HBOOK format
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Oct 2019
     */
  class PAWHandler : public EventExporter {
  public:
    /// Class constructor
    explicit PAWHandler(const ParametersList&);
    ~PAWHandler();

    void setCrossSection(double, double) override;

    void initialise(const Parameters&) override;
    /// Writer operator
    void operator<<(const Event&) override;

  private:
    const std::string filename_;
    const std::string tree_name_;
    const bool compress_;
    int ntuple_size_;
    int fd_{99}, ev_ntup_{4444}, run_ntup_{5555};
    paw::CepGenEvent cepgen_event_;
    paw::CepGenRun cepgen_run_;
  };

  PAWHandler::PAWHandler(const ParametersList& params)
      : EventExporter(params),
        filename_(steer<std::string>("filename")),
        tree_name_(steer<std::string>("treeName")),
        compress_(steer<bool>("compress")),
        ntuple_size_(steer<int>("ntupleSize")) {
    /*int status;
      std::string chopt{"BSIZE"};
      //hbset_(chopt.data(), ntuple_size_, status, chopt.size());
      if (status != 0)
        throw CG_FATAL("PAWHandler") << "Failed to set the size to " << ntuple_size_ << "!\n\t"
                                     << "HBSET returned " << status << ".";*/
    CG_WARNING("");
    //quest_.iq[9] = 256000;
  }

  PAWHandler::~PAWHandler() {
    int icycle;
    HROUT(0, icycle, (char*)" ");
    HREND((char*)tree_name_.data());
    KUCLOS(1, (char*)" ", 1);
  }

  ParametersDescription PAWHandler::description() {
    auto desc = EventExporter::description();
    desc.add<std::string>("filename", "output.hbook");
    desc.add<std::string>("treeName", "NTUPLE");
    desc.add<bool>("compress", false);
    desc.add<int>("ntupleSize", 65536);
    return desc;
  }

  void PAWHandler::initialise(const Parameters&) {
    static bool kHLIMITcalled = false;
    if (!kHLIMITcalled) {
      //HLIMIT(PAWC_SIZE);
      kHLIMITcalled = true;
    }
    QUEST.iquest[9] = 500000;
    CG_WARNING("") << 10;
    int status;
    char tree_name[500];
    char filename[500];
    tree_name_.copy(tree_name, tree_name_.size());
    filename_.copy(filename, filename_.size());
    char mode[] = "NX";
    HROPEN(fd_, tree_name, filename, mode, ntuple_size_, status);
    CG_WARNING("") << 11;
    if (status != 0)
      throw CG_FATAL("PAWHandler") << "Failed to initialise the file \"" << filename_ << "\"!\n\t"
                                   << "HROPEN returned " << status << ".";
    CG_WARNING("") << 0;
    /*const std::string chopt = " "; // disk-resident ntuple
      hbnt_( ev_ntup_, tree_name_.c_str(), chopt.c_str(), tree_name_.size(), chopt.size() );*/
    /*char tags[1024];
      HBOOKN(ev_ntup_, (char*)tree_name_.data(), 17, (char*)"TOPDIR", 500, tags);
      const std::string evt_content =
          "gen_time:R,tot_time:R,weight:R,"
          "np[0,5000]:I,pt(np):R,eta(np):R,phi(np):R,rapidity(np):R,"
          "E(np):R,m(np):R,charge(np):R,"
          "pdg_id(np):I,parent1(np):I,parent2(np):I,stable(np):I,"
          "role(np):I,status(np):I";
      HBNAME(ev_ntup_, (char*)cepgen_event_.TREE_NAME, (void*)&cepgen_event_, (char*)evt_content.data());
      const std::string run_content = "sqrt_s:R,xsect:R,errxsect:R,num_events:I,litigious_events:I";
      HBNAME(run_ntup_, (char*)cepgen_run_.TREE_NAME, (void*)&cepgen_run_, (char*)run_content.data());*/
    HBNT(ev_ntup_, (char*)tree_name_.data(), (char*)" ");
    HBNAME(ev_ntup_, (char*)"BLOCK", cepgen_event_.gen_time, (char*)"R:R*4");

    CG_WARNING("") << 1;
    HBOOK1(1, (char*)"some random distribution", 1000, -4., 4., 0.);
    CG_WARNING("") << 2;
    for (size_t i = 0; i < 1000; ++i)
      HFILL(1, rand() * 1. / RAND_MAX, 0., 1.);
    HPRINT(1);
  }

  void PAWHandler::operator<<(const Event& ev) {
    CG_WARNING("") << ev;
    cepgen_event_.clear();
    cepgen_event_.gen_time = ev.time_generation;
    cepgen_event_.tot_time = ev.time_total;
    const auto& parts = compress_ ? ev.compress().particles() : ev.particles();
    //--- loop over all particles in event
    for (const auto& part : parts) {
      const auto& mom = part.momentum();
      cepgen_event_.rapidity[cepgen_event_.np] = mom.rapidity();
      cepgen_event_.pt[cepgen_event_.np] = mom.pt();
      cepgen_event_.eta[cepgen_event_.np] = mom.eta();
      cepgen_event_.phi[cepgen_event_.np] = mom.phi();
      cepgen_event_.E[cepgen_event_.np] = part.energy();
      cepgen_event_.m[cepgen_event_.np] = part.mass();
      cepgen_event_.pdg_id[cepgen_event_.np] = part.integerPdgId();
      cepgen_event_.parent1[cepgen_event_.np] = part.mothers().size() > 0 ? *part.mothers().begin() : -1;
      cepgen_event_.parent2[cepgen_event_.np] = part.mothers().size() > 1 ? *part.mothers().rbegin() : -1;
      cepgen_event_.status[cepgen_event_.np] = (int)part.status();
      cepgen_event_.stable[cepgen_event_.np] = (short)part.status() > 0;
      cepgen_event_.charge[cepgen_event_.np] = part.charge();
      cepgen_event_.role[cepgen_event_.np] = part.role();
      cepgen_event_.np++;
    }
    HFNT(ev_ntup_);
  }
}  // namespace cepgen

REGISTER_EXPORTER("paw", PAWHandler)
