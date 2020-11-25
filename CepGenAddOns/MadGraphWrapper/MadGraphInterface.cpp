#ifndef MADGRAPH_BIN
#error "*** MADGRAPH_BIN variable not set! ***"
#endif
#ifndef MADGRAPH_PROC_TMPL
#error "*** MADGRAPH_PROC_TMPL variable not set! ***"
#endif
#ifndef CC_CFLAGS
#error "*** CC_CFLAGS variable not set! ***"
#endif

#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Core/Exception.h"

#include <fstream>
#include <array>

#if __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
#elif __cplusplus >= 201103L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "*** no support for filesystem! ***"
#endif

namespace cepgen
{
  MadGraphInterface::MadGraphInterface( const ParametersList& params ) :
    proc_( params.get<std::string>( "process" ) ),
    model_( params.get<std::string>( "model" ) ),
    card_path_( params.get<std::string>( "cardPath", "/tmp/cepgen_mg5_input.dat" ) ),
    tmp_dir_( params.get<std::string>( "tmpDir", "/tmp/cepgen_mg_amc" ) ),
    log_filename_( params.get<std::string>( "logFile", "/tmp/cepgen_mg_amc.log" ) )
  {
    if ( proc_.empty() )
      throw CG_FATAL( "MadGraphInterface" )
        << "'process' keyword not set to the parameters!\n"
        << params;
    std::ofstream log( log_filename_, std::ios::trunc ); // clearing the log
  }

  std::string
  MadGraphInterface::run() const
  {
    prepareCard();

    std::ofstream log( log_filename_, std::ios::app ); // clearing the log
    log << "\n\n*** mg5_aMC process generation ***\n\n";

    CG_INFO( "MadGraphInterface:run" )
      << "Running the mg5_aMC process generation.";
    std::string cmd = MADGRAPH_BIN;
    cmd += " -f "+card_path_;
    log << generateProcess();

    log << "\n\n*** mg5_aMC process library compilation ***\n\n";

#ifdef _WIN32
    const std::string lib_path = "CepGenMadGraphProcess.dll";
#else
    const std::string lib_path = "libCepGenMadGraphProcess.so";
#endif

    CG_INFO( "MadGraphInterface:run" )
      << "Preparing the mg5_aMC process library.";
    log << generateLibrary( lib_path );

    return lib_path;
  }

  std::string
  MadGraphInterface::generateProcess() const
  {
    std::string cmd = MADGRAPH_BIN;
    cmd += " -f "+card_path_;
    return runCommand( cmd );
  }

  std::string
  MadGraphInterface::generateLibrary( const std::string& out_lib ) const
  {
    std::vector<std::string> src_files;

    src_files.emplace_back( prepareMadGraphProcess() );

    const fs::path tmp_path( tmp_dir_ );

    //--- find all processes registered
    std::vector<std::string> processes;
    for ( const auto& p : fs::directory_iterator( tmp_path/"SubProcesses" ) )
      if ( p.path().filename().string()[0] == 'P' ) {
        processes.emplace_back( p.path() );
        for ( const auto& f : fs::directory_iterator( p ) )
          if ( f.path().extension() == ".cc" )
            src_files.emplace_back( f.path() );
      }

    if ( processes.size() != 1 )
      throw CG_FATAL( "MadGraphInterface:generateLibrary" )
        << "Currently only single-process cases are supported!";

    //--- find all model libraries
    std::vector<std::string> models;
    for ( const auto& p : fs::directory_iterator( tmp_path/"lib" ) )
      models.emplace_back( p.path() );

    std::string cmd = CC_CFLAGS;
    cmd += " -fPIC -shared";
    cmd += " -Wno-unused-variable -Wno-int-in-bool-context";
    cmd += " -I"+( tmp_path/"src" ).string();
    cmd += " -I"+processes.at( 0 );
    //cmd += " -L"+( tmp_path/"lib" ).string();
    for ( const auto& m : models )
      cmd += " "+m;
    cmd += " "+utils::merge( src_files, " " );
    cmd += " -o "+out_lib;
    return runCommand( cmd );
  }

  void
  MadGraphInterface::prepareCard() const
  {
    std::ofstream card( card_path_ );
    if ( !model_.empty() )
      card << "import model " << model_ << "\n";
    card << "generate " << proc_ << "\n";
    card << "output standalone_cpp " << tmp_dir_ << "\n";
    card << "exit\n";
    card.close();
  }

  std::string
  MadGraphInterface::prepareMadGraphProcess() const
  {
    std::ifstream tmpl_file( MADGRAPH_PROC_TMPL );
    std::string tmpl = std::string(
      std::istreambuf_iterator<char>( tmpl_file ),
      std::istreambuf_iterator<char>() );
    std::string descr = proc_;
    if ( !model_.empty() )
      descr += " (model: "+model_+")";
    utils::replace_all( tmpl, "XXX_PROC_DESCRIPTION_XXX", descr );
    utils::replace_all( tmpl, "XXX_PROC_NAME_XXX", proc_ );

    std::string src_filename = tmp_dir_+"/cepgen_proc_interface.cpp";
    std::ofstream src_file( src_filename );
    src_file << tmpl;
    src_file.close();
    return src_filename;
  }

  std::string
  MadGraphInterface::runCommand( const std::string& cmd ) const
  {
    std::array<char,256> buffer;
    std::unique_ptr<FILE, decltype(&pclose)> pipe( popen( cmd.c_str(), "r" ), pclose );

    CG_DEBUG( "MadGraphInterface:runCommand" )
      << "Running\n\t" << cmd;

    std::string result;
    while ( fgets( buffer.data(), buffer.size(), pipe.get() ) )
      result += buffer.data();
    return result;
  }
}
