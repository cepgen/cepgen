#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#ifdef PYTHON
#include <Python.h>

#include "Handler.h"

namespace cepgen
{
  namespace strfun { class Parameterisation; }
  class ParametersList;
  namespace card
  {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler : public Handler
    {
      public:
        /// Read a standard configuration card
        explicit PythonHandler( const char* file );
        ~PythonHandler();
        static PyObject* getElement( PyObject* obj, const char* key );
        static PyObject* encode( const char* str );

      private:
        static constexpr const char* MODULE_NAME = "mod_name";
        static constexpr const char* PROCESS_NAME = "process";

        static void throwPythonError( const std::string& message );
        static std::string getPythonPath( const char* file );

        template<typename T> bool is( PyObject* obj ) const;
        template<typename T> T get( PyObject* obj ) const;

        void fillLimits( PyObject* obj, const char* key, Limits& lim );
        void fillParameter( PyObject* parent, const char* key, bool& out );
        void fillParameter( PyObject* parent, const char* key, int& out );
        void fillParameter( PyObject* parent, const char* key, unsigned long& out );
        void fillParameter( PyObject* parent, const char* key, unsigned int& out );
        void fillParameter( PyObject* parent, const char* key, double& out );
        void fillParameter( PyObject* parent, const char* key, std::string& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<int>& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<double>& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<std::string>& out );
        void fillParameter( PyObject* parent, const char* key, ParametersList& out );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseParticlesCuts( PyObject* );
        void parseLogging( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseTamingFunctions( PyObject* );
        void parseHadroniser( PyObject* );
        void parseStructureFunctions( PyObject*, std::shared_ptr<strfun::Parameterisation>& sf_handler );
    };
    template<> bool PythonHandler::is<int>( PyObject* obj ) const;
    template<> int PythonHandler::get<int>( PyObject* obj ) const;
    template<> unsigned long PythonHandler::get<unsigned long>( PyObject* obj ) const;
    template<> bool PythonHandler::is<ParametersList>( PyObject* obj ) const;
    template<> ParametersList PythonHandler::get<ParametersList>( PyObject* obj ) const;
    template<> bool PythonHandler::is<double>( PyObject* obj ) const;
    template<> double PythonHandler::get<double>( PyObject* obj ) const;
    template<> bool PythonHandler::is<std::string>( PyObject* obj ) const;
    template<> std::string PythonHandler::get<std::string>( PyObject* obj ) const;
  }
}

#endif

#endif
