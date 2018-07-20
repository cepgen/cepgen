#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#ifdef PYTHON
#include <Python.h>

#include "Handler.h"

namespace CepGen
{
  class StructureFunctions;
  namespace Cards
  {
    /// CepGen Python configuration cards reader/writer
    class PythonHandler : public Handler
    {
      public:
        /// Read a standard configuration card
        explicit PythonHandler( const char* file );
        ~PythonHandler();
        static PyObject* getElement( PyObject* obj, const char* key );
        static std::string decode( PyObject* obj );
        static PyObject* encode( const char* str );

      private:
        static constexpr const char* MODULE_NAME = "mod_name";
        static constexpr const char* PROCESS_NAME = "process";

        static void throwPythonError( const std::string& message );
        static std::string getPythonPath( const char* file );
        static bool isInteger( PyObject* obj );
        static int asInteger( PyObject* obj );

        void fillLimits( PyObject* obj, const char* key, Limits& lim );
        void fillParameter( PyObject* parent, const char* key, bool& out );
        void fillParameter( PyObject* parent, const char* key, int& out );
        void fillParameter( PyObject* parent, const char* key, unsigned long& out );
        void fillParameter( PyObject* parent, const char* key, unsigned int& out );
        void fillParameter( PyObject* parent, const char* key, double& out );
        void fillParameter( PyObject* parent, const char* key, std::string& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<double>& out );
        void fillParameter( PyObject* parent, const char* key, std::vector<std::string>& out );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseParticlesCuts( PyObject* );
        void parseLogging( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseTamingFunctions( PyObject* );
        void parseHadroniser( PyObject* );
        void parseStructureFunctions( PyObject*, std::shared_ptr<StructureFunctions>& sf_handler );
    };
  }
}

#endif

#endif

