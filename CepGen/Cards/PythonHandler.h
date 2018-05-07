#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#ifdef PYTHON
#include <Python.h>

#include "Handler.h"

namespace CepGen
{
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
        static constexpr const char* module_name_ = "mod_name";

        static void throwPythonError( const std::string& message );
        static std::string getPythonPath( const char* file );
        static bool isInteger( PyObject* obj );
        static int asInteger( PyObject* obj );

        void getLimits( PyObject* obj, const char* key, Limits& lim );
        void getParameter( PyObject* parent, const char* key, bool& out );
        void getParameter( PyObject* parent, const char* key, int& out );
        void getParameter( PyObject* parent, const char* key, unsigned long& out );
        void getParameter( PyObject* parent, const char* key, unsigned int& out );
        void getParameter( PyObject* parent, const char* key, double& out );
        void getParameter( PyObject* parent, const char* key, std::string& out );
        void getParameter( PyObject* parent, const char* key, std::vector<std::string>& out );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseParticlesCuts( PyObject* );
        void parseLogging( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseTamingFunctions( PyObject* );
        void parseHadroniser( PyObject* );

#ifdef PYTHON2
        char* sfilename_;
#else
        wchar_t* sfilename_;
#endif
        PyObject* cfg_;
    };
  }
}

#endif

#endif
