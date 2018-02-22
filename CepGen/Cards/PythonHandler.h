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
        PythonHandler( const char* file );
        ~PythonHandler() {}
        static PyObject* getElement( PyObject* obj, const char* key );
        static const char* decode( PyObject* obj );
        static PyObject* encode( const char* str );

      private:
        static constexpr const char* module_name_ = "mod_name";

        static void throwPythonError( const char* message, const ExceptionType& type = FatalError );
        static std::string getPythonPath( const char* file );
        static bool isInteger( PyObject* obj );
        static int asInteger( PyObject* obj );

        void getLimits( PyObject* obj, const char* key, Kinematics::Limits& lim );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseParticlesCuts( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseTamingFunctions( PyObject* );
        void parseHadroniser( PyObject* );
    };
  }
}

#endif

#endif
