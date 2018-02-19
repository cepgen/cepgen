#ifndef CepGen_Cards_PythonHandler_h
#define CepGen_Cards_PythonHandler_h

#include "Handler.h"
#include <Python.h>

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
        ~PythonHandler();

        /// Store a configuration into a steering card
        //static void store( const Parameters*, const char* file );

      private:
        void throwPythonError( const char* message, const ExceptionType& type = FatalError ) const;
        static const char* decode( PyObject* obj );
        static PyObject* encode( const char* str );
        PyObject* getElement( PyObject* obj, const char* key );
        void getLimits( PyObject* obj, const char* key, Kinematics::Limits& lim );

        void parseIncomingKinematics( PyObject* );
        void parseOutgoingKinematics( PyObject* );
        void parseParticlesCuts( PyObject* );
        void parseIntegrator( PyObject* );
        void parseGenerator( PyObject* );
        void parseTamingFunctions( PyObject* );
        void parseHadroniser( PyObject* );

        /*static void writeProcess( const Parameters*, PyObject* );
        static void writeIncomingKinematics( const Parameters*, PyObject* );
        static void writeOutgoingKinematics( const Parameters*, PyObject* );
        static void writeTamingFunctions( const Parameters*, PyObject* );
        static void writeIntegrator( const Parameters*, PyObject* );
        static void writeHadroniser( const Parameters*, PyObject* );
        static void writeGenerator( const Parameters*, PyObject* );*/

        PyObject* cfg_;
    };
  }
}

#endif
