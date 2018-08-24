#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Processes/Parameters.h"
#include "CepGen/Core/utils.h"

#ifdef PYTHON
#if PY_MAJOR_VERSION < 3
#  define PYTHON2
#endif

namespace CepGen
{
  namespace Cards
  {
    //------------------------------------------------------------------
    // typed retrieval helpers
    //------------------------------------------------------------------

    template<> bool
    PythonHandler::is<int>( PyObject* obj ) const
    {
#ifdef PYTHON2
      return ( PyInt_Check( obj ) || PyBool_Check( obj ) );
#else
      return ( PyLong_Check( obj ) || PyBool_Check( obj ) );
#endif
    }

    template<> int
    PythonHandler::get<int>( PyObject* obj ) const
    {
#ifdef PYTHON2
      return PyInt_AsLong( obj );
#else
      return PyLong_AsLong( obj );
#endif
    }

    template<> unsigned long
    PythonHandler::get<unsigned long>( PyObject* obj ) const
    {
#ifdef PYTHON2
      return PyInt_AsUnsignedLongMask( obj );
#else
      if ( !PyLong_Check( obj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, obj->ob_type->tp_name ) );
      return PyLong_AsUnsignedLong( obj );
#endif
    }

    template<> bool
    PythonHandler::is<double>( PyObject* obj ) const
    {
      return PyFloat_Check( obj );
    }

    template<> double
    PythonHandler::get<double>( PyObject* obj ) const
    {
      return PyFloat_AsDouble( obj );
    }

    template<> bool
    PythonHandler::is<std::string>( PyObject* obj ) const
    {
#ifdef PYTHON2
      return PyString_Check( obj );
#else
      return PyUnicode_Check( obj );
#endif
    }

    template<> std::string
    PythonHandler::get<std::string>( PyObject* obj ) const
    {
      std::string out;
#ifdef PYTHON2
      out = PyString_AsString( obj ); // deprecated in python v3+
#else
      PyObject* pstr = PyUnicode_AsEncodedString( obj, "utf-8", "strict" ); // new
      if ( !pstr )
        throwPythonError( "Failed to decode a Python object!" );
      out = PyBytes_AS_STRING( pstr );
      Py_CLEAR( pstr );
#endif
      return out;
    }

    template<> bool
    PythonHandler::is<Process::Parameters>( PyObject* obj ) const
    {
      return PyDict_Check( obj );
    }

    template<> Process::Parameters
    PythonHandler::get<Process::Parameters>( PyObject* obj ) const
    {
      Process::Parameters out;
      PyObject* pkey = nullptr, *pvalue = nullptr;
      Py_ssize_t pos = 0;
      while ( PyDict_Next( obj, &pos, &pkey, &pvalue ) ) {
        const std::string skey = get<std::string>( pkey );
        if ( is<int>( pvalue ) )
          out.set<int>( skey, get<int>( pvalue ) );
        else if ( is<double>( pvalue ) )
          out.set<double>( skey, get<double>( pvalue ) );
        else if ( is<std::string>( pvalue ) )
          out.set<std::string>( skey, get<std::string>( pvalue ) );
        else if ( is<Process::Parameters>( pvalue ) )
          out.set<Process::Parameters>( skey, get<Process::Parameters>( pvalue ) );
        else if ( PyTuple_Check( pvalue ) ) { // vector
          PyObject* pfirst = PyTuple_GetItem( pvalue, 0 );
          if ( is<int>( pfirst ) ) {
            std::vector<int> vec;
            for ( Py_ssize_t i = 0; i < PyTuple_Size( pvalue ); ++i ) {
              PyObject* pit = PyTuple_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector %s", skey ) );
              vec.emplace_back( get<int>( pit ) );
            }
            out.set<std::vector<int> >( skey, vec );
          }
          else if ( is<double>( pfirst ) ) {
            std::vector<double> vec;
            for ( Py_ssize_t i = 0; i < PyTuple_Size( pvalue ); ++i ) {
              PyObject* pit = PyTuple_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector %s", skey ) );
              vec.emplace_back( get<double>( pit ) );
            }
            out.set<std::vector<double> >( skey, vec );
          }
          else if ( is<std::string>( pfirst ) ) {
            std::vector<std::string> vec;
            for ( Py_ssize_t i = 0; i < PyTuple_Size( pvalue ); ++i ) {
              PyObject* pit = PyTuple_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector %s", skey ) );
              vec.emplace_back( get<std::string>( pit ) );
            }
            out.set<std::vector<std::string> >( skey, vec );
          }
          else if ( is<Process::Parameters>( pfirst ) ) {
            std::vector<Process::Parameters> vec;
            for ( Py_ssize_t i = 0; i < PyTuple_Size( pvalue ); ++i ) {
              PyObject* pit = PyTuple_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector %s", skey ) );
              vec.emplace_back( get<Process::Parameters>( pit ) );
            }
            out.set<std::vector<Process::Parameters> >( skey, vec );
          }
        }
      }
      return out;
    }
  }
}

#endif

