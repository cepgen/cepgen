#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#ifdef PYTHON
#if PY_MAJOR_VERSION < 3
#  define PYTHON2
#endif

namespace cepgen
{
  namespace card
  {
    //------------------------------------------------------------------
    // typed retrieval helpers
    //------------------------------------------------------------------

    template<> bool
    PythonHandler::is<int>( PyObject* obj ) const
    {
#ifdef PYTHON2
      return PyInt_Check( obj );
#else
      return PyLong_Check( obj );
#endif
    }

    template<> bool
    PythonHandler::is<bool>( PyObject* obj ) const
    {
      return PyBool_Check( obj );
    }

    template<> int
    PythonHandler::get<int>( PyObject* obj ) const
    {
      if ( !is<int>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
#ifdef PYTHON2
      return PyInt_AsLong( obj );
#else
      return PyLong_AsLong( obj );
#endif
    }

    template<> unsigned long
    PythonHandler::get<unsigned long>( PyObject* obj ) const
    {
      if ( !is<int>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
#ifdef PYTHON2
      return PyInt_AsUnsignedLongMask( obj );
#else
      if ( !PyLong_Check( obj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, obj->ob_type->tp_name ) );
      return PyLong_AsUnsignedLong( obj );
#endif
    }

    template<> long long
    PythonHandler::get<long long>( PyObject* obj ) const
    {
      if ( !is<int>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
      return PyLong_AsLongLong( obj );
    }

    template<> bool
    PythonHandler::is<double>( PyObject* obj ) const
    {
      return PyFloat_Check( obj );
    }

    template<> double
    PythonHandler::get<double>( PyObject* obj ) const
    {
      if ( !is<double>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
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
      if ( !is<std::string>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
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
    PythonHandler::is<Limits>( PyObject* obj ) const
    {
      if ( !PyTuple_Check( obj ) || PyTuple_Size( obj ) < 1 || PyTuple_Size( obj ) > 2 )
        return false;
      return true;
    }

    template<> Limits
    PythonHandler::get<Limits>( PyObject* obj ) const
    {
      if ( !is<Limits>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
      Limits out;
      double min = get<double>( PyTuple_GetItem( obj, 0 ) );
      out.min() = min;
      if ( PyTuple_Size( obj ) > 1 ) {
        double max = get<double>( PyTuple_GetItem( obj, 1 ) );
        if ( max != -1 )
          out.max() = max;
      }
      return out;
    }

    template<> bool
    PythonHandler::is<ParametersList>( PyObject* obj ) const
    {
      return PyDict_Check( obj );
    }

    template<> ParametersList
    PythonHandler::get<ParametersList>( PyObject* obj ) const
    {
      if ( !is<ParametersList>( obj ) )
        throwPythonError( Form( "Object has invalid type %s", obj->ob_type->tp_name ) );
      ParametersList out;
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
        else if ( is<ParametersList>( pvalue ) )
          out.set<ParametersList>( skey, get<ParametersList>( pvalue ) );
        else if ( is<Limits>( pvalue ) )
          out.set<Limits>( skey, get<Limits>( pvalue ) );
        else if ( PyTuple_Check( pvalue ) || PyList_Check( pvalue ) ) { // vector
          PyObject* pfirst = PyTuple_GetItem( pvalue, 0 );
          PyObject* pit = nullptr;
          const bool tuple = PyTuple_Check( pvalue );
          const Py_ssize_t num_entries = ( tuple )
            ? PyTuple_Size( pvalue )
            : PyList_Size( pvalue );
          if ( is<int>( pfirst ) ) {
            std::vector<int> vec;
            for ( Py_ssize_t i = 0; i < num_entries; ++i ) {
              pit = ( tuple ) ? PyTuple_GetItem( pvalue, i ) : PyList_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector '%s'", skey.c_str() ) );
              vec.emplace_back( get<int>( pit ) );
            }
            out.set<std::vector<int> >( skey, vec );
          }
          else if ( is<double>( pfirst ) ) {
            std::vector<double> vec;
            for ( Py_ssize_t i = 0; i < num_entries; ++i ) {
              pit = ( tuple ) ? PyTuple_GetItem( pvalue, i ) : PyList_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector '%s'", skey.c_str() ) );
              vec.emplace_back( get<double>( pit ) );
            }
            out.set<std::vector<double> >( skey, vec );
          }
          else if ( is<std::string>( pfirst ) ) {
            std::vector<std::string> vec;
            for ( Py_ssize_t i = 0; i < num_entries; ++i ) {
              pit = ( tuple ) ? PyTuple_GetItem( pvalue, i ) : PyList_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector '%s'", skey.c_str() ) );
              vec.emplace_back( get<std::string>( pit ) );
            }
            out.set<std::vector<std::string> >( skey, vec );
          }
          else if ( is<ParametersList>( pfirst ) ) {
            std::vector<ParametersList> vec;
            for ( Py_ssize_t i = 0; i < num_entries; ++i ) {
              pit = ( tuple ) ? PyTuple_GetItem( pvalue, i ) : PyList_GetItem( pvalue, i );
              if ( pit->ob_type != pfirst->ob_type )
                throwPythonError( Form( "Mixed types detected in vector '%s'", skey.c_str() ) );
              vec.emplace_back( get<ParametersList>( pit ) );
            }
            out.set<std::vector<ParametersList> >( skey, vec );
          }
        }
      }
      return out;
    }
  }
}

#endif

