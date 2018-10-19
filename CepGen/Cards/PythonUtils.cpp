#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#ifdef PYTHON

#include <string>
#include <algorithm>
#include <frameobject.h>

#if PY_MAJOR_VERSION < 3
#  define PYTHON2
#endif

namespace cepgen
{
  namespace card
  {
    //------------------------------------------------------------------
    // Python API helpers
    //------------------------------------------------------------------

    std::string
    PythonHandler::pythonPath( const char* file )
    {
      std::string s_filename = file;
      s_filename = s_filename.substr( 0, s_filename.find_last_of( "." ) ); // remove the extension
      std::replace( s_filename.begin(), s_filename.end(), '/', '.' ); // replace all '/' by '.'
      return s_filename;
    }

    void
    PythonHandler::throwPythonError( const std::string& message )
    {
      PyObject* ptype = nullptr, *pvalue = nullptr, *ptraceback_obj = nullptr;
      // retrieve error indicator and clear it to handle ourself the error
      PyErr_Fetch( &ptype, &pvalue, &ptraceback_obj );
      PyErr_Clear();
      // ensure the objects retrieved are properly normalised and point to compatible objects
      PyErr_NormalizeException( &ptype, &pvalue, &ptraceback_obj );
      std::ostringstream oss; oss << message;
      if ( ptype != nullptr ) { // we can start the traceback
        oss << "\n\tError: "
#ifdef PYTHON2
            << PyString_AsString( PyObject_Str( pvalue ) ); // deprecated in python v3+
#else
            << PyUnicode_AsUTF8( PyObject_Str( pvalue ) );
#endif
        PyTracebackObject* ptraceback = (PyTracebackObject*)ptraceback_obj;
        std::string tabul = "↪ ";
        if ( ptraceback != nullptr ) {
          while ( ptraceback->tb_next != nullptr ) {
            PyFrameObject* pframe = ptraceback->tb_frame;
            if ( pframe != nullptr ) {
              int line = PyCode_Addr2Line( pframe->f_code, pframe->f_lasti );
#ifdef PYTHON2
              const char* filename = PyString_AsString( pframe->f_code->co_filename );
              const char* funcname = PyString_AsString( pframe->f_code->co_name );
#else
              const char* filename = PyUnicode_AsUTF8( pframe->f_code->co_filename );
              const char* funcname = PyUnicode_AsUTF8( pframe->f_code->co_name );
#endif
              oss << Form( "\n\t%s%s on %s (line %d)", tabul.c_str(), boldify( funcname ).c_str(), filename, line );
            }
            else
              oss << Form( "\n\t%s issue in line %d", tabul.c_str(), ptraceback->tb_lineno );
            tabul = std::string( "  " )+tabul;
            ptraceback = ptraceback->tb_next;
          }
        }
      }
      Py_Finalize();
      throw CG_FATAL( "PythonHandler:error" ) << oss.str();
    }

    PyObject*
    PythonHandler::encode( const char* str )
    {
      PyObject* obj = PyUnicode_FromString( str ); // new
      if ( !obj )
        throwPythonError( Form( "Failed to encode the following string:\n\t%s", str ) );
      return obj;
    }

    PyObject*
    PythonHandler::element( PyObject* obj, const char* key )
    {
      PyObject* pout = nullptr, *nink = encode( key );
      if ( !nink )
        return pout;
      pout = PyDict_GetItem( obj, nink ); // borrowed
      Py_CLEAR( nink );
      if ( pout )
        CG_DEBUG( "PythonHandler:element" )
          << "retrieved " << pout->ob_type->tp_name << " element \"" << key << "\" "
          << "from " << obj->ob_type->tp_name << " object\n\t"
          << "new reference count: " << pout->ob_refcnt;
      else
        CG_DEBUG( "PythonHandler:element" )
          << "did not retrieve a valid element \"" << key << "\"";
      return pout;
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, bool& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = (bool)get<int>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, int& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = get<int>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, unsigned long& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = get<unsigned long>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, unsigned int& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = get<unsigned long>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, double& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = get<double>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::string& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out = get<std::string>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* obj, const char* key, Limits& out )
    {
      PyObject* pobj = element( obj, key ); // borrowed
      if ( pobj )
        out = get<Limits>( pobj );
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<double>& out )
    {
      out.clear();
      PyObject* pobj = element( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i ); // borrowed
        if ( is<double>( pit ) )
          out.emplace_back( get<double>( pit ) );
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<std::string>& out )
    {
      out.clear();
      PyObject* pobj = element( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type %s", key, pobj->ob_type->tp_name ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i ); // borrowed
        if ( is<std::string>( pit ) )
          out.emplace_back( get<std::string>( PyTuple_GetItem( pobj, i ) ) );
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, std::vector<int>& out )
    {
      out.clear();
      PyObject* pobj = element( parent, key ); // borrowed
      if ( !pobj )
        return;
      if ( !PyTuple_Check( pobj ) )
        throwPythonError( Form( "Object \"%s\" has invalid type", key ) );
      for ( Py_ssize_t i = 0; i < PyTuple_Size( pobj ); ++i ) {
        PyObject* pit = PyTuple_GetItem( pobj, i );
        if ( !is<int>( pit ) )
          throwPythonError( Form( "Object %d has invalid type", i ) );
        out.emplace_back( get<int>( pit ) );
      }
    }

    void
    PythonHandler::fillParameter( PyObject* parent, const char* key, ParametersList& out )
    {
      PyObject* pobj = element( parent, key ); // borrowed
      if ( pobj )
        out += get<ParametersList>( pobj );
    }
  }
}

#endif

