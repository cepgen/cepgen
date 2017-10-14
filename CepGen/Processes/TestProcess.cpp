#include "TestProcess.h"

using namespace CepGen::Process;

TestProcess::TestProcess() :
  GenericProcess( "test", ".oO TEST PROCESS Oo.", false ),
  use_default_formula_( true ), num_args_( 3 )
{}

TestProcess::TestProcess( const char* formula, std::vector<std::string> args ) :
  GenericProcess( "test", Form( ".oO TEST PROCESS (%s) Oo.", formula ), false ),
  use_default_formula_( false ), num_args_( args.size() )
{
  if ( num_args_ == 2 ) {
    std::array<std::string,2> args_arr;
    std::copy_n( args.begin(), num_args_, args_arr.begin() );
    funct2_ = Functional<2>( formula, args_arr );
  }
  else if ( num_args_ == 3 ) {
    std::array<std::string,3> args_arr;
    std::copy_n( args.begin(), num_args_, args_arr.begin() );
    funct3_ = Functional<3>( formula, args_arr );
  }
}

double
TestProcess::computeWeight()
{
  if ( use_default_formula_ )
    return 1./( 1.-cos( x( 0 )*M_PI )*cos( x( 1 )*M_PI )*cos( x( 2 )*M_PI ) );
  if ( num_args_ == 2 )
    return funct2_.eval( { { x( 0 ), x( 1 ) } } );
  if ( num_args_ == 3 )
    return funct3_.eval( { { x( 0 ), x( 1 ), x( 2 ) } } );
  return 0.;
}
