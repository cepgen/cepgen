#ifndef CepGen_Physics_Limits_h
#define CepGen_Physics_Limits_h

#include <utility>
#include <iostream>

namespace CepGen
{
  /// Validity interval for a variable
  class Limits : private std::pair<double,double>
  {
    public:
      /// Define lower and upper limits on a quantity
      Limits( double min = kInvalid, double max = kInvalid );
      Limits( const Limits& );

      /// Lower limit to apply on the variable
      double min() const { return first; }
      /// Lower limit to apply on the variable
      double& min() { return first; }
      /// Upper limit to apply on the variable
      double max() const { return second; }
      /// Upper limit to apply on the variable
      double& max() { return second; }
      /// Find the [0,1] value scaled between minimum and maximum
      double x( double v ) const;
      /// Specify the lower and upper limits on the variable
      void in( double low, double up );
      /// Full variable range allowed
      double range() const;
      /// Have a lower limit?
      bool hasMin() const;
      /// Have an upper limit?
      bool hasMax() const;
      /// Check if the value is inside limits' boundaries
      bool passes( double val ) const;
      /// Is there a lower and upper limit?
      bool valid() const;

      /// Human-readable expression of the limits
      friend std::ostream& operator<<( std::ostream&, const Limits& );

    private:
      static constexpr double kInvalid = -999.999;
  };
}

#endif

