#ifndef CepGen_Physics_Limits_h
#define CepGen_Physics_Limits_h

#include <utility>
#include <iosfwd>

namespace cepgen
{
  /// Validity interval for a variable
  class Limits : private std::pair<double,double>
  {
    public:
      /// Define lower and upper limits on a quantity
      Limits( double min = INVALID, double max = INVALID );
      /// Copy constructor
      Limits( const Limits& );

      Limits operator-() const; ///< Invert this limit
      Limits& operator+=( double c ); ///< Add a constant to this limit
      Limits& operator-=( double c ); ///< Subtract a constant to this limit
      Limits& operator*=( double c ); ///< Multiply this limit by a constant
      friend Limits operator+( Limits lim, double c ); ///< Add a constant to a limit
      friend Limits operator-( Limits lim, double c ); ///< Subtract a constant to a limit
      friend Limits operator*( Limits lim, double c ); ///< Multiply a limit by a constant

      /// Lower limit to apply on the variable
      double min() const { return first; }
      /// Lower limit to apply on the variable
      double& min() { return first; }
      /// Upper limit to apply on the variable
      double max() const { return second; }
      /// Upper limit to apply on the variable
      double& max() { return second; }
      /// Export the limits into external variables
      void save( bool& on, double& lmin, double& lmax ) const;
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
      bool contains( double val ) const;
      /// Is there a lower and upper limit?
      bool valid() const;
      /// Raw value of the limits
      const std::pair<double,double> raw() const { return *this; }

      /// Human-readable expression of the limits
      friend std::ostream& operator<<( std::ostream&, const Limits& );

      /// Placeholder for an invalid value in a limit (for single-edged or invalid limits)
      static constexpr double INVALID = -999.999;
  };
}

#endif

