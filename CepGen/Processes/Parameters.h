#ifndef CepGen_Processes_Parameters_h
#define CepGen_Processes_Parameters_h

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace CepGen
{
  namespace Process
  {
    class Parameters
    {
      public:
        Parameters() {}
        template<typename T> T get( const char* key, T def ) const;
        template<typename T> void set( const char* key, const T& value );

        friend std::ostream& operator<<( std::ostream& os, const Parameters& );

      private:
        std::map<std::string,Parameters> param_values_;
        std::unordered_map<std::string,int> int_values_;
        std::unordered_map<std::string,double> dbl_values_;
        std::unordered_map<std::string,std::string> str_values_;
        std::unordered_map<std::string,std::vector<Parameters> > vec_param_values_;
        std::unordered_map<std::string,std::vector<int> > vec_int_values_;
        std::unordered_map<std::string,std::vector<double> > vec_dbl_values_;
        std::unordered_map<std::string,std::vector<std::string> > vec_str_values_;
    };
    template<> int Parameters::get<int>( const char* key, int def ) const;
    template<> void Parameters::set<int>( const char* key, const int& value );
    template<> std::vector<int> Parameters::get<std::vector<int> >( const char* key, std::vector<int> def ) const;
    template<> void Parameters::set<std::vector<int> >( const char* key, const std::vector<int>& value );

    template<> double Parameters::get<double>( const char* key, double def ) const;
    template<> void Parameters::set<double>( const char* key, const double& value );
    template<> std::vector<double> Parameters::get<std::vector<double> >( const char* key, std::vector<double> def ) const;
    template<> void Parameters::set<std::vector<double> >( const char* key, const std::vector<double>& value );

    template<> std::string Parameters::get<std::string>( const char* key, std::string def ) const;
    template<> void Parameters::set<std::string>( const char* key, const std::string& value );
    template<> std::vector<std::string> Parameters::get<std::vector<std::string> >( const char* key, std::vector<std::string> def ) const;
    template<> void Parameters::set<std::vector<std::string> >( const char* key, const std::vector<std::string>& value );

    template<> Parameters Parameters::get<Parameters>( const char* key, Parameters def ) const;
    template<> void Parameters::set<Parameters>( const char* key, const Parameters& value );
    template<> std::vector<Parameters> Parameters::get<std::vector<Parameters> >( const char* key, std::vector<Parameters> def ) const;
    template<> void Parameters::set<std::vector<Parameters> >( const char* key, const std::vector<Parameters>& value );
  }
}

#endif
