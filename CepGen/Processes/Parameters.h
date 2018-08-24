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
        template<typename T> T get( std::string key, T def ) const;
        template<typename T> void set( std::string key, const T& value );

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
    template<> int Parameters::get<int>( std::string key, int def ) const;
    template<> void Parameters::set<int>( std::string key, const int& value );
    template<> std::vector<int> Parameters::get<std::vector<int> >( std::string key, std::vector<int> def ) const;
    template<> void Parameters::set<std::vector<int> >( std::string key, const std::vector<int>& value );

    template<> double Parameters::get<double>( std::string key, double def ) const;
    template<> void Parameters::set<double>( std::string key, const double& value );
    template<> std::vector<double> Parameters::get<std::vector<double> >( std::string key, std::vector<double> def ) const;
    template<> void Parameters::set<std::vector<double> >( std::string key, const std::vector<double>& value );

    template<> std::string Parameters::get<std::string>( std::string key, std::string def ) const;
    template<> void Parameters::set<std::string>( std::string key, const std::string& value );
    template<> std::vector<std::string> Parameters::get<std::vector<std::string> >( std::string key, std::vector<std::string> def ) const;
    template<> void Parameters::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value );

    template<> Parameters Parameters::get<Parameters>( std::string key, Parameters def ) const;
    template<> void Parameters::set<Parameters>( std::string key, const Parameters& value );
    template<> std::vector<Parameters> Parameters::get<std::vector<Parameters> >( std::string key, std::vector<Parameters> def ) const;
    template<> void Parameters::set<std::vector<Parameters> >( std::string key, const std::vector<Parameters>& value );
  }
}

#endif
