#ifndef CepGen_Core_ParametersList_h
#define CepGen_Core_ParametersList_h

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace CepGen
{
  class ParametersList
  {
    public:
      ParametersList() {}
      ~ParametersList() {} // required for unique_ptr initialisation!
      template<typename T> T get( std::string key, T def ) const;
      template<typename T> T& operator[]( std::string key );
      template<typename T> void set( std::string key, const T& value );
      ParametersList& operator+=( const ParametersList& oth );

      friend std::ostream& operator<<( std::ostream& os, const ParametersList& );

    private:
      std::map<std::string,ParametersList> param_values_;
      std::unordered_map<std::string,int> int_values_;
      std::unordered_map<std::string,double> dbl_values_;
      std::unordered_map<std::string,std::string> str_values_;
      std::unordered_map<std::string,std::vector<ParametersList> > vec_param_values_;
      std::unordered_map<std::string,std::vector<int> > vec_int_values_;
      std::unordered_map<std::string,std::vector<double> > vec_dbl_values_;
      std::unordered_map<std::string,std::vector<std::string> > vec_str_values_;
  };
  template<> int ParametersList::get<int>( std::string key, int def ) const;
  template<> int& ParametersList::operator[]<int>( std::string key );
  template<> void ParametersList::set<int>( std::string key, const int& value );
  template<> std::vector<int> ParametersList::get<std::vector<int> >( std::string key, std::vector<int> def ) const;
  template<> std::vector<int>& ParametersList::operator[]<std::vector<int> >( std::string key );
  template<> void ParametersList::set<std::vector<int> >( std::string key, const std::vector<int>& value );

  template<> double ParametersList::get<double>( std::string key, double def ) const;
  template<> double& ParametersList::operator[]<double>( std::string key );
  template<> void ParametersList::set<double>( std::string key, const double& value );
  template<> std::vector<double> ParametersList::get<std::vector<double> >( std::string key, std::vector<double> def ) const;
  template<> std::vector<double>& ParametersList::operator[]<std::vector<double> >( std::string key );
  template<> void ParametersList::set<std::vector<double> >( std::string key, const std::vector<double>& value );

  template<> std::string ParametersList::get<std::string>( std::string key, std::string def ) const;
  template<> std::string& ParametersList::operator[]<std::string>( std::string key );
  template<> void ParametersList::set<std::string>( std::string key, const std::string& value );
  template<> std::vector<std::string> ParametersList::get<std::vector<std::string> >( std::string key, std::vector<std::string> def ) const;
  template<> std::vector<std::string>& ParametersList::operator[]<std::vector<std::string> >( std::string key );
  template<> void ParametersList::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value );

  template<> ParametersList ParametersList::get<ParametersList>( std::string key, ParametersList def ) const;
  template<> ParametersList& ParametersList::operator[]<ParametersList>( std::string key );
  template<> void ParametersList::set<ParametersList>( std::string key, const ParametersList& value );
  template<> std::vector<ParametersList> ParametersList::get<std::vector<ParametersList> >( std::string key, std::vector<ParametersList> def ) const;
  template<> std::vector<ParametersList>& ParametersList::operator[]<std::vector<ParametersList> >( std::string key );
  template<> void ParametersList::set<std::vector<ParametersList> >( std::string key, const std::vector<ParametersList>& value );
}

#endif
