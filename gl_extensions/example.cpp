#include <graphlab/sdk/toolkit_function_macros.hpp>

using namespace graphlab;
int add_integers(int a, int b) {
      return a + b;
}

//note the use of flexible_type
std::vector<flexible_type> get_value(
        const std::map<flexible_type, flexible_type> &dict, 
        const std::vector<flexible_type> &keys) {
    std::vector<flexible_type> values;
    for( auto key : keys){
        auto iter = dict.find(key);
        if( iter != dict.end()) values.push_back(iter->second);
        else values.push_back(FLEX_UNDEFINED);
    }
    return values;
}

BEGIN_FUNCTION_REGISTRATION
REGISTER_FUNCTION(add_integers, "a", "b"); // provide named parameters
REGISTER_FUNCTION(get_value, "dict", "keys");
END_FUNCTION_REGISTRATION


