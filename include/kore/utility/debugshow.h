#pragma once

#include <iostream>
#include <array>
#include <typeinfo>

#include "typedefs.h"

namespace kore {
namespace debug {

template <typename T>
void debugshow(const char* name, const T& value, std::ostream& os = std::cout, size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(value).name() << "] = " << value << std::endl;
}

template <typename T, size_t R>
void debugshow(const char* name, const std::array<T, R>& value, std::ostream& os = std::cout, size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(value).name() << "]" << std::endl;

  if(0 < R) {
    debugshow("[0]", value[0], os, indent_level+1);
  }
  for (size_t i = 1 ; i < R ; ++i) {
    std::stringstream ss; ss << "[" << i << "]";
    debugshow(ss.str().c_str(), value[i], os, indent_level+1);
  }
}

}
} // namespace kore
