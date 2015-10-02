#pragma once

#include "array.h"
#include <sstream>

namespace kore {
namespace array {

template <std::size_t Rank>
std::ostream & operator<<(std::ostream & os, ArrayStructure<Rank> const & as)
{
  using SizeType = typename ArrayStructure<Rank>::SizeType;
  os << "ArrayStructure<" << Rank << ">[";
  os << "shape :{" << as.shape_[0];
  for (SizeType d = 1; d < Rank; ++d) {
    os << ", " << as.shape_[d];
  }
  os << "}, ";
  os << "stride: {" << as.stride_[0];
  for (SizeType d = 1; d < Rank; ++d) {
    os << ", " << as.stride_[d];
  }
  os << "}, ";
  os << "length: " << as.length_;
  os << "]";
  return os;
}

template <typename ValueType, std::size_t Rank>
std::ostream & operator<< (std::ostream& os, Array<ValueType, Rank> const & arr)
{
  os << "Array<" << typeid(ValueType).name() << ", " << Rank << ">[";
  os << arr.structure_;
  os << ", data(";
  os << arr.data_.use_count();
  os << "): {";
  if (arr.length() > 0) {
    os << (arr.data_.get())[0];
  }
  for (size_t i = 1, n = arr.length(); i < n; ++i) {
    os << ", " << (arr.data_.get())[i];
  }
  os << "}";
  os << "]";
  return os;
}

template <typename ValueType, std::size_t Rank>
std::ostream & operator<< (std::ostream& os, ConstArray<ValueType, Rank> const & arr)
{
  os << "ConstArray<" << typeid(ValueType).name() << ", " << Rank << ">[";
  os << arr.structure_;
  os << ", data(";
  os << arr.data_.use_count();
  os << "): {";
  if (arr.length() > 0) {
    os << (arr.data_.get())[0];
  }
  for (size_t i = 1, n = arr.length(); i < n; ++i) {
    os << ", " << (arr.data_.get())[i];
  }
  os << "}";
  os << "]";
  return os;
}



} // namespace array


namespace debug {
template <size_t Rank>
void debugshow(const char* name, const array::ArrayStructure<Rank>& arr, std::ostream& os, size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(arr).name() << "]" << std::endl;
  debugshow("length", arr.length(), os, indent_level+1);
  debugshow("shape",  arr.shape(), os, indent_level+1);
  debugshow("stride", arr.stride(), os, indent_level+1);
}

template <typename ValueType, size_t Rank>
void debugshow(const char* name, const array::Array<ValueType, Rank>& arr, std::ostream& os, size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(arr).name() << "]" << std::endl;
  debugshow("structure", arr.structure(), os, indent_level+1);
  debugshow("data",  arr.data(), os, indent_level+1);

  for (size_t i = 0 ; i < indent_level+1 ; ++i) { os << "  "; }
  os << "*data" << std::endl;

  for (int i = 0 ; i < arr.length() ; ++i) {
    for (size_t i = 0 ; i < indent_level+2 ; ++i) {
      os << "  ";
    }
    os << arr.data().get()[i] << std::endl;
  }
  debugshow("own", arr.own(), os, indent_level+1);
}


} // namespace debug


} // namespace kore
