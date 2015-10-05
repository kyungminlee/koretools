#pragma once

#include "array.h"

namespace kore {
namespace array {

template <std::size_t Rank>
std::ostream &
operator<<(std::ostream & os, ArrayStructure<Rank> const & as)
{
  using SizeType = typename ArrayStructure<Rank>::SizeType;
  os << "ArrayStructure<" << Rank << ">[";
  os << "shape: {" << as.shape_[0];
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
std::ostream &
operator<< (std::ostream& os, Array<ValueType, Rank> const & arr)
{
  os << "Array<" << typeid(ValueType).name();
  if (std::is_const<ValueType>::value) {
    os << " const ";
  }
  os << "," << Rank << ">";
  os << "[";
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
} // namespace kore
