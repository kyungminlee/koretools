#pragma once

#include <cassert>
#include <cstddef>
#include <iterator>
#include <array>
#include <vector>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <type_traits>
#include <typeinfo>

#include "../typedefs.h"
#include "../utility/debug.h"
//#include "../debug.h"

namespace kore {
namespace array {

template <std::size_t Rank> class ArrayStructure;
template <typename ValueType, std::size_t Rank> class Array;//, typename AllocatorType = std::allocator<ValueType>> class Array;

template <std::size_t Rank>
std::ostream & operator<<(std::ostream & os, ArrayStructure<Rank> const & as);

template <typename ValueType, std::size_t Rank>
std::ostream & operator<<(std::ostream & os, Array<ValueType, Rank> const & arr);

#if 0
template <typename ValueType, std::size_t Rank> class ConstArray;
template <typename ValueType, std::size_t Rank>
std::ostream & operator<<(std::ostream & os, ConstArray<ValueType, Rank> const & arr);
#endif

//! Structure of array
template <std::size_t _Rank>
class ArrayStructure
{
public:
  using SizeType = std::size_t;
  using DifferenceType = std::ptrdiff_t;

  static const std::size_t Rank = _Rank;

  static_assert(Rank > 0, "Rank should be greater than 1");

  ArrayStructure(ArrayStructure&& as) noexcept
    : length_(as.length_)
    , shape_(as.shape_)    // move not necessary
    , stride_(as.stride_)
  {
  }

  ArrayStructure(const ArrayStructure& as) noexcept
    : length_(as.length_)
    , shape_(as.shape_)
    , stride_(as.stride_)
  {
  }

  //! Construct array structure with given shape.
  template <typename ... Args>
  ArrayStructure(SizeType size, Args ... args) noexcept
  { init_(size, args...); }

  ArrayStructure& operator=(const ArrayStructure & as) noexcept
  {
    length_ = as.length_;
    shape_ = as.shape_;
    stride_ = as.stride_;
    return *this;
  }


  //! Ravel multiindex to index
  //!
  template <typename ...Args>
  SizeType index(SizeType s, Args ... args) const {
    static_assert(Rank == 1 + sizeof...(Args), "Number of arguments should match Rank");
    return ravel_(s, args...);
  }

  //! Ravel multiindex to index
  //!
  template <typename ...Args>
  SizeType ravel(SizeType s, Args ... args) const {
    static_assert(Rank == 1 + sizeof...(Args), "Number of arguments should match Rank");
    return ravel_(s, args...);
  }

  //! Ravel multiindex to index
  //!
  SizeType ravel(const std::array<SizeType, Rank>& mi) const
  {
    SizeType idx = 0;
    for (SizeType d = 0u; d < Rank && d >= 0; ++d) {
      assert(0 <= mi[d] && mi[d] < shape_[d]);
      idx += stride_[d] * mi[d];
    }
    return idx;
  }

  //! Unravel index to multiindex
  //!
  //!
  std::array<SizeType, Rank> unravel(SizeType s) const
  {
    std::array<SizeType, Rank> mi;
    assert(s < size());
    for (SizeType d = Rank - 1; d < Rank && d >= 0; --d) {
      mi[d] = s % shape_[d];
      s /= shape_[d];
    }
    return mi;
  }

  void show() const
  {
    std::cout << "length_ : " << length_ << std::endl;
    std::cout << "shape_ :";
    for (int d = 0; d < Rank; ++d) {
      std::cout << " " << shape_[d];
    }
    std::cout << std::endl;

    std::cout << "stride_ :";
    for (int d = 0; d < Rank; ++d) {
      std::cout << " " << stride_[d];
    }
    std::cout << std::endl;
  }

  bool operator==(const ArrayStructure & structure) {
    return (shape_ == structure.shape_) && (stride_ == structure.stride_);
  }


  SizeType length() const noexcept { return length_; }
  SizeType size()   const noexcept { return length_; }
  std::array<SizeType, Rank> shape()  const noexcept { return shape_; }
  std::array<SizeType, Rank> stride() const noexcept { return stride_; }
  SizeType shape(size_t d)  const { assert(d < Rank); return shape_[d]; }
  SizeType stride(size_t d) const { assert(d < Rank); return stride_[d]; }

  friend std::ostream& operator<< <Rank>(std::ostream & os, ArrayStructure<Rank> const &);

private:
  //! Initialize array structure.
  //! \param size
  //! \param args
  template <typename ... Args>
  void init_(SizeType size, Args ... args) noexcept
  {
    init_(args...);
    static const std::size_t D = Rank - 1 - sizeof...(Args);
    static_assert(D < Rank, "D should be smaller than Rank");
    shape_[D] = size;
    stride_[D] = shape_[D + 1] * stride_[D + 1];
    length_ *= shape_[D];
  }

  //! Initialize array structure.
  //! \param size
  void init_(SizeType size) noexcept
  {
    static const std::size_t D = Rank - 1;
    static_assert(D + 1 == Rank, "last D+1 should match Rank");
    shape_[D] = size;
    stride_[D] = 1;
    length_ = size;
  }

  //! internal function for raveling multiindex to index
  //!
  template<typename ...Args>
  SizeType ravel_(SizeType idx, Args... args) const
  {
    static const size_t D = Rank - 1 - sizeof...(Args);
    static_assert(D < Rank, "D should be smaller than Rank");
    assert(idx < shape_[D]);
    return idx * stride_[D] + ravel_(args...);
  }

  //! internal function for raveling multiindex to index
  //!
  SizeType ravel_(SizeType idx) const
  {
    static const size_t D = Rank - 1;
    static_assert(D + 1 == Rank, "last D+1 should match Rank");
    assert(idx < shape_[D]);
    return idx * stride_[D];
  }

private:
  SizeType length_; //!< Total number of elements.
  std::array<SizeType, Rank> shape_; //!< Shape
  std::array<SizeType, Rank> stride_; //!< Stride
}; // class ArrayStructure



//! Array<ValueType, Rank>
template <typename _ValueType, std::size_t _Rank>
class Array
{
public:
  static const std::size_t Rank = _Rank;
  using ValueType = _ValueType;

  using ConstValueType     = ValueType const;
  using NonConstValueType  = typename std::remove_const<ValueType>::type;
  using ReferenceType      = ValueType       &;
  using ConstReferenceType = ValueType const &;
  using PointerType        = ValueType       *;
  using ConstPointerType   = ValueType const *;
  using ArrayStructureType = ArrayStructure<Rank>;
  using SizeType           = typename ArrayStructure<Rank>::SizeType;

  using ConstArray = Array<ConstValueType, Rank>;
  using NonConstArray = Array<NonConstValueType, Rank>;
  //using AllocatorType = _AllocatorType;

  //! \name copy & move constructors and assignments
  //@{
  Array() = delete;
  //Array(Array const & arr) = delete;

#if 0
  Array(Array & arr) noexcept
    : structure_(arr.structure_)
    , data_(arr.data_)
    , own_(arr.own_)
  {
  }
#endif

  Array(Array && arr) noexcept
    : structure_(std::move(arr.structure_))
    , data_(std::move(arr.data_))
    , own_(arr.own_)
  {
  }

#if 0
  Array& operator=(Array & arr) noexcept
  {
    structure_ = arr.structure_;
    data_ = arr.data_;
    own_ = arr.own_;
    return *this;
  }
#endif

  Array& operator=(Array && arr) noexcept
  {
    structure_ = std::move(arr.structure_);
    data_ = std::move(arr.data_);
    own_ = arr.own_;
    return *this;
  }
  //@}


  //! \name Constructors with shape
  //@{
  template <typename ... Args>
  explicit Array(SizeType s, Args ... args)
    : structure_(s, args...)
    , data_(new ValueType[structure_.size()], // exception can only occur here.
            std::default_delete<ValueType[]>()) 
    , own_(true)
  {
  }

  explicit Array(const ArrayStructureType& structure)
    : structure_(structure)
    , data_(new ValueType[structure_.size()], // exception can only occur here.
            std::default_delete<ValueType[]>())
    , own_(true)
  {
  }
  //@}

  //! \name Constructors with pointer and shape
  //@{
  
  //!
  //!
  //!
  template<typename ... Args>
  explicit Array(std::shared_ptr<ValueType> data,
                 bool own, 
                 SizeType s, Args ... args) noexcept
    : structure_(s, args...)
    , data_(data)
    , own_(own)
  {
  }

  //!
  //!
  //!
  template<typename ... Args>
  explicit Array(std::shared_ptr<ValueType> data,
                 bool own,
                 ArrayStructureType const & structure) noexcept
    : structure_(structure)
    , data_(data)
    , own_(own)
  {
  }
  

  //! Construct from raw pointer
  //!
  //! Use null deleter
  template <typename ... Args>
  explicit Array(PointerType data, SizeType s, Args ... args) noexcept // MAYBE?
    : structure_(s, args...)
    , data_(data, [](void const *) {})
    , own_(false)
  {
  }
  //@}


  //! Reshape
  //!
  template <typename ... Args>
  Array<ValueType, 1u+sizeof...(Args)>
    reshape(SizeType s, Args ... args) const
  {
    Array<ValueType, 1u+sizeof...(Args)> ret(data_, own_, s, args...);
    KORE_ASSERT(ret.size() == size(), "size of the array should not change during reshape");
    return ret;
  }

  //! Return a view of the array.
  //!
  //! The returned view can potentially temper with data, 
  //! so this should be a non-const member.
  Array view() 
  {
    // TODO: Safety check??
    Array ret(data_, own_, structure_);
    return ret;
  }

  //! Return a constant view of the array.
  //!
  //! The returned view cannot temper with data.
  //! so this should be a const member.
  ConstArray constant() const
  {
    Array<ConstValueType, Rank> ret(data_, own_, structure_);
    return ret;
  }


  //! Return a clone of the array.
  Array<NonConstValueType, Rank> clone() const
  {
    // TODO: Safety check??
    Array<NonConstValueType, Rank> ret(structure_);
    std::copy(cbegin(), cend(), ret.begin());
    return ret;
  }

  // TODO: MAKE SURE THIS IS A GOOD IDEA
  operator ConstArray const &() const
  {
    return *reinterpret_cast<ConstArray const *>(this);
  }

  template <typename ... Args>
  ReferenceType operator()(SizeType midx, Args ... args)
  {
    // TODO: Safety check??
    static_assert(Rank == 1u + sizeof...(Args), "Number of arguments should match Rank");
    auto p = data_.get();
    KORE_ASSERT(p != nullptr, "p should not be a nullptr (memory deallocated)");
    SizeType idx = structure_.index(midx, args...);
    return p[idx];
  }

  template <typename ... Args>
  ConstReferenceType operator()(SizeType midx, Args ... args) const
  {
    static_assert(Rank == 1u + sizeof...(Args), "Number of arguments should match Rank");
    auto p = data_.get();
    KORE_ASSERT(p != nullptr, "p should not be a nullptr (memory deallocated)");
    SizeType idx = structure_.index(midx, args...);
    return p[idx];
  }

  template <typename ... Args>
  ReferenceType at(SizeType midx, Args ... args)
  {
    // TODO: Safety check??
    static_assert(Rank == 1u + sizeof...(Args), "Number of arguments should match Rank");
    auto p = data_.get();
    if (p == nullptr) { throw std::out_of_range("range"); }
    SizeType idx = structure_.index(midx, args...);
    return p[idx];
  }

  template <typename ... Args>
  ConstReferenceType at(SizeType midx, Args ... args) const
  {
    static_assert(Rank == 1u + sizeof...(Args), "Number of arguments should match Rank");
    auto p = data_.get();
    if (p == nullptr) { throw std::out_of_range("range"); }
    SizeType idx = structure_.index(midx, args...);
    return p[idx];
  }

  void fill(ConstReferenceType v)
  {
    // TODO: Safety check??
    std::fill(begin(), end(), v);
  }


  //TODO BOUND CHECKING
  template <typename SrcType>
  Array& operator+=(Array<SrcType, Rank> const & arr)
  {
    if (shape() != arr.shape()) { throw std::length_error("Array::operator+=(). shapes don't match"); }
    assert(shape() == arr.shape());
    auto ip = arr.cbegin();
    auto op = begin();
    for (; op != end(); ++op, ++ip) { *op += *ip; }
    return *this;
  }

  template <typename SrcType>
  Array& operator-=(Array<SrcType, Rank> const & arr)
  {
    if (shape() != arr.shape()) { throw std::length_error("Array::operator-=(). shapes don't match"); }
    assert(shape() == arr.shape());
    auto ip = arr.cbegin();  auto op = begin();
    for (; op != end(); ++op, ++ip) { *op -= *ip; }
    return *this;
  }

  template <typename SrcType>
  Array& operator*=(Array<SrcType, Rank> const & arr)
  {
    if (shape() != arr.shape()) { throw std::length_error("Array::operator*=(). shapes don't match"); }
    assert(shape() == arr.shape());
    auto ip = arr.cbegin();  auto op = begin();
    for (; op != end(); ++op, ++ip) { *op *= *ip; }
    return *this;
  }


  template <typename SrcType>
  Array& operator/=(Array<SrcType, Rank> const & arr)
  {
    if (shape() != arr.shape()) { throw std::length_error("Array::operator/=(). shapes don't match"); }
    assert(shape() == arr.shape());
    auto ip = arr.cbegin();  auto op = begin();
    for (; op != end(); ++op, ++ip) { *op /= *ip; }
    return *this;
  }

  // Use Template to be More Generic
  Array& operator+=(ConstReferenceType v)
  {
    for (auto op = begin(); op != end(); ++op) { *op += v; }
    return *this;
  }

  Array& operator-=(ConstReferenceType v)
  {
    for (auto op = begin(); op != end(); ++op) { *op -= v; }
    return *this;
  }

  Array& operator*=(ConstReferenceType v)
  {
    for (auto op = begin(); op != end(); ++op) { *op *= v; }
    return *this;
  }

  Array& operator/=(ConstReferenceType v)
  {
    for (auto op = begin(); op != end(); ++op) { *op /= v; }
    return *this;
  }


  const ArrayStructureType& structure() const { return structure_; }

  std::shared_ptr<ValueType> data() const { return data_; }

  //! \name Iterator-like stuff
  //@{
  PointerType begin()                { return data_.get(); }
  PointerType end()                  { return data_.get() + length(); }
  ConstPointerType begin()   const { return data_.get(); }
  ConstPointerType end()     const { return data_.get() + length(); }
  ConstPointerType cbegin()    const { return data_.get(); }
  ConstPointerType cend()      const { return data_.get() + length(); }
  //@}
  
  //! \name Getter
  //@{
  bool own()                          const { return own_; }
  SizeType length()                   const { return structure_.length(); }
  SizeType size()                     const { return structure_.size(); }
  std::array<SizeType, Rank> shape()  const { return structure_.shape(); }
  std::array<SizeType, Rank> stride() const { return structure_.stride(); }
  SizeType shape(size_t d)            const { return structure_.shape(d); }
  SizeType stride(size_t d)           const { return structure_.stride(d); }
  //@}
  
  friend std::ostream& operator<< <ValueType, Rank>(std::ostream & os, Array<ValueType, Rank> const &);
private:
  ArrayStructureType structure_;
  std::shared_ptr<ValueType> data_;
  // PointerType start_, finish_;
  bool own_;
}; // class Array



} // namespace array
} // namespace kore

#include "array_impl.h"
