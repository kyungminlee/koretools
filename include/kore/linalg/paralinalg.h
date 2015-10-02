#ifndef KORE_PARALINALG_H_
#define KORE_PARALINALG_H_

#include <tbb/tbb.h>
#include <Eigen/Eigen>
#include "../typedefs.h"

namespace kore {
namespace paralinalg {

// UnitTranspose
//
//
template <typename Scalar, typename SizeType=size_t>
struct UnitTranspose {
 public:
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixMapType = Eigen::Map<MatrixType>;
  
  explicit UnitTranspose(SizeType n1, SizeType n2, Scalar* in, Scalar* out) :
      _n1(n1), _n2(n2), _in(in), _out(out ? out : in)  { }
  explicit UnitTranspose(const UnitTranspose& s) : _n1(s._n1), _n2(s._n2), _in(s._in), _out(s._out) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      MatrixMapType map_in(&_in[_n1*_n2*i], _n1, _n2);
      MatrixMapType map_out(&_out[_n1*_n2*i], _n2, _n1);
      map_out = map_in.transpose();
    }
  }
 private:
  SizeType _n1, _n2;
  Scalar *_in, *_out;
};

template <typename Scalar, typename SizeType=size_t>
void transpose(SizeType n_item, SizeType n1, SizeType n2, Scalar* in, Scalar* out)
{
  UnitTranspose<Scalar, SizeType> rs(n1, n2, in, out);
  tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
} // transpose


// UnitDot
//
//
template <typename Scalar, typename SizeType=size_t>
struct UnitDot {
 public:
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixMapType = Eigen::Map<MatrixType>;

  explicit UnitDot(SizeType n1, SizeType n2, SizeType n3, Scalar* in1, Scalar* in2, Scalar* out) :
      _n1(n1), _n2(n2), _n3(n3), _in1(in1), _in2(in2), _out(out)  { }
  explicit UnitDot(const UnitDot& s) : _n1(s._n1), _n2(s._n2), _n3(s._n3), _in1(s._in1), _in2(s._in2), _out(s._out) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      MatrixMapType map_in1(&_in1[_n1*_n2*i], _n1, _n2);
      MatrixMapType map_in2(&_in2[_n2*_n3*i], _n2, _n3);
      MatrixMapType map_out(&_out[_n1*_n3*i], _n1, _n3);
      map_out = map_in1 * map_in2;
    }
  }
 private:
  SizeType _n1, _n2, _n3;
  Scalar *_in1, *_in2, *_out;
};

template <typename Scalar, typename SizeType=size_t>
void dot(SizeType n_item, SizeType n1, SizeType n2, SizeType n3, Scalar* in1, Scalar* in2, Scalar* out)
{
  UnitDot<Scalar, SizeType> rs(n1, n2, n3, in1, in2, out);
  tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
} // dot


// UnitInverse
//
//
template <typename Scalar, typename SizeType=size_t>
class UnitInverse {
 public:
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixMapType = Eigen::Map<MatrixType>;

  explicit UnitInverse(SizeType n, Scalar* in, Scalar* out)
      : _n(n), _in(in), _out(out ? out : in)  {  }
  explicit UnitInverse(const UnitInverse& s)
      : _n(s._n), _in(s._in), _out(s._out) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      MatrixMapType map_in(&_in[_n*_n*i], _n, _n);
      MatrixMapType map_out(&_out[_n*_n*i], _n, _n);
      map_out = map_in.inverse();
    }
  }
 private:
  SizeType _n;
  Scalar *_in, *_out;
};


template <typename Scalar, typename SizeType=size_t>
void inverse(SizeType n_item, SizeType n, Scalar* in, Scalar* out)
{
  UnitInverse<Scalar, SizeType> rs(n, in, out);
  tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs); 
} // inv


// UnitEigenSolverWithEigenvector
//
//
template <typename Scalar, typename RealScalar,  typename SizeType=size_t>
struct UnitEigenSolverWithEigenvector {
 public:
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixMapType = Eigen::Map<MatrixType>;
  using RealVectorMapType = Eigen::Map<Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> >;

  explicit UnitEigenSolverWithEigenvector(SizeType n, Scalar* mats, RealScalar* eigvals, Scalar* eigvecs)
      : _n(n), _mats(mats), _eigvals(eigvals), _eigvecs(eigvecs)  {  }
  
  explicit UnitEigenSolverWithEigenvector(const UnitEigenSolverWithEigenvector& s)
      : _n(s._n), _mats(s._mats), _eigvals(s._eigvals), _eigvecs(s._eigvecs) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    Eigen::SelfAdjointEigenSolver<MatrixType> solver(_n);
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      MatrixMapType map(&_mats[_n*_n*i], _n, _n);
      RealVectorMapType eival(&_eigvals[_n*i], _n);
      MatrixMapType eivec(&_eigvecs[_n*_n*i], _n, _n);
      solver.compute(map, Eigen::ComputeEigenvectors);
      eival = solver.eigenvalues();
      eivec = solver.eigenvectors();
    }
  }
 private:
  SizeType _n;
  Scalar *_mats;
  RealScalar *_eigvals;
  Scalar *_eigvecs;
};

template <typename Scalar, typename RealScalar,  typename SizeType=size_t>
struct UnitEigenSolver {
 public:
  using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using MatrixMapType = Eigen::Map<MatrixType>;
  using RealVectorMapType = Eigen::Map<Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> >;

  explicit UnitEigenSolver(SizeType n, Scalar* mats, RealScalar* eigvals)
          : _n(n), _mats(mats), _eigvals(eigvals)  {  }
  explicit UnitEigenSolver(const UnitEigenSolver& s)
      : _n(s._n), _mats(s._mats), _eigvals(s._eigvals) {  }
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    Eigen::SelfAdjointEigenSolver<MatrixType> solver(_n);
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      MatrixMapType map(&_mats[_n*_n*i], _n, _n);
      RealVectorMapType eival(&_eigvals[_n*i], _n);
      solver.compute(map);
      eival = solver.eigenvalues();
    }
  }
 private:
  SizeType _n;
  Scalar *_mats;
  RealScalar *_eigvals;
};


template <typename Scalar, typename RealScalar,  typename SizeType=size_t>
void eigh(SizeType n_item, SizeType n, Scalar* mats, RealScalar* eigvals, Scalar* eigvecs)
{
  if (eigvecs) {
    UnitEigenSolverWithEigenvector<Scalar, RealScalar, SizeType> rs(n, mats, eigvals, eigvecs);
    tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
  } else {
    UnitEigenSolver<Scalar, RealScalar, SizeType> rs(n, mats, eigvals);
    tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
  }
} // eigh

} // namespace paralinalg
} // namespace kore


#endif // KORE_PARALINALG_H_
