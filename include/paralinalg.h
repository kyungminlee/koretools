#ifndef _KORE_PARALINALG_H_
#define _KORE_PARALINALG_H_

#include <tbb/tbb.h>
#include <Eigen/Eigen>

namespace kore {
namespace paralinalg {


template <typename Scalar, typename SizeType=size_t>
struct UnitTranspose {
 public:
  explicit UnitTranspose(SizeType n1, SizeType n2, Scalar* in, Scalar* out) :
      _n1(n1), _n2(n2), _in(in), _out(out ? out : in)  { }
  explicit UnitTranspose(const UnitTranspose& s) : _n1(s._n1), _n2(s._n2), _in(s._in), _out(s._out) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          map_in(&_in[_n1*_n2*i], _n1, _n2);
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          map_out(&_out[_n1*_n2*i], _n2, _n1);
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

template <typename Scalar, typename SizeType=size_t>
class UnitInverse {
 public:
  explicit UnitInverse(SizeType n, Scalar* in, Scalar* out)
      : _n(n), _in(in), _out(out ? out : in)  {  }
  explicit UnitInverse(const UnitInverse& s)
      : _n(s._n), _in(s._in), _out(s._out) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          map_in(&_in[_n*_n*i], _n, _n);
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          map_out(&_out[_n*_n*i], _n, _n);
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


template <typename Scalar, typename RealScalar,  typename SizeType=size_t>
struct UnitEigenSolverWithEigenvector {
 public:
  explicit UnitEigenSolverWithEigenvector(SizeType n, Scalar* mats, RealScalar* eigvals, Scalar* eigvecs)
      : _n(n), _mats(mats), _eigvals(eigvals), _eigvecs(eigvecs)  {  }
  
  explicit UnitEigenSolverWithEigenvector(const UnitEigenSolverWithEigenvector& s)
      : _n(s._n), _mats(s._mats), _eigvals(s._eigvals), _eigvecs(s._eigvecs) {  }
  
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
    Eigen::SelfAdjointEigenSolver<MatrixType> solver(_n);
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> map(&_mats[_n*_n*i], _n, _n);
      Eigen::Map<Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> > eival(&_eigvals[_n*i], _n);
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> eivec(&_eigvecs[_n*_n*i], _n, _n);
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
  explicit UnitEigenSolver(SizeType n, Scalar* mats, RealScalar* eigvals)
          : _n(n), _mats(mats), _eigvals(eigvals)  {  }
  explicit UnitEigenSolver(const UnitEigenSolver& s)
      : _n(s._n), _mats(s._mats), _eigvals(s._eigvals) {  }
  void operator()(const tbb::blocked_range<SizeType>& range) const {
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
    Eigen::SelfAdjointEigenSolver<MatrixType> solver(_n);
    for (SizeType i = range.begin(); i < range.end(); ++i) {
      Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> map(&_mats[_n*_n*i], _n, _n);
      Eigen::Map<Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> > eival(&_eigvals[_n*i], _n);
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
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
  if (eigvecs) {
    UnitEigenSolverWithEigenvector<Scalar, RealScalar, SizeType> rs(n, mats, eigvals, eigvecs);
    tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
  } else {
    UnitEigenSolver<Scalar, RealScalar, SizeType> rs(n, mats, eigvals);
    tbb::parallel_for(tbb::blocked_range<SizeType>(0, n_item), rs);
  }
} // eigh
}
}


#endif // _KORE_PARALINALG_H_
