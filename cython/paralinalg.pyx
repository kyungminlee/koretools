
# distutils: language = c++

cimport numpy as c_np
import numpy as np
#from libcpp cimport bool as cbool
#from cpython cimport bool as pbool
#from cython.operator cimport dereference as deref, preincrement as inc
#from ctypes import c_size_t

cdef extern from "paralinalg.h":
    void _inverse "kore::paralinalg::inverse" [S](ssize_t n_item, ssize_t n, S* in_mat, S* out_mat) nogil
    void _transpose "kore::paralinalg::transpose" [S](ssize_t n_item, ssize_t n1, ssize_t n2, S* in_mat, S* out_mat) nogil
    void _eigh "kore::paralinalg::eigh" [S,RS](ssize_t n_item, ssize_t n, S* in_mats, RS* out_eigvals, S* out_eigvecs) nogil

def inverse(c_np.ndarray in_mat, c_np.ndarray out_mat):
    if out_mat is None:
        out_mat = in_mat
    if in_mat.dtype != out_mat.dtype:
        raise TypeError("dtype of in_mat and out_mat should match.")
    if in_mat.ndim != out_mat.ndim:
        raise TypeError("dimensions of in_mat and out_mat should match.")
    if in_mat.ndim < 2:
        raise TypeError("dimension of in_mat should be no less than two.")
    cdef ssize_t i
    for i in range(in_mat.ndim):
        if in_mat.shape[i] != out_mat.shape[i]:
            raise TypeError("shapes of in_mat and out_mat should match")
    if in_mat.shape[in_mat.ndim-1] != in_mat.shape[in_mat.ndim-2]:
        raise TypeError("last two dimensions of in_mat should be equal.")

    cdef ssize_t n = in_mat.shape[in_mat.ndim-1]
    cdef ssize_t n_item

    cdef:
        c_np.float32_t* in_data_float32
        c_np.float32_t* out_data_float32
        c_np.float64_t* in_data_float64
        c_np.float64_t* out_data_float64
        #c_np.float128_t* in_data_float128
        #c_np.float128_t* out_data_float128

        c_np.complex64_t* in_data_complex64
        c_np.complex64_t* out_data_complex64
        c_np.complex128_t* in_data_complex128
        c_np.complex128_t* out_data_complex128
        #c_np.complex256_t* in_data_complex256
        #c_np.complex256_t* out_data_complex256

    n_item = 1
    for i in range(in_mat.ndim-2):
        n_item *= in_mat.shape[i]


    if in_mat.dtype == np.float32:
        in_data_float32 = <c_np.float32_t*> in_mat.data
        out_data_float32 = <c_np.float32_t*> out_mat.data
        _inverse[c_np.float32_t](n_item, n, in_data_float32, out_data_float32)
    elif in_mat.dtype == np.float64:
        in_data_float64 = <c_np.float64_t*> in_mat.data
        out_data_float64 = <c_np.float64_t*> out_mat.data
        _inverse[c_np.float64_t](n_item, n, in_data_float64, out_data_float64)
        """
    elif in_mat.dtype == np.float128:
        in_data_float128 = <c_np.float128_t*> in_mat.data
        out_data_float128 = <c_np.float128_t*> out_mat.data
        _inverse[c_np.float128_t](n_item, n, in_data_float128, out_data_float128)
        """
    elif in_mat.dtype == np.complex64:
        in_data_complex64 = <c_np.complex64_t*> in_mat.data
        out_data_complex64 = <c_np.complex64_t*> out_mat.data
        _inverse[c_np.complex64_t](n_item, n, in_data_complex64, out_data_complex64)
    elif in_mat.dtype == np.complex128:
        in_data_complex128 = <c_np.complex128_t*> in_mat.data
        out_data_complex128 = <c_np.complex128_t*> out_mat.data
        _inverse[c_np.complex128_t](n_item, n, in_data_complex128, out_data_complex128)
        """
    elif in_mat.dtype == np.complex256:
        in_data_complex256 = <c_np.complex256_t*> in_mat.data
        out_data_complex256 = <c_np.complex256_t*> out_mat.data
        _inverse[c_np.complex256_t](n_item, n, in_data_complex256, out_data_complex256)
        """
    else:
        raise TypeError("dtype of mat not supported")

def transpose(c_np.ndarray in_mat, c_np.ndarray out_mat):
    cdef ssize_t ndim = in_mat.ndim
    cdef ssize_t i

    if out_mat is None:
        out_mat = in_mat
    if in_mat.dtype != out_mat.dtype:
        raise TypeError("dtype of in_mat and out_mat should match.")
    if in_mat.ndim != out_mat.ndim:
        raise TypeError("dimensions of in_mat and out_mat should match.")
    if in_mat.ndim < 2:
        raise TypeError("dimension of in_mat should be no less than two.")
    for i in range(in_mat.ndim-2):
        if in_mat.shape[i] != out_mat.shape[i]:
            raise TypeError("shapes of in_mat and out_mat should match until the last two.")
    if in_mat.shape[ndim-1] != out_mat.shape[ndim-2]:
        raise TypeError("shape of the last dimension of in_mat does not match the second to last dimension of out_mat.")
    if in_mat.shape[ndim-2] != out_mat.shape[ndim-1]:
        raise TypeError("shape of the second to last dimension of in_mat does not match the last dimension of out_mat.")

    cdef ssize_t n1 = in_mat.shape[in_mat.ndim-2]
    cdef ssize_t n2 = in_mat.shape[in_mat.ndim-1]
    cdef ssize_t n_item = 1
    for i in range(in_mat.ndim-2):
        n_item *= in_mat.shape[i]


    cdef:
        c_np.float32_t* in_data_float32
        c_np.float32_t* out_data_float32
        c_np.float64_t* in_data_float64
        c_np.float64_t* out_data_float64
        #c_np.float128_t* in_data_float128
        #c_np.float128_t* out_data_float128

        c_np.complex64_t* in_data_complex64
        c_np.complex64_t* out_data_complex64
        c_np.complex128_t* in_data_complex128
        c_np.complex128_t* out_data_complex128
        #c_np.complex256_t* in_data_complex256
        #c_np.complex256_t* out_data_complex256


    if in_mat.dtype == np.float32:
        in_data_float32 = <c_np.float32_t*> in_mat.data
        out_data_float32 = <c_np.float32_t*> out_mat.data
        _transpose[c_np.float32_t](n_item, n1, n2, in_data_float32, out_data_float32)
    elif in_mat.dtype == np.float64:
        in_data_float64 = <c_np.float64_t*> in_mat.data
        out_data_float64 = <c_np.float64_t*> out_mat.data
        _transpose[c_np.float64_t](n_item, n1, n2, in_data_float64, out_data_float64)
    elif in_mat.dtype == np.complex64:
        in_data_complex64 = <c_np.complex64_t*> in_mat.data
        out_data_complex64 = <c_np.complex64_t*> out_mat.data
        _transpose[c_np.complex64_t](n_item, n1, n2, in_data_complex64, out_data_complex64)
    elif in_mat.dtype == np.complex128:
        in_data_complex128 = <c_np.complex128_t*> in_mat.data
        out_data_complex128 = <c_np.complex128_t*> out_mat.data
        _transpose[c_np.complex128_t](n_item, n1, n2, in_data_complex128, out_data_complex128)
    else:
        raise TypeError("dtype of mat not supported")