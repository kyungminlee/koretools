
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
    void _eigh "kore::paralinalg::eigh" [S, R](ssize_t n_item, ssize_t n, S* in_mats, R* out_eigvals, S* out_eigvecs) nogil

def eigh(c_np.ndarray in_mats not None, c_np.ndarray out_eigvals not None, c_np.ndarray out_eigvecs not None):

    ALLOWED_DTYPES = [(np.float32, np.float32, np.float32),
                      (np.float64, np.float64, np.float64),
                      (np.complex64, np.float32, np.complex64),
                      (np.complex128, np.float64, np.complex128),]

    if (in_mats.dtype, out_eigvals.dtype, out_eigvecs.dtype) not in ALLOWED_DTYPES:
        raise TypeError("dtype does not match")

    if in_mats.ndim != out_eigvecs.ndim:
        raise TypeError("dimensions of in_mats and out_eigvecs should match.")
    if in_mats.ndim != out_eigvals.ndim-1:
        raise TypeError("dimensions of in_mats and out_eigvals should match.")

    if in_mats.ndim < 2:
        raise TypeError("dimension of in_mat should be no less than two.")
    cdef ssize_t i
    for i in range(in_mats.ndim):
        if in_mats.shape[i] != out_eigvecs.shape[i]:
            raise TypeError("shapes of in_mats and out_eigvecs should match")
    for i in range(in_mats.ndim-1):
        if in_mats.shape[i] != out_eigvals.shape[i]:
            raise TypeError("shapes of in_mats and out_eigvals should match")
    if in_mats.shape[in_mats.ndim-1] != in_mats.shape[in_mats.ndim-2]:
        raise TypeError("last two dimensions of in_mats should be equal.")

    cdef ssize_t n = in_mats.shape[in_mats.ndim-1]
    cdef ssize_t n_item

    n_item = 1
    for i in range(in_mats.ndim-2):
        n_item *= in_mats.shape[i]

    cdef:
        c_np.float32_t* in_data_float32
        c_np.float32_t* out_eigvals_float32
        c_np.float32_t* out_eigvecs_float32
        c_np.float64_t* in_data_float64
        c_np.float64_t* out_eigvals_float64
        c_np.float64_t* out_eigvecs_float64

        c_np.complex64_t* in_data_complex64
        c_np.complex64_t* out_eigvecs_complex64
        c_np.complex128_t* in_data_complex128
        c_np.complex128_t* out_eigvecs_complex128

    if in_mats.dtype == np.float32:
        in_data_float32     = <c_np.float32_t*> in_mats.data
        out_eigvals_float32 = <c_np.float32_t*> out_eigvals.data
        out_eigvecs_float32 = <c_np.float32_t*> out_eigvecs.data
        _eigh[c_np.float32_t,c_np.float32_t](n_item, n, in_data_float32, out_eigvals_float32, out_eigvecs_float32)
    elif in_mats.dtype == np.float64:
        in_data_float64     = <c_np.float64_t*> in_mats.data
        out_eigvals_float64 = <c_np.float64_t*> out_eigvals.data
        out_eigvecs_float64 = <c_np.float64_t*> out_eigvecs.data
        _eigh[c_np.float64_t,c_np.float64_t](n_item, n, in_data_float64, out_eigvals_float64, out_eigvecs_float64)
    elif in_mats.dtype == np.complex64:
        in_data_complex64     = <c_np.complex64_t*> in_mats.data
        out_eigvals_float32   = <c_np.float32_t*> out_eigvals.data
        out_eigvecs_complex64 = <c_np.complex64_t*> out_eigvecs.data
        _eigh[c_np.complex64_t,c_np.float32_t](n_item, n, in_data_complex64, out_eigvals_float32, out_eigvecs_complex64)
    elif in_mats.dtype == np.complex128:
        in_data_complex128     = <c_np.complex128_t*> in_mats.data
        out_eigvals_float64    = <c_np.float64_t*> out_eigvals.data
        out_eigvecs_complex128 = <c_np.complex128_t*> out_eigvecs.data
        _eigh[c_np.complex128_t,c_np.float64_t](n_item, n, in_data_complex128, out_eigvals_float64, out_eigvecs_complex128)
    else:
        raise TypeError("dtype of mat not supported")



def inverse(c_np.ndarray in_mat not None, c_np.ndarray out_mat=None):
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

    n_item = 1
    for i in range(in_mat.ndim-2):
        n_item *= in_mat.shape[i]

    cdef:
        c_np.float32_t* in_data_float32
        c_np.float32_t* out_data_float32
        c_np.float64_t* in_data_float64
        c_np.float64_t* out_data_float64

        c_np.complex64_t* in_data_complex64
        c_np.complex64_t* out_data_complex64
        c_np.complex128_t* in_data_complex128
        c_np.complex128_t* out_data_complex128

    if in_mat.dtype == np.float32:
        in_data_float32 = <c_np.float32_t*> in_mat.data
        out_data_float32 = <c_np.float32_t*> out_mat.data
        _inverse[c_np.float32_t](n_item, n, in_data_float32, out_data_float32)
    elif in_mat.dtype == np.float64:
        in_data_float64 = <c_np.float64_t*> in_mat.data
        out_data_float64 = <c_np.float64_t*> out_mat.data
        _inverse[c_np.float64_t](n_item, n, in_data_float64, out_data_float64)
    elif in_mat.dtype == np.complex64:
        in_data_complex64 = <c_np.complex64_t*> in_mat.data
        out_data_complex64 = <c_np.complex64_t*> out_mat.data
        _inverse[c_np.complex64_t](n_item, n, in_data_complex64, out_data_complex64)
    elif in_mat.dtype == np.complex128:
        in_data_complex128 = <c_np.complex128_t*> in_mat.data
        out_data_complex128 = <c_np.complex128_t*> out_mat.data
        _inverse[c_np.complex128_t](n_item, n, in_data_complex128, out_data_complex128)
    else:
        raise TypeError("dtype of mat not supported")

def transpose(c_np.ndarray in_mat not None, c_np.ndarray out_mat=None):
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

        c_np.complex64_t* in_data_complex64
        c_np.complex64_t* out_data_complex64
        c_np.complex128_t* in_data_complex128
        c_np.complex128_t* out_data_complex128

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

