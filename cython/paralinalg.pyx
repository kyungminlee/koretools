# distutils: language = c++
cimport cython
cimport numpy as c_np
import numpy as np

from c_paralinalg cimport inverse as _inverse, transpose as _transpose, dot as _dot, eigh as _eigh
    
ctypedef fused scalar_t:
    c_np.float32_t
    c_np.float64_t
    c_np.complex64_t
    c_np.complex128_t

ctypedef fused real_t:
    c_np.float32_t
    c_np.float64_t

'''
cdef extern from "paralinalg.h" namespace "kore::paralinalg":
    void _inverse "kore::paralinalg::inverse" [S](ssize_t n_item,
                                                  ssize_t n, 
                                                  S* in_mat, 
                                                  S* out_mat) nogil
    void _transpose "kore::paralinalg::transpose" [S](ssize_t n_item, 
                                                      ssize_t n1, 
                                                      ssize_t n2, 
                                                      S* in_mat, 
                                                      S* out_mat) nogil
    void _dot "kore::paralinalg::dot" [S](ssize_t n_item, 
                                          ssize_t n1,
                                          ssize_t n2,
                                          ssize_t n3, 
                                          S* in_mat1, 
                                          S* in_mat2, 
                                          S* out_mat) nogil
    void _eigh "kore::paralinalg::eigh" [S, R](ssize_t n_item,
                                               ssize_t n, 
                                               S* in_mats, 
                                               R* out_eigvals, 
                                               S* out_eigvecs) nogil
'''

def _check_type_eigen(c_np.ndarray in_mats not None,
                      c_np.ndarray out_eigvals not None, 
                      c_np.ndarray out_eigvecs not None):
    ALLOWED_DTYPES = [(np.float32, np.float32, np.float32),
                      (np.float64, np.float64, np.float64),
                      (np.complex64, np.float32, np.complex64),
                      (np.complex128, np.float64, np.complex128),]

    if (in_mats.dtype, out_eigvals.dtype, out_eigvecs.dtype) not in ALLOWED_DTYPES:
        raise TypeError("dtype does not match")

    if in_mats.ndim != out_eigvecs.ndim:
        raise TypeError("dimensions of in_mats and out_eigvecs should match.")
    if in_mats.ndim != out_eigvals.ndim+1:
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

def _check_type_inverse(c_np.ndarray in_mats not None, 
                        c_np.ndarray out_mats not None):
    if in_mats.dtype != out_mats.dtype:
        raise TypeError("dtype of in_mats and out_mats should match.")
    if in_mats.ndim != out_mats.ndim:
        raise TypeError("dimensions of in_mats and out_mats should match.")
    if in_mats.ndim < 2:
        raise TypeError("dimension of in_mats should be no less than two.")
    cdef ssize_t i
    for i in range(in_mats.ndim):
        if in_mats.shape[i] != out_mats.shape[i]:
            raise TypeError("shapes of in_mats and out_mats should match")
    if in_mats.shape[in_mats.ndim-1] != in_mats.shape[in_mats.ndim-2]:
        raise TypeError("last two dimensions of in_mats should be equal.")

def _check_type_transpose(c_np.ndarray in_mats not None,
                          c_np.ndarray out_mats not None):
    cdef ssize_t i
    if in_mats.dtype != out_mats.dtype:
        raise TypeError("dtype of in_mats and out_mats should match.")
    if in_mats.ndim != out_mats.ndim:
        raise TypeError("dimensions of in_mats and out_mats should match.")
    if in_mats.ndim < 2:
        raise TypeError("dimension of in_mats should be no less than two.")
    for i in range(in_mats.ndim-2):
        if in_mats.shape[i] != out_mats.shape[i]:
            raise TypeError("shapes of in_mats and out_mats "
                            "should match until the last two.")
    if in_mats.shape[in_mats.ndim-1] != out_mats.shape[out_mats.ndim-2]:
        raise TypeError("shape of the last dimension of in_mats does not "
                        "match the second to last dimension of out_mats.")
    if in_mats.shape[in_mats.ndim-2] != out_mats.shape[out_mats.ndim-1]:
        raise TypeError("shape of the second to last dimension of in_mats "
                        "does not match the last dimension of out_mats.")

def _check_type_dot(c_np.ndarray in_mats1 not None,
                    c_np.ndarray in_mats2 not None,
                    c_np.ndarray out_mats not None):
    cdef ssize_t i
    if in_mats1.dtype != out_mats.dtype:
        raise TypeError("dtype of in_mats1 and out_mats should match.")
    if in_mats2.dtype != out_mats.dtype:
        raise TypeError("dtype of in_mats2 and out_mats should match.")
    if in_mats1.ndim != out_mats.ndim:
        raise TypeError("dimensions of in_mats1 and out_mats should match.")
    if in_mats2.ndim != out_mats.ndim:
        raise TypeError("dimensions of in_mats2 and out_mats should match.")
    if in_mats1.ndim < 2:
        raise TypeError("dimension of in_mats1 should be no less than two.")
    if in_mats2.ndim < 2:
        raise TypeError("dimension of in_mats2 should be no less than two.")

    for i in range(in_mats1.ndim-2):
        if in_mats1.shape[i] != out_mats.shape[i]:
            raise TypeError("shapes of in_mats1 and out_mats should match "
                            "until the last two.")
        if in_mats2.shape[i] != out_mats.shape[i]:
            raise TypeError("shapes of in_mats2 and out_mats should match "
                            "until the last two.")

    if in_mats1.shape[in_mats1.ndim-2] != out_mats.shape[out_mats.ndim-2]:
        raise TypeError("the second to last dimension of in_mats1 does not "
                        "match the second to last dimension of out_mats.")

    if in_mats1.shape[in_mats1.ndim-1] != in_mats2.shape[in_mats2.ndim-2]:
        raise TypeError("the last dimension of in_mats1 does not match "
                        "the second to last dimension of in_mats2.")
        
    if in_mats2.shape[in_mats2.ndim-1] != out_mats.shape[out_mats.ndim-1]:
        raise TypeError("the last dimension of in_mats2 does not match "
                        "the last dimension of out_mats.")


def eigh(c_np.ndarray in_mats not None, 
         c_np.ndarray out_eigvals not None, 
         c_np.ndarray out_eigvecs not None):
    """Parallel computation of eigenvalues for a list of Hermitian matrices.

    Args:
        in_mats (array_like): input
        out_eigvals (array_like): eigenvalues
        out_eigvecs (array_like): eigenvectors
    """
    _check_type_eigen(in_mats, out_eigvals, out_eigvecs)

    cdef:
        void* in_mats_data = in_mats.data
        void* out_eigvals_data = out_eigvals.data
        void* out_eigvecs_data = out_eigvecs.data

        void (*eigh_f4)(ssize_t, ssize_t,
                        c_np.float32_t*, c_np.float32_t*, c_np.float32_t*) nogil
        void (*eigh_f8)(ssize_t, ssize_t, 
                        c_np.float64_t*, c_np.float64_t*, c_np.float64_t*) nogil
        void (*eigh_c8)(ssize_t, ssize_t, 
                        c_np.complex64_t*, c_np.float32_t*, c_np.complex64_t*) nogil
        void (*eigh_c16)(ssize_t, ssize_t, 
                         c_np.complex128_t*, c_np.float64_t*, c_np.complex128_t*) nogil

    eigh_f4 = _eigh[c_np.float32_t, c_np.float32_t]
    eigh_f8 = _eigh[c_np.float64_t, c_np.float64_t]
    eigh_c8 = _eigh[c_np.complex64_t, c_np.float32_t]
    eigh_c16 = _eigh[c_np.complex128_t, c_np.float64_t]

    cdef:
        ssize_t n = in_mats.shape[in_mats.ndim-1]
        ssize_t n_item = 1
        ssize_t i

    for i in range(in_mats.ndim-2):
        n_item *= in_mats.shape[i]

    if in_mats.dtype == np.float32:
        with nogil:
            eigh_f4(n_item, n,
                    <c_np.float32_t*> in_mats_data,
                    <c_np.float32_t*> out_eigvals_data, 
                    <c_np.float32_t*> out_eigvecs_data)
    elif in_mats.dtype == np.float64:
        with nogil:
            eigh_f8(n_item, n, 
                    <c_np.float64_t*> in_mats_data, 
                    <c_np.float64_t*> out_eigvals_data, 
                    <c_np.float64_t*> out_eigvecs_data)
    elif in_mats.dtype == np.complex64:
        with nogil:
            eigh_c8(n_item, n, 
                    <c_np.complex64_t*> in_mats_data,
                    <c_np.float32_t*> out_eigvals_data,
                    <c_np.complex64_t*> out_eigvecs_data)
    elif in_mats.dtype == np.complex128:
        with nogil:
            eigh_c16(n_item, n,
                     <c_np.complex128_t*> in_mats_data, 
                     <c_np.float64_t*> out_eigvals_data, 
                     <c_np.complex128_t*> out_eigvecs_data)
    else:
        raise TypeError("dtype of mat not supported")

def inverse(c_np.ndarray in_mats not None, c_np.ndarray out_mats):
    """Compute inverses of matrices
    
    Args:
        in_mats (ndarray): input. Dimension should be greater than two.
        out_mats (ndarray): output. Dimensions should match with in_mats.
                            If None, out_mats = in_mats

    Raises:
        TypeError: If shapes or dtypes don't match.

    """
    if out_mats is None:
        out_mats = in_mats

    _check_type_inverse(in_mats, out_mats)

    cdef:
        ssize_t n = in_mats.shape[in_mats.ndim-1]
        ssize_t n_item = 1
        ssize_t i

    for i in range(in_mats.ndim-2):
        n_item *= in_mats.shape[i]

    cdef:
        void* in_mats_data = in_mats.data
        void* out_mats_data = out_mats.data

        void (*inverse_f4)(ssize_t, ssize_t, c_np.float32_t*, c_np.float32_t*) nogil
        void (*inverse_f8)(ssize_t, ssize_t, c_np.float64_t*, c_np.float64_t*) nogil
        void (*inverse_c8)(ssize_t, ssize_t, c_np.complex64_t*, c_np.complex64_t*) nogil
        void (*inverse_c16)(ssize_t, ssize_t, c_np.complex128_t*, c_np.complex128_t*) nogil
    inverse_f4  = _inverse[c_np.float32_t]
    inverse_f8  = _inverse[c_np.float64_t]
    inverse_c8  = _inverse[c_np.complex64_t]
    inverse_c16 = _inverse[c_np.complex128_t]

    if in_mats.dtype == np.float32:
        inverse_f4(n_item, n, <c_np.float32_t*> in_mats_data, <c_np.float32_t*> out_mats_data)
    elif in_mats.dtype == np.float64:
        inverse_f8(n_item, n, <c_np.float64_t*> in_mats_data, <c_np.float64_t*> out_mats_data)
    elif in_mats.dtype == np.complex64:
        inverse_c8(n_item, n, <c_np.complex64_t*> in_mats_data, <c_np.complex64_t*> out_mats_data)
    elif in_mats.dtype == np.complex128:
        inverse_c16(n_item, n, <c_np.complex128_t*> in_mats_data, <c_np.complex128_t*> out_mats_data)
    else:
        raise TypeError("dtype of mat not supported")

def transpose(c_np.ndarray in_mats not None, c_np.ndarray out_mats):
    """Compute transposes of matrices
    
    Args:
        in_mats (ndarray): input. Dimension should be greater than two.
        out_mats (ndarray): output. Dimensions should match with in_mats, with the last two axes swapped.
                            If None, out_mats = in_mats

    Raises:
        TypeError: If shapes or dtypes don't match.

    """
    if out_mats is None:
        out_mats = in_mats

    _check_type_transpose(in_mats, out_mats)

    cdef:
        ssize_t n1 = in_mats.shape[in_mats.ndim-2]
        ssize_t n2 = in_mats.shape[in_mats.ndim-1]
        ssize_t n_item = 1
        ssize_t i

    for i in range(in_mats.ndim-2):
        n_item *= in_mats.shape[i]

    cdef:
        void* in_mats_data = in_mats.data
        void* out_mats_data = out_mats.data

        void (*transpose_f4)(ssize_t, ssize_t, ssize_t, c_np.float32_t*, c_np.float32_t*) nogil
        void (*transpose_f8)(ssize_t, ssize_t, ssize_t, c_np.float64_t*, c_np.float64_t*) nogil
        void (*transpose_c8)(ssize_t, ssize_t, ssize_t, c_np.complex64_t*, c_np.complex64_t*) nogil
        void (*transpose_c16)(ssize_t, ssize_t, ssize_t, c_np.complex128_t*, c_np.complex128_t*) nogil
    transpose_f4  = _transpose[c_np.float32_t]
    transpose_f8  = _transpose[c_np.float64_t]
    transpose_c8  = _transpose[c_np.complex64_t]
    transpose_c16 = _transpose[c_np.complex128_t]

    if in_mats.dtype == np.float32:
        with nogil:
            transpose_f4(n_item, n1, n2, 
                         <c_np.float32_t*> in_mats_data,
                         <c_np.float32_t*> out_mats_data)
    elif in_mats.dtype == np.float64:
        with nogil:
            transpose_f8(n_item, n1, n2,
                         <c_np.float64_t*> in_mats_data,
                         <c_np.float64_t*> out_mats_data)
    elif in_mats.dtype == np.complex64:
        with nogil:
            transpose_c8(n_item, n1, n2,
                         <c_np.complex64_t*> in_mats_data,
                         <c_np.complex64_t*> out_mats_data)
    elif in_mats.dtype == np.complex128:
        with nogil:
            transpose_c16(n_item, n1, n2, 
                          <c_np.complex128_t*> in_mats_data,
                          <c_np.complex128_t*> out_mats_data)
    else:
        raise TypeError("dtype of mat not supported")

def dot(c_np.ndarray in_mats1 not None, 
        c_np.ndarray in_mats2 not None, 
        c_np.ndarray out_mats not None):
    """Compute dot products of matrices
    
    Args:
        in_mats1 (ndarray): input. Dimension should be greater than two.
        in_mats2 (ndarray): input. Dimension should be greater than two.
        out_mats (ndarray): output.

    Raises:
        TypeError: If shapes or dtypes don't match.

    """
    _check_type_dot(in_mats1, in_mats2, out_mats)

    cdef:
        ssize_t n1 = in_mats1.shape[in_mats1.ndim-2]
        ssize_t n2 = in_mats1.shape[in_mats1.ndim-1]
        ssize_t n3 = in_mats2.shape[in_mats2.ndim-1]

        ssize_t n_item = 1
        ssize_t i

    for i in range(in_mats1.ndim-2):
        n_item *= in_mats1.shape[i]

    cdef:
        void* in_mats1_data = in_mats1.data
        void* in_mats2_data = in_mats2.data
        void* out_mats_data = out_mats.data

        void (*dot_f4)(ssize_t, ssize_t, ssize_t, ssize_t, 
                       c_np.float32_t*, c_np.float32_t*, c_np.float32_t*) nogil
        void (*dot_f8)(ssize_t, ssize_t, ssize_t, ssize_t,
                       c_np.float64_t*, c_np.float64_t*, c_np.float64_t*) nogil
        void (*dot_c8)(ssize_t, ssize_t, ssize_t, ssize_t,
                       c_np.complex64_t*, c_np.complex64_t*, c_np.complex64_t*) nogil
        void (*dot_c16)(ssize_t, ssize_t, ssize_t, ssize_t,
                        c_np.complex128_t*, c_np.complex128_t*, c_np.complex128_t*) nogil
    dot_f4  = _dot[c_np.float32_t]
    dot_f8  = _dot[c_np.float64_t]
    dot_c8  = _dot[c_np.complex64_t]
    dot_c16 = _dot[c_np.complex128_t]

    if in_mats1.dtype == np.float32:
        with nogil:
            dot_f4(n_item, n1, n2, n3,
                   <c_np.float32_t*> in_mats1_data,
                   <c_np.float32_t*> in_mats2_data,
                   <c_np.float32_t*> out_mats_data)
    elif in_mats1.dtype == np.float64:
        with nogil:
            dot_f8(n_item, n1, n2, n3,
                   <c_np.float64_t*> in_mats1_data,
                   <c_np.float64_t*> in_mats2_data,
                   <c_np.float64_t*> out_mats_data)
    elif in_mats1.dtype == np.complex64:
        with nogil:
            dot_c8(n_item, n1, n2, n3,
                   <c_np.complex64_t*> in_mats1_data,
                   <c_np.complex64_t*> in_mats2_data,
                   <c_np.complex64_t*> out_mats_data)
    elif in_mats1.dtype == np.complex128:
        with nogil:
            dot_c16(n_item, n1, n2, n3, 
                    <c_np.complex128_t*> in_mats1_data,
                    <c_np.complex128_t*> in_mats2_data,
                    <c_np.complex128_t*> out_mats_data)
    else:
        raise TypeError("dtype of mat not supported")
