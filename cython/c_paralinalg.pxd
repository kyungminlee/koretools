# distutils: language = c++
cimport cython
cimport numpy as c_np


ctypedef fused scalar_t:
    c_np.float32_t
    c_np.float64_t
    c_np.complex64_t
    c_np.complex128_t

ctypedef fused real_t:
    c_np.float32_t
    c_np.float64_t

cdef extern from "paralinalg.h" namespace "kore::paralinalg":
    void inverse "kore::paralinalg::inverse" [S](ssize_t n_item,
                                                 ssize_t n, 
                                                 S* in_mat, 
                                                 S* out_mat) nogil
    void transpose "kore::paralinalg::transpose" [S](ssize_t n_item, 
                                                     ssize_t n1, 
                                                     ssize_t n2, 
                                                     S* in_mat, 
                                                     S* out_mat) nogil
    void dot "kore::paralinalg::dot" [S](ssize_t n_item, 
                                         ssize_t n1,
                                         ssize_t n2,
                                         ssize_t n3, 
                                         S* in_mat1, 
                                         S* in_mat2, 
                                         S* out_mat) nogil
    void eigh "kore::paralinalg::eigh" [S, R](ssize_t n_item,
                                              ssize_t n, 
                                              S* in_mats, 
                                              R* out_eigvals, 
                                              S* out_eigvecs) nogil

