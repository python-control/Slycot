#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from scipy._lib._util import _asarray_validated
from scipy.linalg.misc import LinAlgError, _datacopied
from scipy.linalg.blas import find_best_blas_type
from scipy.linalg.lapack import _compute_lwork

from . import lapack_svd


def svd_full(
    a,
    full_matrices=True,
    compute_u=True,
    compute_v=True,
    overwrite_a=False,
    check_finite=True,
):
    """ """
    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2:
        raise ValueError("expected matrix")
    m, n = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    prefix, dtype, prefer_fortran = find_best_blas_type(arrays=(a,))
    if prefix == "d":
        gesXd, gesXd_lwork = lapack_svd.dgesvd, lapack_svd.dgesvd_lwork
    elif prefix == "s":
        gesXd, gesXd_lwork = lapack_svd.sgesvd, lapack_svd.sgesvd_lwork
    elif prefix == "c":
        gesXd, gesXd_lwork = lapack_svd.cgesvd, lapack_svd.cgesvd_lwork
    else:
        gesXd, gesXd_lwork = lapack_svd.zgesvd, lapack_svd.zgesvd_lwork

    # compute optimal lwork
    lwork = _compute_lwork(
        gesXd_lwork,
        a1,
        compute_u=compute_u,
        compute_v=compute_v,
        full_matrices=full_matrices,
    )

    # perform decomposition
    u, s, v, info = gesXd(
        a1,
        compute_u=compute_u,
        compute_v=compute_v,
        lwork=lwork,
        full_matrices=full_matrices,
        overwrite_a=overwrite_a,
    )

    if info > 0:
        raise LinAlgError("SVD did not converge")
    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal gesdd" % -info)
    if compute_u and compute_v:
        return u, s, v
    elif compute_u:
        return u, s
    elif compute_v:
        return s, v
    else:
        return s


def svd(
    a,
    compute_u=True,
    compute_v=True,
    il=None,
    iu=None,
    overwrite_a=False,
    check_finite=True,
):
    """ """
    # check if the bounds check is being used
    if il is None and iu is None:
        return svd_full(
            a=a,
            compute_u=compute_u,
            compute_v=compute_v,
            overwrite_a=overwrite_a,
            check_finite=check_finite,
        )
    minmn = min(a.shape[0], a.shape[1])

    if il is None:
        il = 1
    elif il < 0:
        if il < -minmn:
            il = -minmn
        il = minmn + il + 1
    else:
        # change to fortran convention
        il += 1
        if il > minmn:
            il = minmn

    if iu is None:
        iu = minmn
    elif iu < 0:
        iu = minmn + iu + 1
    else:
        # change to fortran convention
        iu += 1
        if iu > minmn:
            iu = minmn

    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2:
        raise ValueError("expected matrix")
    m, n = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    prefix, dtype, prefer_fortran = find_best_blas_type(arrays=(a,))
    if prefix == "d":
        gesXd, gesXd_lwork = lapack_svd.dgesvdx, lapack_svd.dgesvdx_lwork
    elif prefix == "s":
        gesXd, gesXd_lwork = lapack_svd.sgesvdx, lapack_svd.sgesvdx_lwork
    elif prefix == "c":
        gesXd, gesXd_lwork = lapack_svd.cgesvdx, lapack_svd.cgesvdx_lwork
    else:
        gesXd, gesXd_lwork = lapack_svd.zgesvdx, lapack_svd.zgesvdx_lwork

    # compute optimal lwork
    lwork = _compute_lwork(
        gesXd_lwork,
        a1,
        compute_u=compute_u,
        compute_v=compute_v,
        il=il,
        iu=iu,
    )

    # perform decomposition
    u, ns, s, v, iwork, info = gesXd(
        a1,
        compute_u=compute_u,
        compute_v=compute_v,
        il=il,
        iu=iu,
        lwork=lwork,
        overwrite_a=overwrite_a,
    )

    if info > 0:
        raise LinAlgError("SVD did not converge")
    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal gesdd" % -info)
    if compute_u and compute_v:
        return u[:ns], s[:ns], v[:ns]
    elif compute_u:
        return u[:ns], s[:ns]
    elif compute_v:
        return s[:ns], v[:ns]
    else:
        return s[:ns]
