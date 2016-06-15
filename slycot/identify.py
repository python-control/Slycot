#!/usr/bin/env python
#
#       identify.py
#
#       Copyright 2010-2011 Enrico Avventi <avventi@Lonewolf>
#       Copyright 2011 Jerker Nordh <jerker.nordh@control.lth.se>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License version 2 as
#       published by the Free Software Foundation.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


from slycot import _wrapper
import numpy as _np
import warnings

def ib01ad(meth, alg, jobd, batch, conct, ctrl, nobr, m,
                        l, nsmp, u, ldu, y, ldy, n, r, ldr, sv, rcond,
                        tol, iwork, dwork, ldwork, iwarn, info)
    """To preprocess the input-output data for estimating the matrices of a linear time-invariant dynamical system and to find an
    estimate of the system order. The input-output data can,
    optionally, be processed sequentially.

    Required arguments
    ------------------

        meth : {'M', 'N'}
             Specifies the subspace identification method to be used,
             as follows:
             = 'M':  MOESP  algorithm with past inputs and outputs;
             = 'N':  N4SID  algorithm.
        alg  : {'C', 'F', 'Q'}
             Specifies the algorithm for computing the triangular
             factor R, as follows:
             = 'C':  Cholesky algorithm applied to the correlation
                     matrix of the input-output data;
             = 'F':  Fast QR algorithm;
             = 'Q':  QR algorithm applied to the concatenated block
                     Hankel matrices.
        jobd : {'M', 'N'}
             Specifies whether or not the matrices B and D should later
             be computed using the MOESP approach, as follows:
             = 'M':  the matrices B and D should later be computed
                     using the MOESP approach;
             = 'N':  the matrices B and D should not be computed using
                     the MOESP approach.
             This parameter is not relevant for METH = 'N'.
        batch : {'F', 'I', 'L', 'O'}
             Specifies whether or not sequential data processing is to
             be used, and, for sequential processing, whether or not
             the current data block is the first block, an intermediate
             block, or the last block, as follows:
             = 'F':  the first block in sequential data processing;
             = 'I':  an intermediate block in sequential data
                     processing;
             = 'L':  the last block in sequential data processing;
             = 'O':  one block only (non-sequential data processing).
             NOTE that when  100  cycles of sequential data processing
                  are completed for  BATCH = 'I',  a warning is
                  issued, to prevent for an infinite loop.
        conct : {'C', 'N'}
             Specifies whether or not the successive data blocks in
             sequential data processing belong to a single experiment,
             as follows:
             = 'C':  the current data block is a continuation of the
                     previous data block and/or it will be continued
                     by the next data block;
             = 'N':  there is no connection between the current data
                     block and the previous and/or the next ones.
             This parameter is not used if BATCH = 'O'.

        ctrl : {'C', 'N'}
 """
