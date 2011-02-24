#!/usr/bin/env python
#
#       math.py
#       
#       Copyright 2010 Enrico Avventi <avventi@Lonewolf>
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

def mc01td(dico,dp,p):
    """ dp,stable,nz = mc01td(dico,dp,p)
    
    To determine whether or not a given polynomial P(x) with real
    coefficients is stable, either in the continuous-time or discrete-
    time case.

    A polynomial is said to be stable in the continuous-time case
    if all its zeros lie in the left half-plane, and stable in the
    discrete-time case if all its zeros lie inside the unit circle.

    
    Required arguments:
        dico : input string(len=1)
            Indicates whether the stability test to be applied to P(x) is in 
            the continuous-time or discrete-time case as follows:
            = 'C':  continuous-time case;
            = 'D':  discrete-time case.
        dp : input int
            The degree of the polynomial P(x).  dp >= 0.
        p : input rank-1 array('d') with bounds (dp + 1)
            This array must contain the coefficients of P(x) in increasing 
            powers of x.
    Return objects:
        dp : int
            If P(dp+1) = 0.0 on entry, then dp contains the index of the highest 
            power of x for which P(dp+1) <> 0.0.
        stable : int
            Equal to 1 if P(x) if stable, 0 otherwise.
        nz : int
            The number of unstable zeros.
    """
    hidden = ' (hidden by the wrapper)' 
    arg_list = ['dico', 'dp', 'P', 'stable', 'nz', 'DWORK', 'IWARN'+hidden, 
        'INFO'+hidden]
    out = _wrapper.mc01td(dico,dp,p)
    if out[-1] < 0:
        error_text = "The following argument had an illegal value: "+arg_list[-out[-1]-1]
        e = ValueError(error_text)
        e.info = out[-1]
        raise e
    if out[-1] == 1:
        warings.warn('entry P(x) is the zero polynomial.')
    if out[-1] == 2:
        warings.warn('P(x) may have zeros very close to stability boundary.')
    if out[-2] > 0:
        warnings.warn('The degree of P(x) has been reduced to %i' %(dp-k))
    return out[:-2]

# to be replaced by python wrappers
