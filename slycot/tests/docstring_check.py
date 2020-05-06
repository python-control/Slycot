"""
docstring_check.py

Copyright 2020 Slycot team

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.
"""
from numpy.testing import assert_raises
from slycot.exceptions import SlycotError, raise_if_slycot_error

def assert_docstring_parse(docstring, erange, checkvars={}):
    """To check that a docstring can be parsed into exceptions
    See also raise_if_slycot_error

    Parameters
    ----------
        docstring: str
            Documentation string with exception definitions
        erange: int or iterable with int
            Error numbers for which the documentation should have
            exception text
        checkvars: dict, optional
            dict of variables for evaluation of <infospec> and formatting the
            exception message
    """

    # if erange is a simple integer, assume a continous range of errors
    try:
        erange = range(1,erange+1)
    except TypeError:
        pass

    for info in erange:
        assert_raises(SlycotError, raise_if_slycot_error, info, [],
                      docstring, checkvars)
