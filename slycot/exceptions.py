"""
exceptions.py

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


class SlycotError(RuntimeError):
    """Slycot exception base class"""

    def __init__(self, message, info):
        super(SlycotError, self).__init__(message)
        self.info = info


class SlycotParameterError(SlycotError, ValueError):
    """A Slycot input parameter had an illegal value.

    In case of a wrong input value, the SLICOT routines return a negative
    info parameter indicating which parameter was illegal.
    """

    def __init__(self, message, info, arg_list=None):
        if not message:
            message = ("The following argument had an illegal value: {}"
                       "".format(arg_list[-info-1]))
        super(SlycotParameterError, self).__init__(message, info)


class SlycotArithmeticError(SlycotError, ArithmeticError):
    """A Slycot computation failed"""

    pass
