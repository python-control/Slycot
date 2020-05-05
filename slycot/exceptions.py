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

    pass

class SlycotArithmeticError(SlycotError, ArithmeticError):
    """A Slycot computation failed"""

    pass

def filter_docstring_exceptions(docstring):
    """Check a docstring to find exception descriptions"""

    # check-count the message indices
    index = 0
    exdict = {}
    msg = []
    for l in docstring.split('\n'):
        l = l.strip()
        if l[:10] == ":e.info = " and l[-1] == ":":
            try:
                idx = int(l[10:-1])
                if msg:
                    exdict[index] = '\n'.join(msg)
                    msg = []
                index = idx
            except ValueError:
                if msg:
                    exdict[index] = '\n'.join(msg)
                    msg = []
                index = 0
        elif index:
            msg.append(l.strip())
    if msg:
        exdict[index] = '\n'.join(msg)
    return exdict

def raise_if_slycot_error(info, arg_list, docstring=None):
    """Raise exceptions if slycot info returned is non-zero

    For negative info, the argument as indicated in arg_list was erroneous

    For positive info, the matching exception text is recovered from
    the docstring, which may in many cases simply be the python interface
    routine docstring
    """

    if info < 0:
        message = ("The following argument had an illegal value: {}"
                    "".format(arg_list[-info-1]))
        raise SlycotParameterError(message, info)
    elif info > 0:
        # process the docstring for the error message
        messages = filter_docstring_exceptions(docstring)
        try:
            raise SlycotParameterError(messages[info-1], info)
        except IndexError:
            raise SlycotParameterError(
                "Slycot returned an unhandled error code {}".format(info),
                info)

