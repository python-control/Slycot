# Generate the routines divide by slicot-chapters for sphinx-doc.
# Only prints out the names, copy & past them into slycot_outer.rst and slycot_inner.rst.
import re
import pandas as pd
import warnings

import slycot
slycot.__version__

def get_slycot_routines(sly):
    all_attributes = dir(sly)
    r = re.compile("[a-z][a-z][0-9][0-9a-z][a-z][a-z]")
    matched_attributes = list(filter(r.match, all_attributes))
    return matched_attributes

def get_slycot_routines_help(file,pre=""):
    textfile = open(file, 'r')
    filetext = textfile.read()
    textfile.close()
    lines = filetext.split("\n")
    res = [ele.replace(" ","") for ele in lines]
    res = [ele.replace("_wrapper.","") for ele in res]

    r = re.compile("[a-z][a-z][0-9][0-9a-z][a-z][a-z]")
    matched_attributes = list(filter(r.match, res))
    return matched_attributes

slycot_wrapper = get_slycot_routines(slycot)
slycot_wrapper.sort()
slycot_f2py_wrapper = get_slycot_routines(slycot._wrapper)
slycot_f2py_wrapper.sort()

slycot_wrapper_help = get_slycot_routines_help("source/reference/slycot_outer.rst")
slycot_wrapper_help.sort()

slycot_f2py_wrapper_help = get_slycot_routines_help("source/reference/slycot_inner.rst")
slycot_f2py_wrapper_help.sort()

missing = list(set(slycot_wrapper) - set(slycot_wrapper_help))
if bool(missing):
    warnings.warn(f"The routines {missing} are missing in 'slycot_outer.rst'.")

missing = list(set(slycot_f2py_wrapper) - set(slycot_f2py_wrapper_help))
if bool(missing):
    warnings.warn(f"The routines _wrapper.{missing} are missing in 'slycot_inner.rst'.")