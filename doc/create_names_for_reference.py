# Generate the routines divide by slicot-chapters for sphinx-doc.
# Only prints out the names, copy & past them into slycot_outer.rst and slycot_inner.rst.
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import slycot
slycot.__version__

def get_slycot_routines(sly):
    all_attributes = dir(sly)
    r = re.compile("[a-z][a-z][0-9][0-9a-z][a-z][a-z]")
    matched_attributes = list(filter(r.match, all_attributes))
    return matched_attributes

slycot_wrapper = get_slycot_routines(slycot)
slycot_wrapper.sort()
slycot_f2py_wrapper = get_slycot_routines(slycot._wrapper)
slycot_f2py_wrapper.sort()

from itertools import groupby
from operator import itemgetter

print(f"\nslycot_wrapper {len(slycot_wrapper)}\n")
for chapter_letter, chapter_routines in groupby(sorted(slycot_wrapper), key=itemgetter(0)):
    print(chapter_letter)
    for routine in chapter_routines:
        print(routine)
    print("\n")

print(f"\nslycot_f2py_wrapper {len(slycot_f2py_wrapper)}\n")
for chapter_letter, chapter_routines in groupby(sorted(slycot_f2py_wrapper), key=itemgetter(0)):
    print(chapter_letter)
    for routine in chapter_routines:
        print("_wrapper."+routine)
    print("\n")
