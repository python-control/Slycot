import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import slycot
slycot.__version__

def get_slycot_routines(sly):
    all_attributes = dir(sly)
    r = re.compile("[a-z][a-z][0-9][0-9a-z][a-z][a-z]")
    matched_attributes = list(filter(r.match, all_attributes)) # Read Note below
    return matched_attributes

slycot_wrapper = get_slycot_routines(slycot)
slycot_wrapper.sort()
slycot_f2py_wrapper = get_slycot_routines(slycot._wrapper)
slycot_f2py_wrapper.sort()

print(f"\nslycot_wrapper {len(slycot_wrapper)}\n")
for routine in slycot_wrapper:
    print(routine)

print(f"\nslycot_f2py_wrapper {len(slycot_f2py_wrapper)}\n")
for routine in slycot_f2py_wrapper:
    print("_wrapper."+routine)

