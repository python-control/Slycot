""" set-conda-test-matrix.py

Create test matrix for conda packages
"""
import json, re
from pathlib import Path

osmap = {'linux': 'ubuntu',
         'osx': 'macos',
         'win': 'windows',
         }

conda_jobs = []
for conda_pkg_file in Path("slycot-conda-pkgs").glob("*/*.tar.bz2"):
    cos = osmap[conda_pkg_file.parent.name.split("-")[0]]
    m = re.search(r'py(\d)(\d+)_', conda_pkg_file.name)
    pymajor, pyminor = int(m[1]), int(m[2])
    cpy = f'{pymajor}.{pyminor}'
    for cbl in  ['unset', 'Generic', 'OpenBLAS', 'Intel10_64lp']:
        cjob = {'packagekey': f'{cos}-{cpy}',
                'os': cos,
                'python': cpy,
                'blas_lib':  cbl}
        if (cos == 'windows' and
            pyminor < 8 and
            cbl not in ['unset', 'Intel10_64lp']):
            # fatal windows error because numpy and matplotlib directly
            # link to mkl on older versions
            cjob['failok'] = "FAILOK"
        conda_jobs.append(cjob)

matrix = { 'include': conda_jobs }
print(json.dumps(matrix))