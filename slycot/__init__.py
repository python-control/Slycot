try:
    __SLYCOT_SETUP__
except NameError:
    __SLYCOT_SETUP__ = False


if __SLYCOT_SETUP__:
    import sys as _sys
    _sys.stderr.write('Running from Slycot source directory.\n')
    del _sys
else:

    # import slycot.examples

    # The Slycot library is organised by 11-chapters. Each chapter can be identified by a single letter.
    # The following chapters are included:
    # A : Analysis Routines (included)
    # B : Benchmark
    # C : Adaptive Control
    # D : Data Analysis
    # F : Filtering
    # I : Identification
    # M : Mathematical Routines (included)
    # N : Nonlinear Systems
    # S : Synthesis Routines (included)
    # T : Transformation Routines (included)
    # U : Utility Routines


    # Analysis routines (18/61 wrapped)
    from .analysis import (ab01nd,
                           ab04md,
                           ab05md, ab05nd,
                           ab07nd,
                           ab08nd, ab08nz,
                           ab09ad, ab09ax, ab09bd, ab09md, ab09nd,
                           ab13bd, ab13dd, ab13ed, ab13fd, ab13md,
                           ag08bd)
    
    # Benchmark routines (0/6 wrapped)

    # Adaptive control routines (0/0 wrapped)

    # Data analysis routines (0/8 wrapped)

    # Filtering routines (0/6 wrapped)

    # Identification routines (0/15 wrapped)

    # Mathematical routines (8/291 wrapped)
    from .math import (mb02ed, mb03rd, mb03vd, mb03vy, mb03wd,
                       mb05md, mb05nd,
                       mc01td)

    # Nonlinear Systems (0/16 wrapped)

    # Synthesis routines ((17+1)/131 wrapped), sb03md57 is not part of slicot
    from .synthesis import (sb01bd,
                            sb02md, sb02mt, sb02od,
                            sb03md, sb03md57, sb03od,
                            sb04md, sb04qd,
                            sb10ad, sb10dd, sb10fd, sb10hd, sb10jd, sb10yd,
                            sg02ad,
                            sg03ad, sg03bd)
                            
    # Transformation routines (12/77 wrapped)
    from .transform import (tb01id, tb01pd,
                            tb03ad,
                            tb04ad,
                            tb05ad,
                            tc01od,
                            tc04ad,
                            td04ad,
                            tf01md, tf01rd,
                            tg01ad, tg01fd)

    # Utility routines (0/7 wrapped)


    from .version import __version__

    __all__ = [
        ab01nd, ab04md, ab05md, ab05nd, ab07nd, ab08nd, ab08nz,
        ab09ad, ab09ax, ab09bd, ab09md, ab09nd, ab13bd, ab13dd,
        ab13ed, ab13fd, ab13md, ag08bd, mb02ed, mb03rd, mb03vd,
        mb03vy, mb03wd, mb05md, mb05nd, mc01td, sb01bd, sb02md,
        sb02mt, sb02od, sb03md, sb03md57, sb03od, sb04md, sb04qd,
        sb10ad, sb10dd, sb10fd, sb10hd, sb10jd, sb10yd, sg02ad,
        sg03ad, sg03bd, tb01id, tb01pd, tb03ad, tb04ad, tb05ad,
        tc01od, tc04ad, td04ad, tf01md, tf01rd, tg01ad, tg01fd,
        __version__
        ]

def test():
    import pytest
    pytest.main(['--pyargs', 'slycot'])
