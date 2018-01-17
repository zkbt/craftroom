__version__ = '0.1'

# specify whether we're calling this from within setup.py
try:
    __CRAFTROOM_SETUP__
except NameError:
    __CRAFTROOM_SETUP__ = False

if not __CRAFTROOM_SETUP__:
    # (run this stuff if it's not form within setup.py)
    pass
