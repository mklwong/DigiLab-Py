from py_digilab.read_model import read_sbml

try:
    import py_digilab.version as __version__
except:
    from warnings import warn as _warn
    _warn('version not found in module. Has py_digilab been installed?')