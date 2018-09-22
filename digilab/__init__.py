from digilab.read_model import read_sbml

try:
    import digilab.version as __version__
except:
    from warnings import warn as _warn
    _warn('version not found in module. Has digilab been installed?')