from os import listdir, path
from re import compile

from ._hammer_config import _use_root as _HAMMER_USE_ROOT, _install_prefix as _HAMMER_INSTALL_PREFIX  # type: ignore[import-not-found]

if _HAMMER_USE_ROOT:
    try:
        import ROOT  # type: ignore[import-not-found]
    except ImportError:
        raise ImportError(
            "This build of hammer requires PyROOT. "
            "Install ROOT from https://root.cern or via conda: "
            "'conda install -c conda-forge root'"
        ) from None
else:
    try:
        import cppyy  # type: ignore[import-not-found]
    except ImportError:
        raise ImportError(
            "This build of hammer requires standalone cppyy. "
            "Install it with: 'pip install cppyy'"
        ) from None


def _get_library_name() -> str:
    lib_pattern = compile(r".*\.(so|dylib|DLL)$")
    dir_name = path.dirname(path.abspath(__file__))
    for filename in listdir(dir_name):
        if lib_pattern.match(filename):
            return path.join(dir_name, filename)
    raise RuntimeError(
        "Cannot find the Hammer shared library (.so/.dylib/.DLL) in the package directory. "
        "The package may be incomplete or incorrectly installed."
    )


__all__ = [
    "BinContents", "FourMomentum", "Hammer", "HistoInfo",
    "IOBuffer", "IOBuffers", "Log", "PAction", "PID",
    "Particle", "Process", "RecordType", "RootIOBuffer",
    "WTerm",
]

_lib = _get_library_name()

if _HAMMER_USE_ROOT:
    try:
        _include_dir = path.join(_HAMMER_INSTALL_PREFIX, "include")
        if path.isdir(_include_dir):
            ROOT.gInterpreter.AddIncludePath(_include_dir)  # type: ignore[possibly-unbound]
        ROOT.gSystem.Load(_lib)  # type: ignore[possibly-unbound]
        from ROOT.Hammer import (  # type: ignore[import-not-found]
            PID, BinContents, FourMomentum, Hammer,
            HistoInfo, IOBuffer, IOBuffers, Log, PAction,
            Particle, Process, RecordType, RootIOBuffer,
            WTerm)
    except Exception as e:
        raise RuntimeError(f"Failed to load Hammer ROOT bindings: {e}") from e
    from . import root_pythonizer  # noqa: F401
else:
    try:
        _include_dir = path.join(_HAMMER_INSTALL_PREFIX, "include")
        if path.isdir(_include_dir):
            cppyy.add_include_path(_include_dir)  # type: ignore[possibly-unbound]
        cppyy.load_reflection_info(_lib)  # type: ignore[possibly-unbound]
        _Hammer_ns = cppyy.gbl.Hammer  # type: ignore[possibly-unbound]
        PID = _Hammer_ns.PID
        BinContents = _Hammer_ns.BinContents
        FourMomentum = _Hammer_ns.FourMomentum
        Hammer = _Hammer_ns.Hammer
        HistoInfo = _Hammer_ns.HistoInfo
        IOBuffer = _Hammer_ns.IOBuffer
        IOBuffers = _Hammer_ns.IOBuffers
        Log = _Hammer_ns.Log
        Particle = _Hammer_ns.Particle
        Process = _Hammer_ns.Process
    except Exception as e:
        raise RuntimeError(f"Failed to load Hammer cppyy bindings: {e}") from e
    from .cppyy_pythonizer import RecordType, WTerm, PAction  # noqa: F401 (also registers pythonizations)
