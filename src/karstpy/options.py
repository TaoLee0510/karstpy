"""Context-local options used by the sparse joint fit."""

from contextlib import contextmanager
from contextvars import ContextVar


_DEFAULTS = {
    "convergence_restarts": False,
    "choose_accepted": False,
    "sparse_joint_template": None,
    "pair_history": None,
    "state_space_drift": None,
}
_overrides = ContextVar("karstpy_options", default={})


def get_option(name):
    """Return the option in effect for the current context."""
    return _overrides.get().get(name, _DEFAULTS[name])


@contextmanager
def option_context(**options):
    """Temporarily override one or more fitting options."""
    unknown = options.keys() - _DEFAULTS.keys()
    if unknown:
        raise KeyError(f"Unknown karstpy options: {', '.join(sorted(unknown))}")
    token = _overrides.set(_overrides.get() | options)
    try:
        yield
    finally:
        _overrides.reset(token)
