"""Public fitting entry point."""

from .sparse_joint import fit_sparse_joint


def fit_karst(counts, dt=None, beta=None, sparse_joint_control=None, *, labels=None):
    """Fit the sparse joint model, forwarding named controls to fit_sparse_joint."""
    if sparse_joint_control is not None and not isinstance(sparse_joint_control, dict):
        raise ValueError("`sparse_joint_control` must be a dictionary of named arguments.")
    control = dict(sparse_joint_control or {})
    if "counts" in control:
        raise ValueError("Give `counts` only as the first argument.")
    for name, value in (("dt", dt), ("beta", beta), ("labels", labels)):
        if value is not None:
            if name in control:
                raise ValueError(f"Give `{name}` either directly or in `sparse_joint_control`.")
            control[name] = value
    return fit_sparse_joint(counts, **control)
