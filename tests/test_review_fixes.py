from karstpy import fit_sparse_joint, prepare_karst_data


def test_explicit_beta_overrides_prepared_data_without_mutating_it(counts):
    data = prepare_karst_data(counts, dt=3, beta=0.001)
    fit = fit_sparse_joint(data, beta=0.0002, history=9)
    assert fit.data.beta == 0.0002
    assert data.beta == 0.001
