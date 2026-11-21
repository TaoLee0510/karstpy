# karst 1.1.3 函数与测试索引

由 R 解析器从 `_reference/karst-1.1.3/` 自动生成（2026-11-20），行号即快照中的行号。

列的含义：

- `导出`：NAMESPACE 中导出，或是 S3 方法。
- `标记`：函数体里用到的东西。TMB = TMB 对象或其内部状态；nlminb；Matrix = Matrix 包；RNG = 随机数；option = `getOption()`；C++ = 调用 Rcpp 函数。
- `调用者`：快照中调用它的其他顶层函数，按函数名加左括号的文本匹配得到，注释里的提及也算在内，所以可能多列；以字符串或 `do.call` 方式的调用会漏掉。

`R/sparse-joint-cohort.R` 不迁移（D10 = b），列在这里只是为了完整。

## R/alfak-compat.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `alfak_neighbours` | function | 18-30 |  |  | sparse_joint_alfak_xval, as_alfak |
| `alfak_r2r` | function | 33-37 |  |  | sparse_joint_alfak_xval |
| `sparse_joint_relative_draws` | function | 48-82 |  | RNG | as_alfak |
| `alfak_symmetric_columns` | function | 89-99 |  |  | as_alfak |
| `alfak_efflux_constant` | function | 112-128 |  |  | as_alfak |
| `alfak_row_summary` | function | 130-135 |  |  | as_alfak |
| `sparse_joint_alfak_xval` | function | 150-239 |  |  | as_alfak |
| `as_alfak` | function | 337-465 | 是 |  | write_alfak_outputs, alfak_compat, predict.karst_alfak_landscape, fit_sparse_joint |
| `write_alfak_outputs` | function | 480-493 | 是 |  | alfak_compat |
| `alfak_compat` | function | 534-581 | 是 |  |  |
| `predict.karst_alfak_landscape` | function | 591-612 | 是 |  |  |
| `print.karst_alfak_output` | function | 615-621 | 是 |  |  |

## R/count-check.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `check_sparse_joint_counts` | function | 34-64 | 是 |  |  |
| `count_check_nearest_distance` | function | 70-83 |  |  | check_sparse_joint_counts |

## R/data-classes.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `new_karst_data` | function | 3-13 |  |  | prepare_karst_data, sparse_joint_keep_samples |
| `print.karst_data` | function | 16-23 | 是 |  |  |

## R/data-validate.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `validate_count_matrix` | function | 1-21 |  |  | prepare_karst_data |
| `validate_scalar` | function | 23-28 |  |  | as_alfak, check_sparse_joint_counts, resolve_sample_times, fit_sparse_joint_cohort, fit_sparse_joint_rates, fit_sparse_joint |
| `validate_beta` | function | 34-39 |  |  | prepare_karst_data, fit_sparse_joint_cohort, fit_sparse_joint |

## R/fit-karst.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `fit_karst` | function | 24-58 | 是 |  |  |

## R/forward-closed-form.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `sparse_division_kernel` | function | 26-59 |  | C++ | fit_sparse_joint_cohort, sparse_joint_series_method_base, sparse_joint_problem |
| `log_phi1` | function | 63-67 |  |  | closed_form_forward_r |
| `dyson_g2` | function | 72-78 |  |  | log_phi2 |
| `log_phi2` | function | 81-86 |  |  | closed_form_forward_r |
| `row_log_sum_exp` | function | 88-96 |  |  | closed_form_forward_r |
| `closed_form_forward_r` | function | 100-124 |  |  |  |
| `exact_forward_reference` | function | 129-153 |  | Matrix |  |

## R/karyotype-utils.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `parse_karyotypes` | function | 9-11 |  | C++ | alfak_neighbours, sparse_joint_relative_draws, as_alfak, predict.karst_alfak_landscape, check_sparse_joint_counts, prepare_karst_data, coerce_grf_karyotypes, fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, sparse_joint_cohort_histories, sparse_joint_pair_domains, sparse_joint_major, sparse_joint_domain, sparse_joint_domain_coverage, fit_sparse_joint, predict.karst_sparse_joint_fit |
| `format_karyotypes` | function | 18-20 |  | C++ | alfak_neighbours, predict.karst_alfak_landscape, sparse_division_kernel, format_grf_labels, predict.karst_sparse_joint_cohort_fit, one_step_neighbours, sparse_joint_major, sparse_joint_domain, predict.karst_sparse_joint_fit |
| `canonical_karyotype_labels` | function | 26-34 |  |  | prepare_karst_data, predict.karst_sparse_joint_cohort_fit, predict.karst_sparse_joint_fit |

## R/prepare-data.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `prepare_karst_data` | function | 21-40 | 是 |  | fit_sparse_joint_cohort, sparse_joint_cohort_histories, fit_sparse_joint |
| `resolve_sample_times` | function | 45-63 |  |  | prepare_karst_data |
| `karst_data_times` | function | 71-74 | 是 |  | as_alfak, sparse_joint_series_method_base, sparse_joint_keep_samples, sparse_joint_problem, fit_sparse_joint |

## R/simulate-l1-gp-landscape.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `validate_grf_integer` | function | 1-10 |  |  | simulate_grf_landscape, simulate_sparse_counts |
| `validate_grf_scalar` | function | 12-25 |  |  | simulate_grf_landscape, simulate_sparse_counts |
| `validate_grf_vector` | function | 27-33 |  |  | simulate_grf_landscape |
| `with_grf_seed` | function | 35-47 |  | RNG | simulate_grf_landscape, simulate_sparse_counts |
| `format_grf_labels` | function | 49-54 |  |  | simulate_grf_landscape, predict_landscape_fitness, grf_one_step_neighbours, grf_initial_karyotypes, simulate_sparse_counts |
| `coerce_grf_karyotypes` | function | 56-79 |  |  | predict_landscape_fitness, grf_one_step_neighbours, simulate_sparse_counts |
| `raw_grf_values` | function | 81-90 |  |  | simulate_grf_landscape, predict_landscape_fitness |
| `simulate_grf_landscape` | function | 123-206 | 是 | RNG | simulate_l1_gp_landscape |
| `grf_cache_signature` | function | 211-214 |  |  | simulate_grf_landscape, predict_landscape_fitness |
| `simulate_l1_gp_landscape` | function | 231-268 | 是 |  |  |
| `predict_landscape_fitness` | function | 281-318 | 是 |  | simulate_sparse_counts |

## R/simulate.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `grf_one_step_neighbours` | function | 1-27 |  |  | grf_initial_karyotypes |
| `grf_initial_karyotypes` | function | 29-43 |  |  | simulate_sparse_counts |
| `grf_multinomial_draw` | function | 49-62 |  | RNG | grf_downsample_population, grf_sample_counts, simulate_sparse_counts |
| `grf_population_env` | function | 64-70 |  |  | simulate_sparse_counts |
| `grf_env_population` | function | 72-78 |  |  | simulate_sparse_counts |
| `grf_add_count` | function | 80-89 |  |  | simulate_sparse_counts |
| `grf_viable` | function | 91-93 |  |  | simulate_sparse_counts |
| `grf_missegregate_daughters` | function | 95-109 |  | RNG | simulate_sparse_counts |
| `grf_downsample_population` | function | 111-117 |  |  | simulate_sparse_counts |
| `grf_sample_counts` | function | 119-124 |  |  | simulate_sparse_counts |
| `simulate_sparse_counts` | function | 175-332 | 是 | RNG |  |

## R/sparse-joint-cohort.R（不迁移）

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `fit_sparse_joint_cohort` | function | 71-422 | 是 | TMB |  |
| `predict.karst_sparse_joint_cohort_fit` | function | 446-549 | 是 |  |  |
| `print.karst_sparse_joint_cohort_fit` | function | 552-565 | 是 |  |  |
| `sparse_joint_cohort_coverage` | function | 569-580 |  |  | fit_sparse_joint_cohort, sparse_joint_cohort_histories |
| `sparse_joint_cohort_check_pairs` | function | 586-598 |  |  | fit_sparse_joint_cohort |
| `sparse_joint_cohort_histories` | function | 600-640 |  |  | fit_sparse_joint_cohort |

## R/sparse-joint-convergence.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `sparse_joint_convergence_tolerance` | value | 49-49 |  |  |  |
| `sparse_joint_tmb_state_names` | value | 52-53 |  | TMB |  |
| `sparse_joint_tmb_state` | function | 55-57 |  |  | sparse_joint_check_scales, sparse_joint_optimise_scales |
| `sparse_joint_restore_tmb_state` | function | 59-62 |  |  | sparse_joint_check_scales, sparse_joint_optimise_scales |
| `sparse_joint_fix_inner_start` | function | 66-70 |  | TMB | sparse_joint_check_scales, sparse_joint_optimise_scales, sparse_joint_inner_mode |
| `sparse_joint_inner_min_eigen` | function | 78-80 |  | TMB | sparse_joint_guarded_objective |
| `sparse_joint_inner_flattest` | function | 87-120 |  | TMB Matrix | sparse_joint_inner_min_eigen, sparse_joint_check_scales |
| `sparse_joint_latent_labels` | function | 128-145 |  |  | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_flat_direction` | function | 150-158 |  |  | sparse_joint_check_scales |
| `sparse_joint_flat_composition` | function | 163-175 |  |  | sparse_joint_check_scales |
| `sparse_joint_describe_flat` | function | 181-208 |  |  | sparse_joint_convergence_warning |
| `sparse_joint_projected_gradient` | function | 211-217 |  |  | sparse_joint_check_scales |
| `sparse_joint_fd_gradient` | function | 220-230 |  |  | sparse_joint_check_scales |
| `sparse_joint_check_scales` | function | 239-274 |  | TMB | sparse_joint_optimise_scales |
| `sparse_joint_guarded_objective` | function | 278-290 |  | TMB | sparse_joint_optimise_scales |
| `sparse_joint_compass_search` | function | 295-322 |  |  | sparse_joint_optimise_scales |
| `sparse_joint_optimise_scales` | function | 330-432 |  | TMB nlminb option | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_choose_candidate` | function | 441-449 |  | option | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_convergence_warning` | function | 453-466 |  |  | fit_sparse_joint_cohort, fit_sparse_joint |

## R/sparse-joint-rates.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `fit_sparse_joint_rates` | function | 67-92 | 是 |  | as_alfak, fit_sparse_joint |
| `predict.karst_sparse_joint_rate_fit` | function | 112-143 | 是 |  |  |
| `print.karst_sparse_joint_rate_fit` | function | 146-152 | 是 |  |  |

## R/sparse-joint-series-methods.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `sparse_joint_series_methods` | value | 8-8 |  |  |  |
| `sparse_joint_pair_domains` | function | 14-35 |  |  | sparse_joint_series_method_base |
| `sparse_joint_pair_coverage` | function | 39-47 |  |  | sparse_joint_check_pair_coverage, sparse_joint_pair_left_out |
| `sparse_joint_series_method_base` | function | 51-71 |  |  | sparse_joint_pairwise_problem, sparse_joint_state_space_problem |
| `sparse_joint_pairwise_problem` | function | 76-116 |  |  | sparse_joint_problem |
| `sparse_joint_state_space_problem` | function | 122-169 |  |  | sparse_joint_problem |
| `sparse_joint_series_method_expected` | function | 174-201 |  |  | sparse_joint_posterior |
| `sparse_joint_fill_latent` | function | 205-217 |  |  | fit_sparse_joint |
| `sparse_joint_first_sample_latent` | function | 222-231 |  |  | as_alfak |
| `sparse_joint_check_pair_coverage` | function | 236-245 |  |  | fit_sparse_joint |
| `sparse_joint_pair_left_out` | function | 260-282 |  |  | fit_sparse_joint |
| `sparse_joint_latent_states` | function | 286-299 |  |  | fit_sparse_joint |
| `sparse_joint_first_sample_zeta` | function | 303-307 |  |  | as_alfak |

## R/sparse-joint.R

| 名称 | 类型 | 行 | 导出 | 标记 | 调用者 |
|---|---|---|---|---|---|
| `one_step_neighbours` | function | 8-23 |  |  | sparse_joint_major, sparse_joint_domain |
| `sparse_joint_major` | function | 35-50 |  |  | sparse_joint_domain |
| `sparse_joint_major_distance` | function | 53-59 |  |  | fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, sparse_joint_domain_coverage, sparse_joint_design, fit_sparse_joint, predict.karst_sparse_joint_fit |
| `sparse_joint_domain` | function | 64-79 |  |  | fit_sparse_joint_cohort, sparse_joint_cohort_histories, sparse_joint_pair_domains, sparse_joint_pair_coverage, sparse_joint_check_coverage, sparse_joint_problem |
| `sparse_joint_domain_coverage` | function | 87-108 |  |  | sparse_joint_cohort_coverage, sparse_joint_pair_coverage, sparse_joint_check_coverage |
| `sparse_joint_coverage_error` | function | 115-140 |  |  | sparse_joint_cohort_coverage, sparse_joint_check_pair_coverage, sparse_joint_check_coverage |
| `sparse_joint_check_coverage` | function | 144-153 |  |  | fit_sparse_joint |
| `sparse_joint_format_coverage` | function | 156-161 |  |  | print.karst_sparse_joint_fit |
| `sparse_joint_design` | function | 168-173 |  |  | sparse_joint_relative_draws, sparse_joint_alfak_xval, fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, sparse_joint_series_method_base, sparse_joint_problem, predict.karst_sparse_joint_fit |
| `sparse_joint_initial_prior` | function | 179-181 |  |  | fit_sparse_joint_cohort, sparse_joint_pairwise_problem |
| `sparse_joint_initial_prior_matrix` | function | 183-196 |  |  | sparse_joint_initial_prior, sparse_joint_problem |
| `sparse_joint_inner_mode` | function | 214-234 |  | TMB | fit_sparse_joint_cohort, sparse_joint_posterior, sparse_joint_rebuild, fit_sparse_joint |
| `sparse_joint_latent_vector` | function | 242-247 |  |  | sparse_joint_rebuild |
| `sparse_joint_posterior` | function | 270-325 |  | TMB Matrix RNG | fit_sparse_joint_cohort, sparse_joint_rebuild, fit_sparse_joint |
| `sparse_joint_coef_covariance` | function | 329-333 |  | Matrix | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_default_calibration` | value | 362-362 |  |  |  |
| `sparse_joint_validate_calibration` | function | 366-378 |  |  | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_calibrate` | function | 382-387 |  |  | fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, fit_sparse_joint, predict.karst_sparse_joint_fit |
| `sparse_joint_relative` | function | 395-410 |  |  | fit_sparse_joint_cohort, fit_sparse_joint |
| `sparse_joint_relative_outside` | function | 414-418 |  |  | predict.karst_sparse_joint_fit |
| `sparse_joint_observed_level` | function | 426-432 |  |  | fit_sparse_joint |
| `sparse_joint_summary_total` | function | 436-440 |  |  | sparse_joint_observed_level |
| `sparse_joint_fitness_columns` | function | 444-452 |  |  | fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, fit_sparse_joint, predict.karst_sparse_joint_fit |
| `sparse_joint_self_probability` | function | 455-458 |  | C++ | as_alfak, predict.karst_alfak_landscape, fit_sparse_joint_cohort, predict.karst_sparse_joint_cohort_fit, predict.karst_sparse_joint_fit |
| `sparse_joint_history_candidates` | function | 466-469 |  |  | fit_sparse_joint |
| `sparse_joint_template` | function | 475-478 |  |  | sparse_joint_problem |
| `sparse_joint_major_columns` | function | 484-487 |  |  | sparse_joint_problem, fit_sparse_joint |
| `sparse_joint_uncounted_samples` | function | 496-502 |  |  | fit_sparse_joint |
| `sparse_joint_keep_samples` | function | 507-514 |  |  | fit_sparse_joint |
| `sparse_joint_intermediate_dispersion` | function | 520-531 |  |  | fit_sparse_joint |
| `sparse_joint_used_counts` | function | 541-555 |  |  | sparse_joint_problem, fit_sparse_joint |
| `sparse_joint_problem` | function | 560-626 |  |  | sparse_joint_alfak_xval, sparse_joint_first_sample_latent, sparse_joint_first_sample_zeta, sparse_joint_rebuild, fit_sparse_joint |
| `sparse_joint_rebuild` | function | 636-698 |  | TMB | sparse_joint_alfak_xval, as_alfak, fit_sparse_joint |
| `fit_sparse_joint` | function | 989-1401 | 是 | TMB nlminb option | as_alfak, alfak_compat, fit_karst, fit_sparse_joint_cohort |
| `predict.karst_sparse_joint_fit` | function | 1427-1474 | 是 |  |  |
| `print.karst_sparse_joint_fit` | function | 1477-1491 | 是 |  |  |

## 测试

`移植` 一列按 HANDOFF.md 4.2 节去掉 cohort 部分。

### tests/testthat/helper-sparse-joint.R

| 测试 | 行 | 移植 |
|---|---|---|
| (helper) labelled_counts | 3-8 | yes |
| (helper) simulate_sparse_joint_counts | 14-30 | yes |
| (helper) simulate_balanced_counts | 35-50 | yes |
| (helper) simulate_cohort_counts | 56-76 | no (cohort) |
| (helper) simulate_series_counts | 81-99 | yes |
| (helper) renamed_like | 104-107 | yes |

### tests/testthat/test-alfak-compat.R

| 测试 | 行 | 移植 |
|---|---|---|
| neighbours follow alfak's order and rules | 1-6 | yes |
| a rebuilt posterior reproduces the fit | 8-22 | yes |
| as_alfak returns alfak's objects for the model's posterior | 24-75 | yes |
| as_alfak keeps the model's level or reports birth rates on request | 77-97 | yes |
| as_alfak reports frequencies at time 0 as alfak does | 99-128 | yes |
| alfak_compat() takes alfak's input and writes its files | 130-166 | yes |
| as_alfak() and alfak_compat() take fits of more than two time points | 168-205 | yes |

### tests/testthat/test-api-cleanup.R

| 测试 | 行 | 移植 |
|---|---|---|
| legacy exports are absent | 1-21 | yes |

### tests/testthat/test-forward-closed-form.R

| 测试 | 行 | 移植 |
|---|---|---|
| (helper) small_domain | 1-4 | yes |
| (helper) dense_reference | 6-15 | yes |
| sparse kernel reproduces the exact single-daughter kernel on one-step edges | 17-39 | yes |
| phi functions equal the time-ordered integrals they stand for | 41-62 | yes |
| uniformisation reference agrees with a dense matrix exponential | 64-75 | yes |
| closed form converges to the exact solution as the flux vanishes | 77-96 | yes |
| a deleterious neighbour settles at mutation-selection balance with its parent | 98-112 | yes |

### tests/testthat/test-history-choice.R

| 测试 | 行 | 移植 |
|---|---|---|
| every fit records the objective at the posterior mode it reports | 1-13 | yes |
| candidate histories are compared on the objective at the posterior mode | 15-45 | yes |
| a candidate whose scale optimisation gave no finite objective is not compared | 47-68 | yes |
| the trends of a cohort are compared on the objective at the posterior mode | 70-107 | no (cohort) |

### tests/testthat/test-karyotype-utils.R

| 测试 | 行 | 移植 |
|---|---|---|
| karyotype parsing round-trips | 1-5 | yes |

### tests/testthat/test-rebuild-stored-mode.R

| 测试 | 行 | 移植 |
|---|---|---|
| a fit's own posterior is rebuilt at its stored mode, with no search | 6-45 | yes |
| the posterior of a fit at a stored mode that is not exact is the fit's own | 47-58 | yes |

### tests/testthat/test-review-fixes-2.R

| 测试 | 行 | 移植 |
|---|---|---|
| fit_karst() fits prepared data at the data's own missegregation rate | 3-15 | yes |
| karyotype labels are read one way only | 17-43 | yes |
| the missegregation rate of a fit must be positive | 45-51 | yes |
| a free balance level is refused for the pairwise and state-space methods | 53-60 | yes |
| (helper) mock_inner_objective | 64-76 | yes |
| the posterior mode is searched once more from the initial values before it fails | 78-95 | yes |
| a candidate history whose posterior mode cannot be found is not chosen | 97-129 | yes |
| the second search uses TMB's real objective and the problem's own initial values | 131-190 | yes |
| the initial latent values are ordered as in TMB's parameter vector | 192-215 | yes |
| candidates that failed beside ones whose optimisation failed, and the restart record | 217-240 | yes |
| a posterior mode found by the second search is recorded | 242-262 | yes |
| the first fits of the intermediate weight and of initial_state = 'auto' fall back on failure | 264-289 | yes |
| a fit with no usable candidate history stops | 291-297 | yes |
| time points left out keep their names, their karyotypes and alfak's outputs | 299-332 | yes |
| karyotypes seen only between the first and the last time point count as unobserved | 334-351 | yes |
| the efflux correction is alfak's: the x0-weighted mean at time 0 is g0 | 353-388 | yes |
| cross-validation holds out each fold's prior and drift weight as the fit treats unseen states | 390-429 | yes |
| the simulators keep their parameters and their samples | 431-454 | yes |

### tests/testthat/test-review-fixes.R

| 测试 | 行 | 移植 |
|---|---|---|
| an explicit missegregation rate replaces the one stored in prepared data | 3-18 | yes |

### tests/testthat/test-simulate-l1-gp-landscape.R

| 测试 | 行 | 移植 |
|---|---|---|
| paper-style GRF landscape is lightweight and reproducible | 1-11 | yes |
| GRF oracle anchors founder and is deterministic | 13-33 | yes |
| legacy landscape constructor now returns lazy GRF | 35-49 | yes |
| lazy GRF sparse simulator returns valid observed truth | 51-77 | yes |

### tests/testthat/test-simulate.R

| 测试 | 行 | 移植 |
|---|---|---|
| GRF oracle is reproducible by seed | 1-6 | yes |
| sparse simulator returns valid two-timepoint counts | 8-25 | yes |

### tests/testthat/test-sparse-joint-cohort.R

| 测试 | 行 | 移植 |
|---|---|---|
| a cohort of one dataset with full sharing is the single-dataset fit | 1-15 | no (cohort) |
| datasets share the karyotypes they have in common | 17-48 | no (cohort) |
| prediction takes the best-informed source | 50-72 | no (cohort) |
| invalid cohort inputs are rejected | 74-84 | no (cohort) |
| the trend is per dataset, shared, or chosen by marginal likelihood | 86-119 | no (cohort) |
| with balance, each dataset's history is chosen as for a single dataset | 121-139 | no (cohort) |
| a cohort dataset with a timepoint outside its domain stops with an error naming it | 141-155 | no (cohort) |
| a cohort fit reports the domain coverage of each dataset | 157-179 | no (cohort) |

### tests/testthat/test-sparse-joint-convergence.R

| 测试 | 行 | 移植 |
|---|---|---|
| (helper) fake_hessian_env | 1-6 | yes |
| inverse iteration finds the smallest eigenvalue of a sparse Hessian | 8-37 | yes |
| inverse iteration also gives the flat direction | 39-62 | yes |
| the flat direction is named by karyotype, cells and share | 64-103 | yes, without lines 92-102 (cohort) |
| finite-difference and projected gradients respect the bounds | 105-114 | yes |
| compass search finds a bounded minimum and reports improvements | 116-129 | yes |
| the candidate kept is the lowest finite objective, or the best accepted one on request | 131-140 | yes |
| checking the scales leaves TMB where the fit left it | 142-166 | yes |
| a fit reports how its scale optimisation was checked | 168-184 | yes |
| a fit that is not accepted warns once, with its reasons | 186-204 | yes |
| a fit with a near-singular inner Hessian names the flat direction | 206-241 | yes |
| when nlminb fails, the prior-centred scales are kept, as before the check | 243-264 | yes |
| searching again is available on request | 266-277 | yes |

### tests/testthat/test-sparse-joint-rates.R

| 测试 | 行 | 移植 |
|---|---|---|
| a range of birth rates is fitted at log-spaced rates, with point predictions at the lowest | 1-37 | yes |
| the assumed rate scales the deficit of neighbours and leaves their ranking | 39-55 | yes |
| invalid ranges are rejected | 57-63 | yes |
| a range of missegregation rates is fitted like a range of birth rates | 65-93 | yes |
| exactly one range is given, and not together with its rate | 95-102 | yes |

### tests/testthat/test-sparse-joint-series-methods.R

| 测试 | 行 | 移植 |
|---|---|---|
| two time points give the two-time-point fit with every series method | 1-12 | yes |
| pair domains give each pair of time points its own majors | 14-25 | yes |
| the pairwise method fits each pair as two time points with one landscape | 27-70 | yes, without lines 35-49 (cohort) |
| the state-space method follows the anchored fit when there is no drift | 72-99 | yes |
| predictions and alfak outputs work for every series method | 101-135 | yes |
| a first time point outside every pair's domain stops the fit | 137-145 | yes |
| (helper) replace_series_sample | 148-155 | yes |
| a time point in between with no cell in its pair's domain is left out | 157-181 | yes |
| a last time point with no cell in its pair's domain pairs with a time point further in | 183-195 | yes |

### tests/testthat/test-sparse-joint-series.R

| 测试 | 行 | 移植 |
|---|---|---|
| (helper) series_settings | 1-6 | yes |
| count data take two or more timepoints | 8-29 | yes |
| the series template gives the two-timepoint objective for two timepoints | 31-55 | yes |
| two-timepoint fits are the same with either template | 57-75 | yes |
| compiled series model evaluates the closed form at every time point | 77-102 | yes |
| a series fit uses every time point | 104-143 | yes |
| a karyotype first seen later is centred on its frequency then | 145-152 | yes |
| a timepoint with no cell next to a major names the timepoints the majors come from | 154-163 | yes |
| functions that take two timepoints say so | 165-172 | yes, without lines 167-168 (cohort) |
| intermediate time points can count for the majors only | 174-218 | yes |
| a time point in between with no counted cell is left out | 220-273 | yes |
| with majors only, a time point in between with cells on their neighbours only is left out | 275-294 | yes |

### tests/testthat/test-sparse-joint.R

| 测试 | 行 | 移植 |
|---|---|---|
| the domain is the major karyotypes and all their one-step neighbours | 1-22 | yes |
| majors need both a share of a time point and a number of cells | 24-44 | yes |
| compiled forward model equals the R closed form | 46-77 | yes |
| fit recovers the order of clones and the deficit of their neighbours | 79-103 | yes |
| prediction returns the posterior inside the domain and the trend outside it | 105-120 | yes |
| an observed karyotype outside the domain is placed at the level of observed neighbours | 122-136 | yes |
| invalid inputs are rejected | 138-144 | yes |
| the balance prior is the closed-form influx from the majors | 146-182 | yes |
| with balance, the first time point tells neighbours apart | 184-205 | yes |
| reported uncertainty of relative fitness adds model error to the model's own | 207-252 | yes |
| uncertainty is finite and no wider than the prior where the counts fit poorly | 254-265 | yes |
| fit_karst fits the sparse joint model | 267-285 | yes |
| mutation-selection balance is the default initial state | 287-292 | yes |
| by default the history is chosen among the interval and longer division counts | 294-323 | yes |
| a timepoint with no cell next to a major karyotype stops with a typed error | 325-352 | yes |
| a fit reports how many cells of each timepoint lie inside the domain | 354-372 | yes |
| counts are checked per time point before fitting | 374-390 | yes |
| a free level of the balance prior and initial_state = "auto" | 392-438 | yes |

