# Manuscript-scale settings recorded from the production workflows.
# The portable reference implementation is not the authoritative production
# code; see ../../production/ and ../../data/production_results/raw/.
SIGA_MANUSCRIPT <- list(
  siga_s_continuous = list(
    target = "S", n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260725L,
    calibration_type = "one_path"
  ),
  siga_s_binary_large_n = list(
    target = "S", n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260726L,
    calibration_type = "one_path"
  ),
  siga_s_binary_small_n = list(
    target = "S", n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260728L,
    calibration_type = "one_path"
  ),
  full_grid_r = list(
    target = "R", n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260801L,
    calibration_type = "three_path", epsilon = "1/n"
  ),
  pair_path_stress = list(
    n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260729L,
    epsilon = "1/n"
  ),
  swift_direct_inspired = list(
    n_outer = 100000L, B_randomization = 4999L,
    B_calibration = 100000L, base_seed = 20260724L,
    epsilon = "1/n"
  )
)
