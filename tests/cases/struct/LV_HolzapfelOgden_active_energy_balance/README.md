# Problem Description

This case is a variant of
[LV_HolzapfelOgden_active](../LV_HolzapfelOgden_active/README.md) with
`Energy_balance_face` set to `top` on the endocardial follower-pressure boundary.
The option applies a basal traction to balance the endocardial pressure resultant,
using pressure and geometry from the previous accepted state. The basal load is
held fixed during the timestep's Newton iterations.

The mesh, fibers, and prescribed pressure and active-stress histories are shared
with the original case through relative paths. All other input settings are
unchanged. The regression advances one timestep of 0.001 s and compares the output
with the enabled reference solution in `result_001.vtu`.
