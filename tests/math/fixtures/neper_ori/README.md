# Neper orientation fixtures

The same 20 random orientations written by Neper in every orientation descriptor and both
conventions, one orientation per line (Neper `.ori` plain format, 12 decimals). Used by
`tests/math/test_orientation.f90` (ROADMAP step 1.3): every file must give the same `rot_c2s`,
and FEPX must write the same values back.

Generated on 2026-10-08 with Neper 5.0.1-520-devel:

```
for des in rodrigues euler-bunge euler-kocks euler-roe rotmat axis-angle quaternion; do
  for conv in passive active; do
    neper -T -n 20 -id 7 -ori random -oridescriptor $des:$conv -format ori -o r_${des}_${conv}
  done
done
```

Checked when generated: `r_rotmat_active` is exactly the transpose of `r_rotmat_passive`;
`r_rotmat_passive` rows are orthonormal to 1.2e-12; Bunge angles of `r_euler-bunge_passive` give
the matrices of `r_rotmat_passive` (v2's `euler_bunge_to_rotmat` formula) to 5e-13.
`rotmat` lines are the matrix g row by row; `quaternion` lines are q0 q1 q2 q3; `axis-angle` lines
are the axis and the angle in degrees; Euler angles are in degrees.
