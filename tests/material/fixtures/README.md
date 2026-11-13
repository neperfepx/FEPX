# Crystal fixtures

`schmid_reference.txt`: FEPX's Schmid tensors d ⊗ n of every slip system
(`src/material/crystal_mod.f90::crystal_init`, crystal frame), for fcc, bcc with the
{112}<111> systems (24), hcp with c/a = 1.587 and bct with c/a = 0.5456. Each case is a header line
`<structure> <num_slip> <c_over_a>` followed by one line per slip system (the 9 entries row by row).
`tests/material/test_crystal.f90` checks that the slip systems, and their order, do not change
(1e-15).

Generated on 2026-10-25 (ROADMAP step 9.13, D91) with `dump_schmid.f90` (this directory), built
against the objects of a Release build of `ac211915` (tag `v2-final-comparison`):

```
mpif90 -O2 -I<build dir>/mod dump_schmid.f90 <build dir>/src/libfepx_core.a -o dump_schmid
./dump_schmid > schmid_reference.txt
```

It replaces `v2_schmid.txt`, the tensors built by v2-fixed (`crys_type_mod2.f90::schmidtensors`),
against which FEPX's ported slip systems were checked from step 1.4 to step 9.13: the two files
are identical. The v2 fixture and its program are in the tag `v2-final-comparison`.

`kd25_vertices.txt`: single crystal yield-surface vertices computed by Kasemer & Dawson's own
implementation of their method (KD25, arXiv:2502.19531; github.com/mkasemer/scys, GPL-3, commit
`c70d070`, 2025-01-21; the code is not part of FEPX, user-approved development cross-check,
ROADMAP D57). Five cases: fcc with unit strengths; bcc {110}<111> with system 1 5 % weaker;
hcp c/a 1.587 with basal 1, prismatic 1.25, pyramidal <c+a> 1.7; hcp c/a 1.59 with pyramidal
<c+a> only (basal and prismatic at 100, out of reach); bct with unit strengths (KD25's bct
ignores c/a: c/a = 1). Each case is a header line `<structure> <c_over_a> <num_slip>
<num_vertices>`, a line of the strengths of the systems, and one line per vertex (the deviatoric
5-vector of KD25 Eq. 8, both signs of each pair). KD25's code uses absolute tolerances of 1e-4
and rounds candidates to 4 decimals, so `tests/material/test_yield_vertices.f90` compares the
vertex sets to 1e-3. The systems of KD25's code are FEPX's (and v2's), in the same order, for
these cases, except KD25's bct system 26: its normal (1 2 1) and direction [1 0 1] are not
orthogonal (an error in KD25's table, found in step 6.7; with it KD25's code gives 664 bct
vertices instead of 680). `kd25_dump_vertices.py` replaces it by FEPX's and v2's (-1 2 1)[1 0 1]
before the computation, leaving KD25's code unchanged.

Generated on 2026-10-11 with `kd25_dump_vertices.py` (this directory), Python 3.14, NumPy 2.3.5:

```
git clone https://github.com/mkasemer/scys <clone>      # commit c70d070
python3 -I kd25_dump_vertices.py <clone> > kd25_vertices.txt    (about 2 minutes)
```
