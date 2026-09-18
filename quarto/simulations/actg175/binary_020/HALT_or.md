# HALT — ACTG175 binary/OR Stage 2

- Halted at (UTC): 2026-09-18T00:23:34Z
- Cell: orfs_or150_n500
- Assertion: render batch_1001_2000 failed (rc=1; 124 = timeout 9915s)
- HEAD at halt: 50c25059
- Cumulative render wall: 14447 s (ceiling 105165 s)
- Completed cells stay committed; nothing is rolled back.

## Observed values

```
[31m

processing file: sim_fs_mr_field_or_template.qmd
[39m1/57                                 
2/57 [setup-knobs]                   
3/57                                 
4/57 [build-dgm]                     
5/57                                 
6/57 [machinery]                     
7/57                                 
8/57 [run-batch]                     
[31mError:
! OR positivity violated -- non-positive finite values in: or_H_lo (1)
[39m[31m
Quitting from sim_fs_mr_field_or_template.qmd:919-1171 [run-batch]
Execution halted
[39m[33mWARN: Error encountered when rendering files[39m
WALL_SECONDS=1054 RC=1 PEAK_MB=80219 OUT=fs_effMaxSG_mr_field_or150_n500_nb20_orfs_batch_1001_2000.html
```
