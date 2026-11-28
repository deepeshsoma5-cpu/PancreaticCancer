# Integrity report

Generated 2026-11-28T01:44Z by code/40_integrity_audit.sh from the git history.

## Key commits

~~~
521df7d  2026-11-23T03:29:29Z  Lock signature (fit on TCGA only, before external data is scored)
85ab2ce  2026-11-23T03:29:55Z  External validation on GSE62452
ac526e2  2026-11-24T22:34:59Z  Score locked signatures on GSE21501 (fresh) and GSE62452
de61008  2026-11-25T18:59:04Z  Pre-register fresh-cohort test of locked signature (before any new cohort is scored)
70e8749  2026-11-25T19:14:51Z  Plan additions (exclusions, replication rule) before scoring; inspect GSE57495 and GSE79668
383af9f  2026-11-25T19:18:32Z  Fix duplicate screen thresholds before running the screen
02ba724  2026-11-25T19:25:15Z  Record screen outcome before scoring
e9f8aea  2026-11-25T19:25:31Z  Score locked signature once on GSE57495 (plus recompute GSE62452, GSE21501) under the pre-committed rule
ea823fb  2026-11-26T02:46:13Z  Pre-specify final pathway analysis and close validation (before running)
59a6bfb  2026-11-26T02:46:31Z  Final pathway analysis (patient level, pre-specified) and auto-generated final results
5ec9d27  2026-11-27T19:59:15Z  Pre-specify final pathway analysis and close validation (before running)
27f2c2f  2026-11-27T20:04:26Z  Revert "Pre-specify final pathway analysis and close validation (before running)"
~~~

## Order checks (author timestamps)

~~~
PASS        signature locked before the GSE62452 scoring
PASS        signature locked before the GSE21501 scoring
PASS        signature locked before the GSE57495 scoring
PASS        fresh cohort plan before the GSE57495 scoring
PASS        replication rule and exclusions before the GSE57495 scoring
PASS        duplicate screen thresholds before the GSE57495 scoring
PASS        screen outcome before the GSE57495 scoring
PASS        final pathway plan before the final pathway results
~~~

## Locked signature file

~~~
PASS        results/10_locked_signature.csv is identical at HEAD and at commit 521df7d
blob id at HEAD: a381ab4b6726821541e22de012ef4dcc69b2373f
~~~

## Notes that belong in the paper and the logbook

- GSE62452 (85ab2ce) and GSE21501 (ac526e2) were scored before the fresh cohort plan (de61008) and before the replication rule was written (70e8749). Only GSE57495 was scored after both existed. The criterion is the standard one (the confidence interval excludes the null and the HR is above 1), but for the first two cohorts it was applied after their results were known.
- The decision to stop validation at three cohorts was made after seeing the scores.
- Commit 5ec9d27 accidentally overwrote docs/plan_final.md with a new timestamp and was reverted in 27f2c2f. The original plan commit is ea823fb. The final pathway plan was committed 18 seconds before its results (59a6bfb), so it documents the rules but is not an advance registration.
- Git author timestamps can be edited. The GitHub push history is the stronger record.
