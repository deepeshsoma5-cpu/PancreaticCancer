#!/bin/bash
# checks the order of the key commits and that the locked signature is unchanged
cd "$(dirname "$0")/.." || exit 1
ts() { git show -s --format=%at "$1" 2>/dev/null; }
show() { git show -s --format="%h  %aI  %s" "$1" 2>/dev/null || echo "$1  NOT FOUND"; }
before() {
  a=$(ts "$2"); b=$(ts "$3")
  if [ -n "$a" ] && [ -n "$b" ] && [ "$a" -lt "$b" ]; then
    echo "PASS        $1"
  else
    echo "NOT TRUE    $1 (or a commit was not found)"
  fi
}
{
echo "# Integrity report"
echo
echo "Generated $(date -u +%Y-%m-%dT%H:%MZ) by code/40_integrity_audit.sh from the git history."
echo
echo "## Key commits"
echo
echo '~~~'
for h in 521df7d 85ab2ce ac526e2 de61008 70e8749 383af9f 02ba724 e9f8aea ea823fb 59a6bfb 5ec9d27 27f2c2f; do show "$h"; done
echo '~~~'
echo
echo "## Order checks (author timestamps)"
echo
echo '~~~'
before "signature locked before the GSE62452 scoring" 521df7d 85ab2ce
before "signature locked before the GSE21501 scoring" 521df7d ac526e2
before "signature locked before the GSE57495 scoring" 521df7d e9f8aea
before "fresh cohort plan before the GSE57495 scoring" de61008 e9f8aea
before "replication rule and exclusions before the GSE57495 scoring" 70e8749 e9f8aea
before "duplicate screen thresholds before the GSE57495 scoring" 383af9f e9f8aea
before "screen outcome before the GSE57495 scoring" 02ba724 e9f8aea
before "final pathway plan before the final pathway results" ea823fb 59a6bfb
echo '~~~'
echo
echo "## Locked signature file"
echo
echo '~~~'
if git cat-file -e 521df7d:results/10_locked_signature.csv 2>/dev/null; then
  if git diff --quiet 521df7d HEAD -- results/10_locked_signature.csv; then
    echo "PASS        results/10_locked_signature.csv is identical at HEAD and at commit 521df7d"
  else
    echo "FAIL        results/10_locked_signature.csv differs from commit 521df7d"
  fi
  echo "blob id at HEAD: $(git hash-object results/10_locked_signature.csv)"
else
  echo "NOT TRUE    could not read the signature file at commit 521df7d"
fi
echo '~~~'
echo
echo "## Notes that belong in the paper and the logbook"
echo
echo "- GSE62452 (85ab2ce) and GSE21501 (ac526e2) were scored before the fresh cohort plan (de61008) and before the replication rule was written (70e8749). Only GSE57495 was scored after both existed. The criterion is the standard one (the confidence interval excludes the null and the HR is above 1), but for the first two cohorts it was applied after their results were known."
echo "- The decision to stop validation at three cohorts was made after seeing the scores."
echo "- Commit 5ec9d27 accidentally overwrote docs/plan_final.md with a new timestamp and was reverted in 27f2c2f. The original plan commit is ea823fb. The final pathway plan was committed 18 seconds before its results (59a6bfb), so it documents the rules but is not an advance registration."
echo "- Git author timestamps can be edited. The GitHub push history is the stronger record."
} > docs/INTEGRITY_REPORT.md
cat docs/INTEGRITY_REPORT.md
