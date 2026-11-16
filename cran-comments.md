## Submission

This is an update of ICtest from version 0.3-7 to 0.3-8. It adds the new functions
FOBIaug(), SIRaug(), PCAaugHD() and PCAlrt(), fixes several bugs and adds regression tests.

Note that PCAschott() now computes the high-dimensional test of Schott (2006), as was
always documented. The likelihood ratio test it computed before is available as the
new function PCAlrt(). Results of PCAschott() therefore differ from earlier versions,
which is stated in the ChangeLog.

## Test environments

* local Windows 11, R-devel (2026-07-24 r90297 ucrt)
* GitHub Actions: Ubuntu (R-devel, R-release, R-oldrel-1), macOS (R-release),
  Windows (R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

All five reverse dependencies (KernelICA, singR, ssaBSS, tensorBSS, tsBSS) were
checked with R CMD check against ICtest 0.3-8, and none of them showed any new
problems.
