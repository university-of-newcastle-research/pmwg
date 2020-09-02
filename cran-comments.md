## Resubmission
This is a resubmission. In this version I have:

* Fixed references in the ‘Description’ of the DESCRIPTION file, now in the form
  authors (year) <doi:...> - See R CMD check results below for introduced NOTE

* Examples for unexported functions have now been omitted.

* I have rerun the r-devel build on winbuilder and local checks

## Test environments
* Ubuntu 18.04.5 LTS (local system), R version 3.6.3 (2020-02-29)
* Ubuntu 16.04.6 LTS (on travis-ci), R version 3.6.3 (2017-01-27)
* Ubuntu 16.04.6 LTS (on travis-ci), R version 4.0.0 (2020-04-24)
* Ubuntu 16.04.6 LTS (on travis-ci), R Under development (unstable) (2020-08-04 r78970)
* macOS High Sierra 10.13.6 (on travis-ci), R version 3.6.3 (2020-02-29)
* macOS High Sierra 10.13.6 (on travis-ci), R version 4.0.2 (2020-06-22)
* Windows 10 x64 (build 19041) (local system), R version 4.0.2 (2020-06-22)
* win-builder, R Under development (unstable) (2020-08-23 r79071)
* R-hub (rhub::check_for_cran()): Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* R-hub (rhub::check_for_cran()): Ubuntu Linux 16.04 LTS, R-release, GCC
* R-hub (rhub::check_for_cran()): Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release. (resubmission)
* Possibly mis-spelled words in DESCRIPTION: (Author names for paper reference)
  * Gunawan (15:8)
  * Kohn (15:50)
