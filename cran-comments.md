## R CMD check results

 checking CRAN incoming feasibility ... [7s/32s] NOTE
 Maintainer: ‘Gavin Cooper <gavin@gavincooper.net>’
 
 Found the following (possibly) invalid URLs:
   URL: https://www.pnas.org/content/105/45/17538
     From: man/forstmann.Rd
           man/sampled_forstmann.Rd
     Status: 403
     Message: Forbidden

0 errors ✔ | 0 warnings ✔ | 1 note ✖


The pnas reference site appears to be blocking automated webcrawler requests, as
various canonical links attempted lead to the same error and all links tried
manually through a browser successfully open the article website.

## Test environments
* Ubuntu Ubuntu 20.04.6 LTS (local system), R version 4.3.1 (2023-06-16)
* R-hub (rhub::check_for_cran()): Windows Server 2022, R-devel, 64 bit
* R-hub (rhub::check_for_cran()): Ubuntu Linux 20.04.1 LTS, R-release, GCC
* R-hub (rhub::check_for_cran()): Fedora Linux, R-devel, clang, gfortran

## dependency results

There are 0 reverse dependencies.

Some dependencies have been removed (summary see below)

MCMCpack
