natcpp v0.3.1 is a patch release to fix an installation failure seen on
CRAN's r-devel-linux-x86_64-debian-clang check for v0.3.0.

## Resubmission

This release fixes a package load failure caused by an unresolved
`__atomic_compare_exchange` symbol with Debian clang 22 on CRAN's build
system. The package now uses a configure-time probe and links `libatomic`
only on platforms where the C++ atomic compare/exchange test requires it.

## Test environments
* win-builder (devel), <https://win-builder.r-project.org/FpXFoajGDVzl> (Status: OK)
* local R installation, R 4.6.1
* Continuous Integration via GitHub actions
  * windows-latest (release)
  * macOS-latest (release)
  * ubuntu-latest (release)
  * ubuntu-latest (devel)
  * ubuntu-latest (oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes
