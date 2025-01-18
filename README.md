# R3BRoot Software [![license](https://alfa-ci.gsi.de/shields/badge/license-GPL--3.0-orange.svg)](COPYRIGHT) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5549469.svg)](https://doi.org/10.5281/zenodo.5549469) [![OpenSSF Best Practices](https://www.bestpractices.dev/projects/9851/badge)](https://www.bestpractices.dev/projects/9851) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)

[![CI-CD](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/main.yml/badge.svg)](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/main.yml) [![Static Analysis](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/static_analysis.yml/badge.svg)](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/static_analysis.yml) [![Validate Codemeta](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/codemeta_validate.yaml/badge.svg)](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/codemeta_validate.yaml) [![Cleanup Caches on PR Close](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/cleanup_cache.yml/badge.svg)](https://github.com/R3BRootGroup/R3BRoot/actions/workflows/cleanup_cache.yml)

## The R3BRoot Framework

The R3BRoot software is based on the [FairRoot](https://github.com/FairRootGroup/FairRoot) framework and can be used to perform Monte Carlo simulations and experimental data analysis of the R3B (Reactions with Relativistic Radioactive Beams) nuclear physics experiments at the GSI-FAIR research center (Facility for Antiproton and Ion Research). The user can create simulated data and/or perform analysis with the same framework. Geant3 and [Geant4](https://geant4.web.cern.ch/) transport engines are supported, however, the implemented tasks that create simulated data do not depend on a particular Monte Carlo engine. The framework delivers base classes which enable the users to construct their detectors and/or analysis tasks in a simple way, it also delivers some general functionality like track visualization. Moreover, an interface for reading experimental and/or simulated magnetic field maps is also implemented. R3BRoot software is a source distribution with recurring releases for macOS and Linux.

## Discussion Forum

For the software-related user support you can post a new topic on our [forum](https://forum.gsi.de/index.php?t=index&cat=40&).

## License

R3BRoot is distributed under the terms of the GNU General Public Licence version 3 ([GPLv3](https://github.com/R3BRootGroup/R3BRoot/blob/dev/LICENSE)).

## Code of Conduct

We are committed to fostering a welcoming and inclusive environment in the R3BRoot project. Please take a moment to review our [Code of Conduct](./CODE_OF_CONDUCT.md), which outlines our expectations for all contributors and community members.

## Release Information

Please visit [releases](https://github.com/R3BRootGroup/R3BRoot/releases)

## Download

~~~bash
git clone https://github.com/R3BRootGroup/R3BRoot.git
cd R3BRoot
git clone https://github.com/R3BRootGroup/macros.git
~~~

## Using the Project Template

R3BRoot delivers various project templates that can be used as a starting point for anybody who would like to build simulation and reconstruction algorithms. The project Templates are located in the [R3BRoot/template](template) directory.


## Step by Step Installation

### Required Software

First, you will need to install FairSoft and FairRoot. For more details:

1. Install [FairSoft](https://github.com/FairRootGroup/FairSoft), required version Nov22p1 or later

2. Install [FairRoot](http://fairroot.gsi.de), required version 18.2.1 or later

### Configure and Compile

~~~bash
export SIMPATH=%PATH_TO_FAIRSOFT%
export FAIRROOTPATH=%PATH_TO_FAIRROOT%
# from %R3BRoot_DIRECTORY%
mkdir build
cmake -S . -B ./build
cmake --build ./build -- -j${number_of_threads}
source build/config.sh
~~~

For the building with Conan package manager, see [this instruction](doc/conan_usage.md).

To run the detector tests do:

~~~bash
 # In the same shell from %R3BRoot_DIRECTORY%:
cd build
ctest --parallel ${number_of_threads}
~~~

## Simulations and Data Analysis

This is performed from the GitHub parameter and data analysis repository, which contains all the macros and parameter files needed by the user to carry out the simulations and data analysis of each experiment. There is one repository per experiment, please, visit the [R3B-Wiki](https://wiki.r3b-nustar.de) for more details.

## Contributing

Please ask your questions, request features, and report issues by [creating a github issue](https://github.com/R3BRootGroup/R3BRoot/issues/new/choose).

## Code Formatting

The R3BRoot project uses clang-format-15 to ensure a common code formatting. The script "clang-format-all.sh" can be used for this purpose: 
~~~bash
source util/clang-format-all.sh
~~~

## Security Policy

The R3BRoot project is focused on scientific computing and does not involve typical security concerns like handling sensitive data. However, we are committed to maintaining the integrity and quality of our codebase.

If you discover any potential issues with this project, please refer to our [Security Policy](./SECURITY.md) for instructions on how to report them.

Thank you for helping us keep R3BRoot secure and reliable!

## More Information

* [Static analyzer using Clang-tidy](config/clang_tidy/README.md)
* [CMake build system for R3BRoot](doc/cmake_usage.md)
* [How to use an unmerged pull request](doc/git_usage.md#fetch-the-update-from-an-unmerged-pull-request-pr)
