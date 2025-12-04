sciLoc, a Seiscomp plugin for iLoc
==================================

[![License](https://img.shields.io/badge/License-BSD%203--Clause-3da639?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyMy4yIj48ZGVmcy8+PHBhdGggZmlsbD0iI2VlZSIgZD0iTTAgMTIuMUMuMSA1LjUgNC43LjggMTAuNC4xYzYuNy0uOSAxMi40IDMuNyAxMy41IDkuNyAxIDUuOC0yLjEgMTEuMS03LjMgMTMuMy0uNC4yLS43LjEtLjktLjRMMTMgMTZjLS4xLS40IDAtLjYuMy0uOCAxLjItLjUgMS45LTEuNCAyLjEtMi43LjMtMS45LTEtMy43LTIuOS00aC0uMmMtMS44LS4yLTMuNSAxLjEtMy44IDIuOS0uMyAxLjYuNSAzLjEgMiAzLjguNS4yLjYuNC40LjlsLTIuNiA2LjhjLS4xLjMtLjQuNS0uOC4zLTIuNy0xLjEtNS0zLjEtNi4zLTUuOEMuMSAxNSAuMSAxMy4xIDAgMTIuMXoiLz48L3N2Zz4=)](https://opensource.org/licenses/BSD-3-Clause)
![C++](https://img.shields.io/badge/C++-11+-1069ac?logo=c%2B%2B)
![C](https://img.shields.io/badge/C-99+-7991b5?logo=c&logoColor=eee)

Tested on:

[![MacOS](https://img.shields.io/badge/MacOS-13.7%20%7C%2010.14%20%7C%2010.15-999999?logo=apple&logoColor=eee)](https://www.apple.com/)
[![Ubuntu](https://img.shields.io/badge/Ubuntu-16.04%20%7C%2018.04%20LTS-e95420?logo=ubuntu&logoColor=eee)](https://www.ubuntu.com/)
[![CentOS](https://img.shields.io/badge/CentOS-7%20%7C%208-262577?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA5MSA5MC41Ij48cGF0aCBmaWxsPSIjYjRiZWMxIiBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMi4wMDIiIGQ9Ik00NS44MSA0MS4yMWwtMTkuOC0xOS45IDE5LjgtMTkuOSAxOS44IDE5Ljl6Ii8+PHBhdGggZmlsbD0iIzM5NGQ1NCIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjIuMDAyIiBkPSJNNjkuODEgNjUuMjFsLTE5LjgtMTkuOSAxOS44LTE5LjkgMTkuOCAxOS45eiIvPjxwYXRoIGZpbGw9IiNmY2U5NGYiIGQ9Ik0xNC43MSAxNC4xMWgyN3YyNi45aC0yN3oiLz48cGF0aCBmaWxsPSIjYjRiZWMxIiBkPSJNNDkuNzEgNDkuMjFoMjYuOHYyN2gtMjYuOHoiLz48cGF0aCBmaWxsPSIjNGY2YTc0IiBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMi4wMDIiIGQ9Ik00MS42MSA0NS4zMWwtMjAuMSAxOS44LTIwLjEtMTkuOCAyMC4xLTE5Ljh6Ii8+PHBhdGggZmlsbD0iIzRmNmE3NCIgZD0iTTQ5LjgxIDE0LjMxaDI2Ljl2MjYuOWgtMjYuOXoiLz48ZyBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMi4wMDIiPjxwYXRoIGZpbGw9IiM4YTlhYTAiIGQ9Ik02NS42MSA2OS4zMWwtMjAuMSAxOS44LTIwLjEtMTkuOCAyMC4xLTE5Ljh6TTEzLjcxIDEzLjIxaDI4djI3LjloLTI4eiIvPjxwYXRoIGZpbGw9Im5vbmUiIGQ9Ik00OS43MSA0OS4yMWgyNy44djI4aC0yNy44ek00OS44MSAxMy4zMWgyNy45djI3LjloLTI3Ljl6Ii8+PHBhdGggZmlsbD0iIzM5NGQ1NCIgZD0iTTEzLjcxIDQ5LjIxaDI4djI4aC0yOHoiLz48cGF0aCBmaWxsPSJub25lIiBkPSJNNDEuNjEgNDUuMzFsLTIwLjEgMTkuOC0yMC4xLTE5LjggMjAuMS0xOS44ek00NS44MSA0MS4yMWwtMTkuOC0xOS45IDE5LjgtMTkuOSAxOS44IDE5Ljl6TTY1LjYxIDY5LjMxbC0yMC4xIDE5LjgtMjAuMS0xOS44IDIwLjEtMTkuOHpNNjkuODEgNjUuMjFsLTE5LjgtMTkuOSAxOS44LTE5LjkgMTkuOCAxOS45eiIvPjwvZz48L3N2Zz4=)](https://www.centos.org/)

About sciLoc
----------

sciLoc is a SeisComp plugin for iLoc, a single-event earthquake location
algorithm that tackles the problem that raypaths traversing unmodeled 3D
velocity heterogeneities in the Earth cause correlated travel-time predictions
that can result in location bias and underestimated location uncertainties.

### A brief history of sciLoc

The original algorithm was developed under an AFRL contract (Bondár and McLaughlin, 2009a).
Operational location algorithm at the ISC since 2010 (Bondár and Storchak, 2011) and at the EMSC since 2019, (Steed et al., 2019; Bondár et al., 2020).
SeisComp plugin since 2020 (Weber et al., 2019).

### sciLoc in a nutshell

- SeisComp plugin for iLoc with a SeisComp GUI.
- Assumes that an event is already formed, and the phases associated to the event.
- Accounts for correlated travel-time prediction errors.
- Initial hypocenter guess from Neighbourhood Algorithm search (Sambridge, 1999; Sambridge and Kennett, 2001).
- Linearised inversion using a priori estimate of the full data covariance matrix (Bondár and McLaughlin, 2009a).
- Attempts for free-depth solution only if there is depth resolution, otherwise sets default depth from a global grid of reliable free depth events from historical seismicity.
- Uses seismic, hydroacoustic and infrasound observations.
- Arrival time, slowness and azimuth measurements are used in the location.
- Identifies seismic phases w.r.t the initial guess, then the best hypocenter estimate from the NA grid search.
- Uses all valid ak135 (Kennett et al., 1995) phases in location.
- Elevation and ellipticity corrections (Dziewonski and Gilbert, 1976; Kennett and Gudmundsson, 1996).
- Depth-phase bounce point corrections (Engdahl et al., 1998).
- Uses RSTT travel-time predictions for Pn/Sn and Pg/Lg (Myers et al., 2010).
- RSTT provides its own uncertainty estimates (Begnaud et al., 2020, 2021).
- Optional use of a local velocity model.
- Predictions for local phases are calculated up to 3 degrees, beyond that iLoc switches to RSTT/ak135 predictions.
- Local phase travel time predictions are calculated for Pg/Sg, Pb/Sb, Pn/Sn.
- Performs both the Bondár and McLaughlin (2009b) and the Gallacher et al. (2025) ground truth candidate event selection test for relocated events.

### Citations

- Bondár, I. and K. McLaughlin (2009). Seismic location bias and uncertainty in the presence of correlated and non-Gaussian travel-time errors, Bull. Seism. Soc. Am., 99, 172-193, https://doi.org/10.1785/0120080922.
- Bondár, I., and D. Storchak (2011). Improved location procedures at the International Seismological Centre, Geophys. J. Int., 186, 1220-1244, https://doi.org/10.1111/j.1365-246X.2011.05107.x.
- Weber B., I. Bondár, D. Rößler, J. Becker, SeisComp3 iLoc Integration Applied to Array Processing, CTBT: Science and Technology Conference, Book of Abstracts, T3.5.P54, p.187, Vienna, 24-28 June, 2019.

Compiling iLoc from source
------------------------------

### Environment variables

The environment variable `$ILOCROOT` is required by the *Makefiles* during
the compilation process. Set it to your iLoc directory in your `.bashrc` file:

export ILOCROOT=$HOME/iLoc


### Dependencies


#### MacOS

Below is a list of packages and software required to build iLoc from source.
Xcode or Command Line Tools need to be installed.

Software  | Purpose
:---------|:------------------------------
make      | Running compile scripts
gcc       | Compile C code
g++       | Build GeoTess and core RSTT libraries


#### Linux

Below is a list of packages and software required to build RSTT from source.

Software  | Purpose
:---------|:------------------------------
make      | Running compile scripts
gcc       | Build C library and tests
g++       | Build GeoTess and core RSTT libraries and tests
lapack    | Build lapack and blas libraries

The easiest way to satisfy these dependencies on Linux is, depending on your
distribution and package manager, by running one of these sets of commands
in a termal window:

```bash
# C++, C
$ sudo apt install build-essentials
$ sudo apt install g++-multilib

# lapack and blas libraries
$ sudo apt install liblapack3
```


```bash
# C++, C
$ sudo yum install kernel-devel gcc gcc-c++
$ sudo yum install g++-multilib

# lapack and blas libraries
$ sudo yum groupinstall 'Development Tools' -y
$ sudo yum install lapack.x86_64 -y
```


Contact Information
-------------------

For questions/issues/comments about the software, please contact:

* [Istvan Bondar][slsiloc] (Seismic Location Services)
* [Dirk Roessler][gempa] (Gempa GmbH)

[iloc]:     https://www.slsiloc.eu
[seiscomp]:    https://www.gempa.de
[slsiloc]:  mailto:istvan.bondar@slsiloc.eu
[gempa]:   roessler@gempa.de
