# Antineutrino Analysis in the MINERvA's Target region focused on neutrons in the final state.
AS OF 03/23/2022, THIS HAS BEEN MOVED TO THE MinervaExpt ORGANIZATION ON GITHUB AND THIS IS ARCHIVAL. Code for David Last's thesis analysis which currently looks for CC muon anti-neutrino events with at least one well-identified neutron in the final state in the passive nuclear targets in MINERvA. Uses the MINERvA Analysis Toolkit (MAT).  Starts from slightly-modified (longterm plan is to not have these modifications be separate, but some careful handling of merging them in is needed) MINERvA data preservation "anaTuple" .root files, produced by MasterAnaDev (MAD), and produces all histograms needed to extract a cross section.  The code is based on the tutorial from which it is forked. See [the wiki](https://github.com/MinervaExpt/MINERvA-101-Cross-Section/wiki) for more information on the tutorial.

## Installation Instructions
**First** try [the wiki instructions](https://github.com/MinervaExpt/MINERvA-101-Cross-Section/wiki/Installation/#installing_the_whole_minerva_101_2021_tutorial) with `dlast44/Target-Antinu-Neutrons.git` replacing instances of `MinervaExpt/MINERvA-101-Cross-Section.git` for a comprehensive guide to the dependencies and this package.  If that file no longer exists, you'll have to install the [dependencies](#Dependencies) yourself and follow the instrutions below.

**In case the wiki instructions aren't available:**
1. Install [dependencies](#Dependencies)
2. Make a working directory.  Mine is called `TargetNeutronsAna/`
3. Download the complete source from [github](https://github.com/dlast44/Target-Antinu-Neutrons): `git clone git@github.com:dlast44/Target-Antinu-Neutrons.git #Makes a src subdirectory for you`
4. Make a build directory: `mkdir opt && cd opt && mkdir build && cd build #opt for optimized build as opposed to debug build`
5. Run cmake to generate a build system: ``cmake ../../Target-Antinu-Neutrons/bootstrap -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release``
   If `cmake` fails talking about dependencies, you may have to manually point to them with command line arguments like `-DPlotUtils_DIR=/path/to/PlotUtilsInstallPrefix/lib`
6. Compile and install: `make install`

## Usage
```
runEventLoop <dataPlaylist.txt> <mcPlaylist.txt>

*** Explanation ***
Reduce MasterAnaDev AnaTuples to event selection histograms to extract a
single-differential inclusive cross section for the 2021 MINERvA 101 tutorial.

*** The Input Files ***
Playlist files are plaintext files with 1 file name per line.  Filenames may be
xrootd URLs or refer to the local filesystem.  The first playlist file's
entries will be treated like data, and the second playlist's entries must
have the "Truth" tree to use for calculating the efficiency denominator.

*** Output ***
Produces a single runEventLoop.root file with all histograms needed for the
ExtractCrossSection program also built by this package.  You'll need a
.rootlogon.C that loads ROOT object definitions from PlotUtils to access
systematics information from these files.

*** Return Codes ***
0 indicates success.  All histograms are valid only in this case.  Any other
return code indicates that histograms should not be used.  Error messages
about what went wrong will be printed to stderr.  So, they'll end up in your
terminal, but you can separate them from everything else with something like:
"runEventLoop data.txt mc.txt 2> errors.txt"
```

## Dependencies
0. [git](https://git-scm.com/downloads): version control system.  **You probably already have this**
1. [CMake 2.8](https://cmake.org/install/): build system generator for lots of operating systems.  **You probably already have this**
2. [ROOT](https://root.cern/install/): object-oriented toolkit for high energy physics analysis.  Make sure to enable at least xrootd, kerberos, and Minuit for this tutorial.  **Already installed on MINERvA GPVMs**
3. [MAT](https://github.com/MinervaExpt/MAT): the MINERvA Analysis Toolkit
4. [MAT-MINERvA](https://github.com/MinervaExpt/MAT-MINERvA): MINERvA-specific systematics and other plugins to the MAT
5. [UnfoldUtils](https://github.com/MinervaExpt/UnfoldUtils): MINERvA's fork of [RooUnfold](https://gitlab.cern.ch/RooUnfold/RooUnfold) with a compatibility layer for embedding systematics in histograms and MINERvA-specific systematics tweaks (TODO: via github)
6. [GENIEXSecExtract](https://github.com/MinervaExpt/GENIEXSecExtract): MINERvA's closure test procedure
7. [MParamFiles](https://cdcvs.fnal.gov/redmine/projects/minerva-sw/repository/show/AnalysisFramework/MParamFiles): additional reweighting and calibration parameters for MINERvA systematics
8. [MATFluxAndReweightFiles](https://cdcvs.fnal.gov/redmine/projects/minerva-sw/repository/show/AnalysisFramework/Ana/MATFluxAndReweightFiles): the flux incident on the MINERvA detector during different run periods and some more reweight configurations for systematic uncertainties
