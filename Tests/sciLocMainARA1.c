/*
 * Copyright (c) 2018, Istvan Bondar,
 * Written by Istvan Bondar, ibondar2014@gmail.com
 *
 * BSD Open Source License.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * sciLoc
 *
 * Istvan Bondar
 * ibondar2014@gmail.com
 *
 * SeisComp3 single-event location for with RSTT travel-time predictions
 *     based on ISC and iLoc location algorithms
 *     (Bondar and Storchak, 2011; Bondar et al. 2019)
 *     supports RSTT (Myers et al., 2010) global upper mantle model
 *         RSTT regional/local travel-time predictions (Pg/Lg/Pn/Sn)
 *     supports user-provided local velocity models
 *         local phase travel-time predictions (Pg/Pb/Pn/Lg/Sg/Sb/Sn)
 *     supports ak135 global travel-time predictions
 *         ellipticity (Dziewonski and Gilbert, 1976, Kennett and Gudmundsson,
 *             1996) and elevation corrections
 *         bounce point corrections for depth phases as well as water column
 *             correction for pwP (Engdahl et al., 1998)
 *     supports infrasound and hydroacoustic azimuth measurements in location
 *     neighbourhood algorithm (NA) grid-search for a refined initial
 *         hypocentre guess (Sambridge, 1999; Sambridge and Kennett, 2001)
 *     iterative linearized least-squares inversion that accounts for
 *         correlated model errors (Bondar and McLaughlin, 2009)
 *     nearest-neighbour ordering of stations to block-diagonalize data
 *         covariance matrix (de Hoon et al., 2004)
 *     depth-phase stacking for robust depth estimation
 *         (Murphy and Barker, 2006)
 *     free-depth solution is attempted only if there is depth resolution
 *         (number of defining depth phases > MinDepthPhases and
 *          number of agencies reporting depth phases >= MindDepthPhaseAgencies)
 *          or number of local defining stations >= MinLocalStations
 *          or number of defining P and S pairs >= MinSPpairs
 *          or (number of defining PcP, ScS >= MinCorePhases and
 *              number of agencies reporting core reflections >= MindDepthPhaseAgencies)
 *
 * References
 *   Bondár, I. and K. McLaughlin, (2009).
 *      Seismic location bias and uncertainty in the presence of correlated and
 *      non-Gaussian travel-time errors,
 *      Bull. Seism. Soc. Am., 99, 172-193, doi:10.1785/0120080922.
 *   Bondár, I. and D. Storchak, (2011).
 *     Improved location procedures at the International Seismological Centre,
 *     Geophys. J. Int., 186, 1220-1244, doi:10.1111/j.1365-246X.2011.05107.x
 *   Bondár, I., E.R. Engdahl, A. Villasenor, J.Harris and D. Storchak (2015).
 *      ISC-GEM: Global instrumental earthquake catalogue (1900-2009),
 *      II. Location and seismicity patterns,
 *      Phys. Earth. Planet. Int., doi: 10.1016/j.pepi.2014.06.002, 239, 2-13.
 *   Bondár, I. P. Mónus, Cs. Czanik, M. Kiszely, Z. Gráczer, Z. Wéber,
 *     and the AlpArrayWorking Group, (2019).
 *     Relocation of Seismicity in the Pannonian Basin Using a Global 3D
 *     Velocity Model,
 *     Seism. Res. Let. doi:10.1785/0220180143, 2019.
 *
 * Input files in the auxiliary data directory
 *     Configuration parameters
 *         $ILOCROOT/auxdata/sciLocpars/Config.txt
 *     Phase specific parameters
 *         $ILOCROOT/auxdata/sciLocpars/PhaseConfig.txt
 *     Default RSTT model
 *         $ILOCROOT/auxdata/RSTTmodel/rstt201404um.geotess
 *     Travel-time tables
 *         $ILOCROOT/auxdata/ak135/[*].tab
 *     Ellipticity correction coefficients for ak135
 *         $ILOCROOT/auxdata/ak135/ELCOR.dat
 *     Topography file for bounce point corrections
 *         $ILOCROOT/auxdata/topo/etopo5_bed_g_i2.bin
 *         ETOPO1 resampled to 5'x5' resolution
 *     Flinn-Engdahl regionalization (1995)
 *         $ILOCROOT/auxdata/FlinnEngdahl/FE.dat
 *     Default depth grid
 *         $ILOCROOT/auxdata/FlinnEngdahl/DefaultDepth0.5.grid
 *     Default depth for Flinn-Engdahl regions
 *         $ILOCROOT/auxdata/FlinnEngdahl/GRNDefaultDepth.ak135.dat
 *         For locations where no default depth grid point exists.
 *     Generic variogram model for data covariance matrix
 *         $ILOCROOT/auxdata/variogram/variogram.model
 *
 *
 * For each event try 2 options:
 *     0 - free depth (if there is depth resolution)
 *     1 - fix depth to region-dependent default depth (if option 0 fails)
 * If convergence is reached, no further options are tried.
 *
 */

#include "sciLocInterface.h"

/*
 *
 * main body (should receive iLocConfig, Hypocenter, Assocs and StaLocs)
 *
 */
int main()
{
    ILOC_CONF iLocConfig;
    ILOC_PHASEIDINFO PhaseIdInfo;
    ILOC_TTINFO TTInfo;              /* global velocity model info and phase list */
    ILOC_TTINFO LocalTTInfo;          /* local velocity model info and phase list */
    ILOC_TT_TABLE *TTtables = (ILOC_TT_TABLE *)NULL;      /* global travel-time tables */
    ILOC_TT_TABLE *LocalTTtables = (ILOC_TT_TABLE *)NULL;  /* local travel-time tables */
    ILOC_EC_COEF *ec = (ILOC_EC_COEF *)NULL;    /* ellipticity correction coefficients */
    ILOC_VARIOGRAM variogram;                          /* generic variogram model */
    ILOC_FE fe;                        /* Flinn-Engdahl geographic region numbers */
    ILOC_DEFAULTDEPTH DefaultDepth;         /* default depth and etopo structures */
    ILOC_HYPO Hypocenter;
    ILOC_ASSOC *Assocs = (ILOC_ASSOC *)NULL;;
    ILOC_STA *StaLocs = (ILOC_STA *)NULL;
    int retval = ILOC_UNKNOWN_ERROR;
/*
 *  set timezone to get epoch times right and other inits
 */
    setenv("TZ", "", 1);
    tzset();
/*
 *
 *  set defaults for iLocConfig for tests
 *
 */
/*
 *  directory of auxiliary data files
 */
    strcpy(iLocConfig.auxdir, "/Users/istvanbondar/iLoc4.0/auxdata");
    iLocConfig.Verbose = 2;
/*
 *  Travel time predictions
 */
    strcpy(iLocConfig.TTmodel, "ak135");
    iLocConfig.UseRSTT = 1;
    strcpy(iLocConfig.RSTTmodel, "/Users/istvanbondar/iLoc4.0/auxdata/RSTTmodels/pdu202009Du.geotess");
    iLocConfig.UseRSTTPnSn = 1;
    iLocConfig.UseRSTTPgLg = 1;
    strcpy(iLocConfig.LocalVmodel, "");
    iLocConfig.UseLocalTT = 0;
    iLocConfig.MaxLocalTTDelta = 3.;
/*
 *  ETOPO parameters
 */
    strcpy(iLocConfig.EtopoFile, "etopo5_bed_g_i2.bin");
    iLocConfig.EtopoNlon = 4321;
    iLocConfig.EtopoNlat = 2161;
    iLocConfig.EtopoRes = 0.0833333;
/*
 *  Linearized inversion
 */
    iLocConfig.MinIterations = 4;
    iLocConfig.MaxIterations = 100;
    iLocConfig.MinNdefPhases = 4;
    iLocConfig.SigmaThreshold = 6.;
    iLocConfig.DoCorrelatedErrors = 1;
    iLocConfig.AllowDamping = 1;
    iLocConfig.DoNotRenamePhases = 0;
/*
 *  depth resolution
 */
    iLocConfig.MaxLocalDistDeg = 0.2;
    iLocConfig.MinLocalStations = 1;
    iLocConfig.MaxSPDistDeg = 2.;
    iLocConfig.MinSPpairs = 3;
    iLocConfig.MinCorePhases = 3;
    iLocConfig.MinDepthPhases = 3;
    iLocConfig.MaxShallowDepthError = 30.;
    iLocConfig.MaxDeepDepthError = 60.;
/*
 *  NA search parameters
 */
    iLocConfig.DoGridSearch = 1;
    iLocConfig.NAsearchRadius = 5.;
    iLocConfig.NAsearchDepth = 300.;
    iLocConfig.NAsearchOT = 30.;
    iLocConfig.NAlpNorm = 1.;
    iLocConfig.NAiterMax = 5;
    iLocConfig.NAinitialSample = 1500;
    iLocConfig.NAnextSample = 150;
    iLocConfig.NAcells = 25;
/*
 *
 * set Hypocenter structure for test purposes
 *    slow converging Ukrainian ammunition explosion
 *
 *  anthropogenic event, depth will be fixed to zero!
 */
    Hypocenter.isManMade = 1;
    Hypocenter.numSta = 8;
    Hypocenter.numPhase = 8;
    Hypocenter.Time = 1510753873.27;
    Hypocenter.Lat = -46.1107;
    Hypocenter.Lon = -59.6923;
    Hypocenter.Depth = 0.6;
    Hypocenter.FixOT = 0;
    Hypocenter.FixLat = 0;
    Hypocenter.FixLon = 0;
    Hypocenter.FixDepth = 1;
    Hypocenter.FixHypo = 0;
/*
 *
 *  memory allocation for stations and phases
 *
 */
    Assocs = (ILOC_ASSOC *)calloc(Hypocenter.numPhase, sizeof(ILOC_ASSOC));
    if ((StaLocs = (ILOC_STA *)calloc(Hypocenter.numSta, sizeof(ILOC_STA))) == NULL) {
        iLoc_Free(Assocs);
        fprintf(stderr, "Main: cannot allocate memory\n");
        return ILOC_MEMORY_ALLOCATION_ERROR;
    }
/*
 *  set StaLocs structure
 */
    StaLocs[0].StaLat = -8.95270;    // H10S3
    StaLocs[0].StaLon = -14.66290;
    StaLocs[0].StaElevation = -1703.0;
    StaLocs[1].StaLat = -8.95910;    // H10S2
    StaLocs[1].StaLon = -14.64530;
    StaLocs[1].StaElevation = -1757.0;
    StaLocs[2].StaLat = -7.84090;    // H10N3
    StaLocs[2].StaLon = -14.50170;
    StaLocs[2].StaElevation = -2041.0;
    StaLocs[3].StaLat = -7.84570;    // H10N1
    StaLocs[3].StaLon = -14.48020;
    StaLocs[3].StaElevation = -1925.0;
    StaLocs[4].StaLat = -7.82780;    // H10N2
    StaLocs[4].StaLon = -14.48750;
    StaLocs[4].StaElevation = -2049.0;
    StaLocs[5].StaLat = -46.85360;    // H04S2
    StaLocs[5].StaLon = 51.89350;
    StaLocs[5].StaElevation = -1225.0;
    StaLocs[6].StaLat = -46.83540;    // H04S3
    StaLocs[6].StaLon = 51.88700;
    StaLocs[6].StaElevation = -1175.0;
    StaLocs[7].StaLat = -46.84090;    // H04S1
    StaLocs[7].StaLon = 51.91300;
    StaLocs[7].StaElevation = -1125.0;
/*
 *  set Assocs structure
 */
    Assocs[0].StaInd = 0;    // H10S3
    Assocs[0].arid = 1;
    strcpy(Assocs[0].PhaseHint, "H");
    Assocs[0].phaseFixed = 0;
    Assocs[0].ArrivalTime = 1510757883.423;
    Assocs[0].BackAzimuth = 218.5;
    Assocs[0].Slowness = 73.60;
    Assocs[0].Timedef = 1;
    Assocs[0].Azimdef = 1;
    Assocs[0].Slowdef = 1;
    Assocs[0].Deltim = ILOC_NULLVAL;
    Assocs[1].StaInd = 1;    // H10S2
    Assocs[1].arid = 2;
    strcpy(Assocs[1].PhaseHint, "H");
    Assocs[1].phaseFixed = 0;
    Assocs[1].ArrivalTime = 1510757883.76;
    Assocs[1].BackAzimuth = 1.5;
    Assocs[1].Slowness = 73.60;
    Assocs[1].Timedef = 1;
    Assocs[1].Azimdef = 1;
    Assocs[1].Slowdef = 1;
    Assocs[1].Deltim = ILOC_NULLVAL;
    Assocs[2].StaInd = 2;    // H10N3
    Assocs[2].arid = 3;
    strcpy(Assocs[2].PhaseHint, "H");
    Assocs[2].phaseFixed = 0;
    Assocs[2].ArrivalTime = 1510757956.191;
    Assocs[2].BackAzimuth = 217.0;
    Assocs[2].Slowness = 75.40;
    Assocs[2].Timedef = 1;
    Assocs[2].Azimdef = 1;
    Assocs[2].Slowdef = 1;
    Assocs[2].Deltim = ILOC_NULLVAL;
    Assocs[3].StaInd = 3;    // H10N1
    Assocs[3].arid = 4;
    strcpy(Assocs[3].PhaseHint, "H");
    Assocs[3].phaseFixed = 0;
    Assocs[3].ArrivalTime = 1510757956.878;
    Assocs[3].BackAzimuth = 217.0;
    Assocs[3].Slowness = 75.40;
    Assocs[3].Timedef = 1;
    Assocs[3].Azimdef = 1;
    Assocs[3].Slowdef = 1;
    Assocs[3].Deltim = ILOC_NULLVAL;
    Assocs[4].StaInd = 4;    // H10N2
    Assocs[4].arid = 5;
    strcpy(Assocs[4].PhaseHint, "H");
    Assocs[4].phaseFixed = 0;
    Assocs[4].ArrivalTime = 1510757957.619;
    Assocs[4].BackAzimuth = 217.0;
    Assocs[4].Slowness = 75.40;
    Assocs[4].Timedef = 1;
    Assocs[4].Azimdef = 1;
    Assocs[4].Slowdef = 1;
    Assocs[4].Deltim = ILOC_NULLVAL;
    Assocs[5].StaInd = 5;    // H04S2
    Assocs[5].arid = 6;
    strcpy(Assocs[5].PhaseHint, "H");
    Assocs[5].phaseFixed = 0;
    Assocs[5].ArrivalTime = 1510759176.373;
    Assocs[5].BackAzimuth = 223.5;
    Assocs[5].Slowness = 75.40;
    Assocs[5].Timedef = 1;
    Assocs[5].Azimdef = 1;
    Assocs[5].Slowdef = 1;
    Assocs[5].Deltim = ILOC_NULLVAL;
    Assocs[6].StaInd = 6;    // H04S3
    Assocs[6].arid = 7;
    strcpy(Assocs[6].PhaseHint, "H");
    Assocs[6].phaseFixed = 0;
    Assocs[6].ArrivalTime = 1510759177.622;
    Assocs[6].BackAzimuth = 223.5;
    Assocs[6].Slowness = 75.40;
    Assocs[6].Timedef = 1;
    Assocs[6].Azimdef = 1;
    Assocs[6].Slowdef = 1;
    Assocs[6].Deltim = ILOC_NULLVAL;
    Assocs[7].StaInd = 7;    // H04S1
    Assocs[7].arid = 8;
    strcpy(Assocs[7].PhaseHint, "H");
    Assocs[7].phaseFixed = 0;
    Assocs[7].ArrivalTime = 1510759177.587;
    Assocs[7].BackAzimuth = 223.5;
    Assocs[7].Slowness = 75.40;
    Assocs[7].Timedef = 1;
    Assocs[7].Azimdef = 1;
    Assocs[7].Slowdef = 1;
    Assocs[7].Deltim = ILOC_NULLVAL;
/*
 *  end of test event
 */


/*
 *
 *  Read data files from iLocConfig.auxdir
 *     iLocConfig->auxdir/iLocpars/IASPEIPhaseMap.txt
 *     iLocConfig->auxdir/iLocpars/PhaseConfig.iLocConfig->TTmodel.txt
 *     iLocConfig->auxdir/iLocConfig->TTmodel/iLocConfig->TTmodel.*.tab
 *     iLocConfig->auxdir/iLocConfig->TTmodel/ELCOR.dat
 *     iLocConfig->auxdir/FlinnEngdahl/FE.dat
 *     iLocConfig->auxdir/FlinnEngdahl/DefaultDepth0.5.grid
 *     iLocConfig->auxdir/FlinnEngdahl/GRNDefaultDepth.iLocConfig->TTmodel.dat
 *     iLocConfig->auxdir/topo/etopo5_bed_g_i2.bin
 *     iLocConfig->auxdir/variogram/variogram.dat
 *     iLocConfig->auxdir/iLocpars/PhaseConfig.iLocConfig->LocalVmodel.txt
 *     iLocConfig->auxdir/localmodels/iLocConfig->LocalVmodel.localmodel.dat
 *     iLocConfig->RSTTmodel
 *
 */
    if (iLoc_ReadAuxDataFiles(&iLocConfig, &PhaseIdInfo, &fe, &DefaultDepth,
                              &variogram, &TTInfo, &TTtables, &ec,
                              &LocalTTInfo, &LocalTTtables)) {
        fprintf(stderr, "Couldn't read iLoc auxiliary data files!\n");
        return ILOC_FAILURE;
    }
/*
 *
 *  Locate event
 *     iLocConfig - iLoc control parameters
 *     Hypocenter - initial guess on input, final solution on output
 *     Assocs     - arrival time, slowness azimuth, phasehints on input,
 *                  phaseids, residuals on output
 *     StationLocs - Station codes and coordinates
 */
    retval =iLoc_Locator(&iLocConfig, &PhaseIdInfo, &fe, &DefaultDepth,
                         &variogram, ec, &TTInfo, TTtables, &LocalTTInfo,
                         LocalTTtables, &Hypocenter, Assocs, StaLocs);
    if (retval)
        fprintf(stderr, "CAUTION: No solution found!\n");
/*
 *
 *  Free allocated memory
 *
 */
    iLoc_FreeAuxData(&PhaseIdInfo, &fe, &DefaultDepth, &variogram, &TTInfo,
            TTtables, ec, &LocalTTInfo, LocalTTtables, iLocConfig.UseRSTT);
    iLoc_Free(Assocs);
    iLoc_Free(StaLocs);
    return retval;
}

