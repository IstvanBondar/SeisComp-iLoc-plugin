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
    iLocConfig.Verbose = 1;
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
    iLocConfig.NAiterMax = 4;
    iLocConfig.NAinitialSample = 1500;
    iLocConfig.NAnextSample = 150;
    iLocConfig.NAcells = 25;
/*
 *
 * set Hypocenter structure for test purposes
 *
 */
    Hypocenter.isManMade = 0;
    Hypocenter.numSta = 140;
    Hypocenter.numPhase = 253;
    Hypocenter.Time = 1262324232.75;
    Hypocenter.Lat = 13.720;
    Hypocenter.Lon = 125.625;
    Hypocenter.Depth = 17.1;
    Hypocenter.FixOT = 0;
    Hypocenter.FixLat = 0;
    Hypocenter.FixLon = 0;
    Hypocenter.FixDepth = 0;
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
    StaLocs[0].StaLat = 13.59600;    // PVCP
    StaLocs[0].StaLon = 124.15400;
    StaLocs[0].StaElevation = 50.0;
    StaLocs[1].StaLat = 12.51050;    // CNP
    StaLocs[1].StaLon = 124.66317;
    StaLocs[1].StaElevation = 0;
    StaLocs[2].StaLat = 11.16500;    // PLP
    StaLocs[2].StaLon = 124.97900;
    StaLocs[2].StaElevation = 133.0;
    StaLocs[3].StaLat = 11.05100;    // OCLP
    StaLocs[3].StaLon = 124.60900;
    StaLocs[3].StaElevation = 130.0;
    StaLocs[4].StaLat = 13.32300;    // AUQP
    StaLocs[4].StaLon = 122.67500;
    StaLocs[4].StaElevation = 50.0;
    StaLocs[5].StaLat = 13.90500;    // GQP
    StaLocs[5].StaLon = 122.44600;
    StaLocs[5].StaElevation = 50.0;
    StaLocs[6].StaLat = 11.56000;    // RCP
    StaLocs[6].StaLon = 122.74500;
    StaLocs[6].StaElevation = 140.0;
    StaLocs[7].StaLat = 10.13400;    // MSLP
    StaLocs[7].StaLon = 124.85900;
    StaLocs[7].StaElevation = 50.0;
    StaLocs[8].StaLat = 13.45900;    // BOAC
    StaLocs[8].StaLon = 121.84400;
    StaLocs[8].StaElevation = 29.0;
    StaLocs[9].StaLat = 10.31700;    // LLP
    StaLocs[9].StaLon = 123.96600;
    StaLocs[9].StaElevation = 135.0;
    StaLocs[10].StaLat = 10.62600;    // GUIM
    StaLocs[10].StaLon = 122.58900;
    StaLocs[10].StaElevation = 50.0;
    StaLocs[11].StaLat = 17.06300;    // PALP
    StaLocs[11].StaLon = 122.42500;
    StaLocs[11].StaElevation = 50.0;
    StaLocs[12].StaLat = 14.10219;    // TGY
    StaLocs[12].StaLon = 120.93669;
    StaLocs[12].StaElevation = 600.0;
    StaLocs[13].StaLat = 15.56200;    // PCPH
    StaLocs[13].StaLon = 121.09600;
    StaLocs[13].StaElevation = 7.0;
    StaLocs[14].StaLat = 17.70300;    // CVP
    StaLocs[14].StaLon = 121.82200;
    StaLocs[14].StaElevation = 60.0;
    StaLocs[15].StaLat = 12.00300;    // BUSP
    StaLocs[15].StaLon = 120.19900;
    StaLocs[15].StaElevation = 66.0;
    StaLocs[16].StaLat = 17.55100;    // SZP
    StaLocs[16].StaLon = 120.45600;
    StaLocs[16].StaElevation = 50.0;
    StaLocs[17].StaLat = 7.07000;    // DAV
    StaLocs[17].StaLon = 125.57900;
    StaLocs[17].StaElevation = 145.7;
    StaLocs[18].StaLat = 23.39239;    // YULB
    StaLocs[18].StaLon = 121.29731;
    StaLocs[18].StaElevation = 295.0;
    StaLocs[19].StaLat = 26.83217;    // JOW
    StaLocs[19].StaLon = 128.27450;
    StaLocs[19].StaElevation = 220.0;
    StaLocs[20].StaLat = 0.47710;    // MRSI
    StaLocs[20].StaLon = 121.94060;
    StaLocs[20].StaElevation = 0.0;
    StaLocs[21].StaLat = -0.91080;    // APSI
    StaLocs[21].StaLon = 121.64870;
    StaLocs[21].StaElevation = 0.0;
    StaLocs[22].StaLat = -0.86300;    // SWI
    StaLocs[22].StaLon = 131.25980;
    StaLocs[22].StaElevation = 0.0;
    StaLocs[23].StaLat = 19.02940;    // QIZ
    StaLocs[23].StaLon = 109.84300;
    StaLocs[23].StaElevation = 230.0;
    StaLocs[24].StaLat = -3.23900;    // NLAI
    StaLocs[24].StaLon = 127.09980;
    StaLocs[24].StaElevation = 0.0;
    StaLocs[25].StaLat = -3.34610;    // MSAI
    StaLocs[25].StaLon = 128.92850;
    StaLocs[25].StaElevation = 82.0;
    StaLocs[26].StaLat = 2.45000;    // SBUM
    StaLocs[26].StaLon = 112.22000;
    StaLocs[26].StaElevation = 30.9;
    StaLocs[27].StaLat = -3.04510;    // TTSI
    StaLocs[27].StaLon = 119.81900;
    StaLocs[27].StaElevation = 0.0;
    StaLocs[28].StaLat = -3.96460;    // SPSI
    StaLocs[28].StaLon = 119.76910;
    StaLocs[28].StaElevation = 0.0;
    StaLocs[29].StaLat = -4.40050;    // BNSI
    StaLocs[29].StaLon = 120.10650;
    StaLocs[29].StaElevation = 0.0;
    StaLocs[30].StaLat = 32.05170;    // NJ2
    StaLocs[30].StaLon = 118.85400;
    StaLocs[30].StaElevation = 45.0;
    StaLocs[31].StaLat = 33.12167;    // JNU
    StaLocs[31].StaLon = 130.87833;
    StaLocs[31].StaElevation = 540.0;
    StaLocs[32].StaLat = -3.46250;    // BBKI
    StaLocs[32].StaLon = 114.84110;
    StaLocs[32].StaElevation = 110.0;
    StaLocs[33].StaLat = 26.45860;    // GYA
    StaLocs[33].StaLon = 106.66400;
    StaLocs[33].StaElevation = 1162.0;
    StaLocs[34].StaLat = 21.33383;    // SLVN
    StaLocs[34].StaLon = 103.90500;
    StaLocs[34].StaElevation = 700.0;
    StaLocs[35].StaLat = 30.27180;    // ENH
    StaLocs[35].StaLon = 109.48700;
    StaLocs[35].StaElevation = 487.0;
    StaLocs[36].StaLat = 37.44211;    // KSAR
    StaLocs[36].StaLon = 127.88439;
    StaLocs[36].StaElevation = 109.0;
    StaLocs[37].StaLat = 37.45400;    // KSRS
    StaLocs[37].StaLon = 127.92300;
    StaLocs[37].StaElevation = 0;
    StaLocs[38].StaLat = 25.12333;    // KMI
    StaLocs[38].StaLon = 102.74000;
    StaLocs[38].StaElevation = 1940.0;
    StaLocs[39].StaLat = -8.47020;    // JAGI
    StaLocs[39].StaLon = 114.15210;
    StaLocs[39].StaElevation = 171.0;
    StaLocs[40].StaLat = 34.03944;    // XAN
    StaLocs[40].StaLon = 108.92139;
    StaLocs[40].StaElevation = 630.0;
    StaLocs[41].StaLat = 36.54170;    // MJAR
    StaLocs[41].StaLon = 138.20900;
    StaLocs[41].StaElevation = 422.0;
    StaLocs[42].StaLat = 18.45750;    // CM31
    StaLocs[42].StaLon = 98.94289;
    StaLocs[42].StaElevation = 306.6;
    StaLocs[43].StaLat = 18.45750;    // CMAR
    StaLocs[43].StaLon = 98.94290;
    StaLocs[43].StaElevation = 307.0;
    StaLocs[44].StaLat = 30.91000;    // CD2
    StaLocs[44].StaLon = 103.75800;
    StaLocs[44].StaElevation = 628.0;
    StaLocs[45].StaLat = 40.04030;    // BJI
    StaLocs[45].StaLon = 116.17500;
    StaLocs[45].StaElevation = 43.0;
    StaLocs[46].StaLat = 2.80100;    // PSI
    StaLocs[46].StaLon = 98.92400;
    StaLocs[46].StaElevation = 987.0;
    StaLocs[47].StaLat = 36.08670;    // LZH
    StaLocs[47].StaLon = 103.84400;
    StaLocs[47].StaElevation = 1560.0;
    StaLocs[48].StaLat = 40.84940;    // HHC
    StaLocs[48].StaLon = 111.56400;
    StaLocs[48].StaElevation = 1169.0;
    StaLocs[49].StaLat = 43.80140;    // CN2
    StaLocs[49].StaLon = 125.44800;
    StaLocs[49].StaElevation = 230.0;
    StaLocs[50].StaLat = 44.19981;    // USRK
    StaLocs[50].StaLon = 131.98881;
    StaLocs[50].StaElevation = 170.0;
    StaLocs[51].StaLat = 44.61640;    // MDJ
    StaLocs[51].StaLon = 129.59200;
    StaLocs[51].StaElevation = 250.0;
    StaLocs[52].StaLat = -13.95740;    // COEN
    StaLocs[52].StaLon = 143.17490;
    StaLocs[52].StaElevation = 285.4;
    StaLocs[53].StaLat = 25.56670;    // SHL
    StaLocs[53].StaLon = 91.88330;
    StaLocs[53].StaElevation = 1600.0;
    StaLocs[54].StaLat = 39.41060;    // GTA
    StaLocs[54].StaLon = 99.81440;
    StaLocs[54].StaElevation = 1341.0;
    StaLocs[55].StaLat = -19.93330;    // WRAB
    StaLocs[55].StaLon = 134.35000;
    StaLocs[55].StaElevation = 366.0;
    StaLocs[56].StaLat = -19.94260;    // WRA
    StaLocs[56].StaLon = 134.33900;
    StaLocs[56].StaElevation = 419.0;
    StaLocs[57].StaLat = 48.47300;    // HABR
    StaLocs[57].StaLon = 135.05100;
    StaLocs[57].StaElevation = 81.0;
    StaLocs[58].StaLat = 29.70000;    // LSA
    StaLocs[58].StaLon = 91.15000;
    StaLocs[58].StaElevation = 3789.0;
    StaLocs[59].StaLat = 47.86519;    // ULN
    StaLocs[59].StaLon = 107.05281;
    StaLocs[59].StaElevation = 1615.0;
    StaLocs[60].StaLat = 47.83469;    // SONM
    StaLocs[60].StaLon = 106.39500;
    StaLocs[60].StaElevation = 1415.8;
    StaLocs[61].StaLat = 26.86000;    // ODAN
    StaLocs[61].StaLon = 87.39000;
    StaLocs[61].StaElevation = 2045.0;
    StaLocs[62].StaLat = -23.66511;    // AS31
    StaLocs[62].StaLon = 133.90531;
    StaLocs[62].StaElevation = 627.3;
    StaLocs[63].StaLat = -23.66640;    // ASAR
    StaLocs[63].StaLon = 133.90400;
    StaLocs[63].StaElevation = 607.0;
    StaLocs[64].StaLat = 26.95000;    // RAMN
    StaLocs[64].StaLon = 86.60000;
    StaLocs[64].StaElevation = 2134.0;
    StaLocs[65].StaLat = 27.66000;    // JIRN
    StaLocs[65].StaLon = 86.19000;
    StaLocs[65].StaElevation = 3064.0;
    StaLocs[66].StaLat = -20.08830;    // CTA
    StaLocs[66].StaLon = 146.25400;
    StaLocs[66].StaElevation = 357.0;
    StaLocs[67].StaLat = 27.91060;    // GUN
    StaLocs[67].StaLon = 85.87940;
    StaLocs[67].StaElevation = 2900.0;
    StaLocs[68].StaLat = 27.57100;    // PKI
    StaLocs[68].StaLon = 85.40900;
    StaLocs[68].StaElevation = 2758.0;
    StaLocs[69].StaLat = 27.57700;    // PKIN
    StaLocs[69].StaLon = 85.39650;
    StaLocs[69].StaElevation = 2300.0;
    StaLocs[70].StaLat = 27.79000;    // KKN
    StaLocs[70].StaLon = 85.28000;
    StaLocs[70].StaElevation = 1920.0;
    StaLocs[71].StaLat = 27.60900;    // DMN
    StaLocs[71].StaLon = 85.10600;
    StaLocs[71].StaElevation = 2225.0;
    StaLocs[72].StaLat = 28.00300;    // GKN
    StaLocs[72].StaLon = 84.63700;
    StaLocs[72].StaElevation = 1478.0;
    StaLocs[73].StaLat = 50.38194;    // ZAK
    StaLocs[73].StaLon = 103.28056;
    StaLocs[73].StaElevation = 1200.0;
    StaLocs[74].StaLat = -9.43164;    // HNR
    StaLocs[74].StaLon = 159.94700;
    StaLocs[74].StaElevation = 72.0;
    StaLocs[75].StaLat = 28.35000;    // DANN
    StaLocs[75].StaLon = 83.76000;
    StaLocs[75].StaElevation = 2500.0;
    StaLocs[76].StaLat = 27.77000;    // KOLN
    StaLocs[76].StaLon = 83.60000;
    StaLocs[76].StaElevation = 1830.0;
    StaLocs[77].StaLat = 51.68069;    // TLY
    StaLocs[77].StaLon = 103.64381;
    StaLocs[77].StaElevation = 579.0;
    StaLocs[78].StaLat = 43.82111;    // WMQ
    StaLocs[78].StaLon = 87.69500;
    StaLocs[78].StaElevation = 897.0;
    StaLocs[79].StaLat = -30.77900;    // FORT
    StaLocs[79].StaLon = 128.05900;
    StaLocs[79].StaElevation = 165.0;
    StaLocs[80].StaLat = 51.13580;    // HVS
    StaLocs[80].StaLon = 93.70220;
    StaLocs[80].StaElevation = 1075.0;
    StaLocs[81].StaLat = -25.36911;    // EIDS
    StaLocs[81].StaLon = 151.08169;
    StaLocs[81].StaElevation = 216.0;
    StaLocs[82].StaLat = 53.10819;    // PETK
    StaLocs[82].StaLon = 157.69889;
    StaLocs[82].StaElevation = 400.0;
    StaLocs[83].StaLat = -32.92700;    // NWAO
    StaLocs[83].StaLon = 117.23400;
    StaLocs[83].StaElevation = 265.0;
    StaLocs[84].StaLat = -32.81020;    // BBOO
    StaLocs[84].StaLon = 136.05880;
    StaLocs[84].StaElevation = 321.0;
    StaLocs[85].StaLat = -31.87689;    // STKA
    StaLocs[85].StaLon = 141.59519;
    StaLocs[85].StaElevation = 272.3;
    StaLocs[86].StaLat = 46.79370;    // MK31
    StaLocs[86].StaLon = 82.29040;
    StaLocs[86].StaElevation = 609.4;
    StaLocs[87].StaLat = 46.79369;    // MKAR
    StaLocs[87].StaLon = 82.29039;
    StaLocs[87].StaElevation = 615.4;
    StaLocs[88].StaLat = 39.51667;    // KSH
    StaLocs[88].StaLon = 75.97306;
    StaLocs[88].StaElevation = 1314.0;
    StaLocs[89].StaLat = 42.24561;    // ULHL
    StaLocs[89].StaLon = 76.24169;
    StaLocs[89].StaElevation = 2040.0;
    StaLocs[90].StaLat = 53.94811;    // ZAA0
    StaLocs[90].StaLon = 84.81881;
    StaLocs[90].StaElevation = 229.4;
    StaLocs[91].StaLat = 53.94811;    // ZALV
    StaLocs[91].StaLon = 84.81881;
    StaLocs[91].StaElevation = 229.4;
    StaLocs[92].StaLat = 42.07781;    // KZA
    StaLocs[92].StaLon = 75.24961;
    StaLocs[92].StaElevation = 3520.0;
    StaLocs[93].StaLat = 42.92081;    // TKM2
    StaLocs[93].StaLon = 75.59661;
    StaLocs[93].StaElevation = 2020.0;
    StaLocs[94].StaLat = 42.22750;    // UCH
    StaLocs[94].StaLon = 74.51339;
    StaLocs[94].StaElevation = 3850.0;
    StaLocs[95].StaLat = 42.99861;    // CHMS
    StaLocs[95].StaLon = 74.75131;
    StaLocs[95].StaElevation = 655.0;
    StaLocs[96].StaLat = 43.26689;    // USP
    StaLocs[96].StaLon = 74.49969;
    StaLocs[96].StaElevation = 740.0;
    StaLocs[97].StaLat = 62.93278;    // SEY
    StaLocs[97].StaLon = 152.38222;
    StaLocs[97].StaElevation = 218.0;
    StaLocs[98].StaLat = 42.13111;    // AML
    StaLocs[98].StaLon = 73.69411;
    StaLocs[98].StaElevation = 3400.0;
    StaLocs[99].StaLat = 42.66150;    // EKS2
    StaLocs[99].StaLon = 73.77719;
    StaLocs[99].StaElevation = 1360.0;
    StaLocs[100].StaLat = 43.10340;    // KK31
    StaLocs[100].StaLon = 70.51150;
    StaLocs[100].StaElevation = 520.9;
    StaLocs[101].StaLat = 43.10340;    // KKAR
    StaLocs[101].StaLon = 70.51150;
    StaLocs[101].StaElevation = 520.9;
    StaLocs[102].StaLat = 71.64900;    // TIXI
    StaLocs[102].StaLon = 128.86650;
    StaLocs[102].StaElevation = 50.0;
    StaLocs[103].StaLat = 53.02490;    // BVA0
    StaLocs[103].StaLon = 70.38850;
    StaLocs[103].StaElevation = 420.0;
    StaLocs[104].StaLat = 53.05810;    // BRVK
    StaLocs[104].StaLon = 70.28280;
    StaLocs[104].StaElevation = 315.0;
    StaLocs[105].StaLat = 49.25560;    // ABKAR
    StaLocs[105].StaLon = 59.94310;
    StaLocs[105].StaElevation = 243.0;
    StaLocs[106].StaLat = 56.82700;    // SVE
    StaLocs[106].StaLon = 60.63700;
    StaLocs[106].StaElevation = 275.0;
    StaLocs[107].StaLat = 50.43480;    // AKTO
    StaLocs[107].StaLon = 58.01640;
    StaLocs[107].StaElevation = 379.0;
    StaLocs[108].StaLat = 56.42930;    // ARU
    StaLocs[108].StaLon = 58.56150;
    StaLocs[108].StaElevation = 260.0;
    StaLocs[109].StaLat = 57.22250;    // OHAK
    StaLocs[109].StaLon = -153.28750;
    StaLocs[109].StaElevation = 77.5;
    StaLocs[110].StaLat = 42.78810;    // ZEI
    StaLocs[110].StaLon = 43.90140;
    StaLocs[110].StaElevation = 1926.0;
    StaLocs[111].StaLat = 67.22739;    // COLD
    StaLocs[111].StaLon = -150.20131;
    StaLocs[111].StaElevation = 377.3;
    StaLocs[112].StaLat = 43.72861;    // KBZ
    StaLocs[112].StaLon = 42.89750;
    StaLocs[112].StaElevation = 0;
    StaLocs[113].StaLat = 43.95531;    // KIV
    StaLocs[113].StaLon = 42.68631;
    StaLocs[113].StaElevation = 1054.0;
    StaLocs[114].StaLat = 63.73233;    // MCK
    StaLocs[114].StaLon = -148.93489;
    StaLocs[114].StaElevation = 618.0;
    StaLocs[115].StaLat = 64.87381;    // COLA
    StaLocs[115].StaLon = -147.85111;
    StaLocs[115].StaElevation = 74.0;
    StaLocs[116].StaLat = 64.77139;    // ILAR
    StaLocs[116].StaLon = -146.88661;
    StaLocs[116].StaElevation = 419.0;
    StaLocs[117].StaLat = 63.64870;    // DOT
    StaLocs[117].StaLon = -144.06200;
    StaLocs[117].StaElevation = 671.0;
    StaLocs[118].StaLat = 55.11380;    // OBN
    StaLocs[118].StaLon = 36.56870;
    StaLocs[118].StaElevation = 0;
    StaLocs[119].StaLat = 64.77739;    // EGAK
    StaLocs[119].StaLon = -141.15811;
    StaLocs[119].StaElevation = 296.6;
    StaLocs[120].StaLat = 62.91820;    // JOF
    StaLocs[120].StaLon = 31.31240;
    StaLocs[120].StaElevation = 180.0;
    StaLocs[121].StaLat = 69.53490;    // ARCES
    StaLocs[121].StaLon = 25.50580;
    StaLocs[121].StaElevation = 403.0;
    StaLocs[122].StaLat = 68.30650;    // INK
    StaLocs[122].StaLon = -133.52540;
    StaLocs[122].StaElevation = 44.0;
    StaLocs[123].StaLat = 62.11280;    // KAF
    StaLocs[123].StaLon = 26.30610;
    StaLocs[123].StaElevation = 205.0;
    StaLocs[124].StaLat = 61.44360;    // FINES
    StaLocs[124].StaLon = 26.07710;
    StaLocs[124].StaElevation = 150.0;
    StaLocs[125].StaLat = 39.72500;    // BRTR
    StaLocs[125].StaLon = 33.63900;
    StaLocs[125].StaElevation = 1440.0;
    StaLocs[126].StaLat = 50.70120;    // AKASG
    StaLocs[126].StaLon = 29.22420;
    StaLocs[126].StaElevation = 160.0;
    StaLocs[127].StaLat = 47.64400;    // BUR08
    StaLocs[127].StaLon = 25.20019;
    StaLocs[127].StaElevation = 1214.0;
    StaLocs[128].StaLat = 45.49089;    // MLR
    StaLocs[128].StaLon = 25.94500;
    StaLocs[128].StaElevation = 1360.0;
    StaLocs[129].StaLat = 74.68670;    // RES
    StaLocs[129].StaLon = -94.90000;
    StaLocs[129].StaElevation = 15.0;
    StaLocs[130].StaLat = 48.93330;    // KOLS
    StaLocs[130].StaLon = 22.27310;
    StaLocs[130].StaElevation = 460.0;
    StaLocs[131].StaLat = 48.90219;    // CRVS
    StaLocs[131].StaLon = 21.46139;
    StaLocs[131].StaElevation = 476.0;
    StaLocs[132].StaLat = 61.03972;    // NOA
    StaLocs[132].StaLon = 11.21475;
    StaLocs[132].StaElevation = 717.0;
    StaLocs[133].StaLat = 62.49322;    // YKA
    StaLocs[133].StaLon = -114.60528;
    StaLocs[133].StaElevation = 197.0;
    StaLocs[134].StaLat = 41.30158;    // PHWY
    StaLocs[134].StaLon = -105.45775;
    StaLocs[134].StaElevation = 2645.0;
    StaLocs[135].StaLat = 29.33380;    // TXAR
    StaLocs[135].StaLon = -103.66700;
    StaLocs[135].StaElevation = 1013.0;
    StaLocs[136].StaLat = 13.14769;    // TORD
    StaLocs[136].StaLon = 1.69469;
    StaLocs[136].StaElevation = 214.3;
    StaLocs[137].StaLat = -40.73278;    // PLCA
    StaLocs[137].StaLon = -70.55083;
    StaLocs[137].StaElevation = 956.0;
    StaLocs[138].StaLat = 8.88610;    // SDV
    StaLocs[138].StaLon = -70.63330;
    StaLocs[138].StaElevation = 1580.0;
    StaLocs[139].StaLat = -16.28790;    // LPAZ
    StaLocs[139].StaLon = -68.13070;
    StaLocs[139].StaElevation = 4774.0;
/*
 *  set Assocs structure
 */
    Assocs[0].StaInd = 0;    // PVCP
    Assocs[0].arid = 1;
    strcpy(Assocs[0].PhaseHint, "Pn");
    Assocs[0].phaseFixed = 0;
    Assocs[0].ArrivalTime = 1262324255.04;
    Assocs[0].BackAzimuth = ILOC_NULLVAL;
    Assocs[0].Slowness = ILOC_NULLVAL;
    Assocs[0].Timedef = 1;
    Assocs[0].Azimdef = 0;
    Assocs[0].Slowdef = 0;
    Assocs[0].Deltim = ILOC_NULLVAL;
    Assocs[1].StaInd = 0;    // PVCP
    Assocs[1].arid = 2;
    strcpy(Assocs[1].PhaseHint, "Sn");
    Assocs[1].phaseFixed = 0;
    Assocs[1].ArrivalTime = 1262324272.47;
    Assocs[1].BackAzimuth = ILOC_NULLVAL;
    Assocs[1].Slowness = ILOC_NULLVAL;
    Assocs[1].Timedef = 1;
    Assocs[1].Azimdef = 0;
    Assocs[1].Slowdef = 0;
    Assocs[1].Deltim = ILOC_NULLVAL;
    Assocs[2].StaInd = 1;    // CNP
    Assocs[2].arid = 3;
    strcpy(Assocs[2].PhaseHint, "Pn");
    Assocs[2].phaseFixed = 0;
    Assocs[2].ArrivalTime = 1262324258.18;
    Assocs[2].BackAzimuth = ILOC_NULLVAL;
    Assocs[2].Slowness = ILOC_NULLVAL;
    Assocs[2].Timedef = 1;
    Assocs[2].Azimdef = 0;
    Assocs[2].Slowdef = 0;
    Assocs[2].Deltim = ILOC_NULLVAL;
    Assocs[3].StaInd = 2;    // PLP
    Assocs[3].arid = 4;
    strcpy(Assocs[3].PhaseHint, "Pn");
    Assocs[3].phaseFixed = 0;
    Assocs[3].ArrivalTime = 1262324273.35;
    Assocs[3].BackAzimuth = ILOC_NULLVAL;
    Assocs[3].Slowness = ILOC_NULLVAL;
    Assocs[3].Timedef = 1;
    Assocs[3].Azimdef = 0;
    Assocs[3].Slowdef = 0;
    Assocs[3].Deltim = ILOC_NULLVAL;
    Assocs[4].StaInd = 2;    // PLP
    Assocs[4].arid = 5;
    strcpy(Assocs[4].PhaseHint, "Sn");
    Assocs[4].phaseFixed = 0;
    Assocs[4].ArrivalTime = 1262324304.56;
    Assocs[4].BackAzimuth = ILOC_NULLVAL;
    Assocs[4].Slowness = ILOC_NULLVAL;
    Assocs[4].Timedef = 1;
    Assocs[4].Azimdef = 0;
    Assocs[4].Slowdef = 0;
    Assocs[4].Deltim = ILOC_NULLVAL;
    Assocs[5].StaInd = 3;    // OCLP
    Assocs[5].arid = 6;
    strcpy(Assocs[5].PhaseHint, "Pn");
    Assocs[5].phaseFixed = 0;
    Assocs[5].ArrivalTime = 1262324278.33;
    Assocs[5].BackAzimuth = ILOC_NULLVAL;
    Assocs[5].Slowness = ILOC_NULLVAL;
    Assocs[5].Timedef = 1;
    Assocs[5].Azimdef = 0;
    Assocs[5].Slowdef = 0;
    Assocs[5].Deltim = ILOC_NULLVAL;
    Assocs[6].StaInd = 4;    // AUQP
    Assocs[6].arid = 7;
    strcpy(Assocs[6].PhaseHint, "Pn");
    Assocs[6].phaseFixed = 0;
    Assocs[6].ArrivalTime = 1262324278;
    Assocs[6].BackAzimuth = ILOC_NULLVAL;
    Assocs[6].Slowness = ILOC_NULLVAL;
    Assocs[6].Timedef = 1;
    Assocs[6].Azimdef = 0;
    Assocs[6].Slowdef = 0;
    Assocs[6].Deltim = ILOC_NULLVAL;
    Assocs[7].StaInd = 5;    // GQP
    Assocs[7].arid = 8;
    strcpy(Assocs[7].PhaseHint, "Pn");
    Assocs[7].phaseFixed = 0;
    Assocs[7].ArrivalTime = 1262324278.84;
    Assocs[7].BackAzimuth = ILOC_NULLVAL;
    Assocs[7].Slowness = ILOC_NULLVAL;
    Assocs[7].Timedef = 1;
    Assocs[7].Azimdef = 0;
    Assocs[7].Slowdef = 0;
    Assocs[7].Deltim = ILOC_NULLVAL;
    Assocs[8].StaInd = 5;    // GQP
    Assocs[8].arid = 9;
    strcpy(Assocs[8].PhaseHint, "Sn");
    Assocs[8].phaseFixed = 0;
    Assocs[8].ArrivalTime = 1262324313.11;
    Assocs[8].BackAzimuth = ILOC_NULLVAL;
    Assocs[8].Slowness = ILOC_NULLVAL;
    Assocs[8].Timedef = 1;
    Assocs[8].Azimdef = 0;
    Assocs[8].Slowdef = 0;
    Assocs[8].Deltim = ILOC_NULLVAL;
    Assocs[9].StaInd = 6;    // RCP
    Assocs[9].arid = 10;
    strcpy(Assocs[9].PhaseHint, "Pn");
    Assocs[9].phaseFixed = 0;
    Assocs[9].ArrivalTime = 1262324286.72;
    Assocs[9].BackAzimuth = ILOC_NULLVAL;
    Assocs[9].Slowness = ILOC_NULLVAL;
    Assocs[9].Timedef = 1;
    Assocs[9].Azimdef = 0;
    Assocs[9].Slowdef = 0;
    Assocs[9].Deltim = ILOC_NULLVAL;
    Assocs[10].StaInd = 6;    // RCP
    Assocs[10].arid = 11;
    strcpy(Assocs[10].PhaseHint, "Sn");
    Assocs[10].phaseFixed = 0;
    Assocs[10].ArrivalTime = 1262324326.72;
    Assocs[10].BackAzimuth = ILOC_NULLVAL;
    Assocs[10].Slowness = ILOC_NULLVAL;
    Assocs[10].Timedef = 1;
    Assocs[10].Azimdef = 0;
    Assocs[10].Slowdef = 0;
    Assocs[10].Deltim = ILOC_NULLVAL;
    Assocs[11].StaInd = 7;    // MSLP
    Assocs[11].arid = 12;
    strcpy(Assocs[11].PhaseHint, "Pn");
    Assocs[11].phaseFixed = 0;
    Assocs[11].ArrivalTime = 1262324289.25;
    Assocs[11].BackAzimuth = ILOC_NULLVAL;
    Assocs[11].Slowness = ILOC_NULLVAL;
    Assocs[11].Timedef = 1;
    Assocs[11].Azimdef = 0;
    Assocs[11].Slowdef = 0;
    Assocs[11].Deltim = ILOC_NULLVAL;
    Assocs[12].StaInd = 8;    // BOAC
    Assocs[12].arid = 13;
    strcpy(Assocs[12].PhaseHint, "Pn");
    Assocs[12].phaseFixed = 0;
    Assocs[12].ArrivalTime = 1262324287.69;
    Assocs[12].BackAzimuth = ILOC_NULLVAL;
    Assocs[12].Slowness = ILOC_NULLVAL;
    Assocs[12].Timedef = 1;
    Assocs[12].Azimdef = 0;
    Assocs[12].Slowdef = 0;
    Assocs[12].Deltim = ILOC_NULLVAL;
    Assocs[13].StaInd = 9;    // LLP
    Assocs[13].arid = 14;
    strcpy(Assocs[13].PhaseHint, "Pn");
    Assocs[13].phaseFixed = 0;
    Assocs[13].ArrivalTime = 1262324291.76;
    Assocs[13].BackAzimuth = ILOC_NULLVAL;
    Assocs[13].Slowness = ILOC_NULLVAL;
    Assocs[13].Timedef = 1;
    Assocs[13].Azimdef = 0;
    Assocs[13].Slowdef = 0;
    Assocs[13].Deltim = ILOC_NULLVAL;
    Assocs[14].StaInd = 9;    // LLP
    Assocs[14].arid = 15;
    strcpy(Assocs[14].PhaseHint, "Sn");
    Assocs[14].phaseFixed = 0;
    Assocs[14].ArrivalTime = 1262324327.14;
    Assocs[14].BackAzimuth = ILOC_NULLVAL;
    Assocs[14].Slowness = ILOC_NULLVAL;
    Assocs[14].Timedef = 1;
    Assocs[14].Azimdef = 0;
    Assocs[14].Slowdef = 0;
    Assocs[14].Deltim = ILOC_NULLVAL;
    Assocs[15].StaInd = 10;    // GUIM
    Assocs[15].arid = 16;
    strcpy(Assocs[15].PhaseHint, "Pn");
    Assocs[15].phaseFixed = 0;
    Assocs[15].ArrivalTime = 1262324297.19;
    Assocs[15].BackAzimuth = ILOC_NULLVAL;
    Assocs[15].Slowness = ILOC_NULLVAL;
    Assocs[15].Timedef = 1;
    Assocs[15].Azimdef = 0;
    Assocs[15].Slowdef = 0;
    Assocs[15].Deltim = ILOC_NULLVAL;
    Assocs[16].StaInd = 11;    // PALP
    Assocs[16].arid = 17;
    strcpy(Assocs[16].PhaseHint, "Pn");
    Assocs[16].phaseFixed = 0;
    Assocs[16].ArrivalTime = 1262324295.45;
    Assocs[16].BackAzimuth = ILOC_NULLVAL;
    Assocs[16].Slowness = ILOC_NULLVAL;
    Assocs[16].Timedef = 1;
    Assocs[16].Azimdef = 0;
    Assocs[16].Slowdef = 0;
    Assocs[16].Deltim = ILOC_NULLVAL;
    Assocs[17].StaInd = 11;    // PALP
    Assocs[17].arid = 18;
    strcpy(Assocs[17].PhaseHint, "Pn");
    Assocs[17].phaseFixed = 0;
    Assocs[17].ArrivalTime = 1262324298.25;
    Assocs[17].BackAzimuth = ILOC_NULLVAL;
    Assocs[17].Slowness = ILOC_NULLVAL;
    Assocs[17].Timedef = 1;
    Assocs[17].Azimdef = 0;
    Assocs[17].Slowdef = 0;
    Assocs[17].Deltim = ILOC_NULLVAL;
    Assocs[18].StaInd = 12;    // TGY
    Assocs[18].arid = 19;
    strcpy(Assocs[18].PhaseHint, "Pn");
    Assocs[18].phaseFixed = 0;
    Assocs[18].ArrivalTime = 1262324299.94;
    Assocs[18].BackAzimuth = ILOC_NULLVAL;
    Assocs[18].Slowness = ILOC_NULLVAL;
    Assocs[18].Timedef = 1;
    Assocs[18].Azimdef = 0;
    Assocs[18].Slowdef = 0;
    Assocs[18].Deltim = ILOC_NULLVAL;
    Assocs[19].StaInd = 12;    // TGY
    Assocs[19].arid = 20;
    strcpy(Assocs[19].PhaseHint, "Sn");
    Assocs[19].phaseFixed = 0;
    Assocs[19].ArrivalTime = 1262324352.79;
    Assocs[19].BackAzimuth = ILOC_NULLVAL;
    Assocs[19].Slowness = ILOC_NULLVAL;
    Assocs[19].Timedef = 1;
    Assocs[19].Azimdef = 0;
    Assocs[19].Slowdef = 0;
    Assocs[19].Deltim = ILOC_NULLVAL;
    Assocs[20].StaInd = 12;    // TGY
    Assocs[20].arid = 21;
    strcpy(Assocs[20].PhaseHint, "Pn");
    Assocs[20].phaseFixed = 0;
    Assocs[20].ArrivalTime = 1262324299.937;
    Assocs[20].BackAzimuth = 77.8;
    Assocs[20].Slowness = 0.30;
    Assocs[20].Timedef = 1;
    Assocs[20].Azimdef = 1;
    Assocs[20].Slowdef = 1;
    Assocs[20].Deltim = ILOC_NULLVAL;
    Assocs[21].StaInd = 12;    // TGY
    Assocs[21].arid = 22;
    strcpy(Assocs[21].PhaseHint, "Sn");
    Assocs[21].phaseFixed = 0;
    Assocs[21].ArrivalTime = 1262324352.787;
    Assocs[21].BackAzimuth = 42.6;
    Assocs[21].Slowness = 21.50;
    Assocs[21].Timedef = 1;
    Assocs[21].Azimdef = 1;
    Assocs[21].Slowdef = 1;
    Assocs[21].Deltim = ILOC_NULLVAL;
    Assocs[22].StaInd = 13;    // PCPH
    Assocs[22].arid = 23;
    strcpy(Assocs[22].PhaseHint, "Pn");
    Assocs[22].phaseFixed = 0;
    Assocs[22].ArrivalTime = 1262324299.43;
    Assocs[22].BackAzimuth = ILOC_NULLVAL;
    Assocs[22].Slowness = ILOC_NULLVAL;
    Assocs[22].Timedef = 1;
    Assocs[22].Azimdef = 0;
    Assocs[22].Slowdef = 0;
    Assocs[22].Deltim = ILOC_NULLVAL;
    Assocs[23].StaInd = 13;    // PCPH
    Assocs[23].arid = 24;
    strcpy(Assocs[23].PhaseHint, "Sn");
    Assocs[23].phaseFixed = 0;
    Assocs[23].ArrivalTime = 1262324348.9;
    Assocs[23].BackAzimuth = ILOC_NULLVAL;
    Assocs[23].Slowness = ILOC_NULLVAL;
    Assocs[23].Timedef = 1;
    Assocs[23].Azimdef = 0;
    Assocs[23].Slowdef = 0;
    Assocs[23].Deltim = ILOC_NULLVAL;
    Assocs[24].StaInd = 14;    // CVP
    Assocs[24].arid = 25;
    strcpy(Assocs[24].PhaseHint, "Pn");
    Assocs[24].phaseFixed = 0;
    Assocs[24].ArrivalTime = 1262324310.7;
    Assocs[24].BackAzimuth = ILOC_NULLVAL;
    Assocs[24].Slowness = ILOC_NULLVAL;
    Assocs[24].Timedef = 1;
    Assocs[24].Azimdef = 0;
    Assocs[24].Slowdef = 0;
    Assocs[24].Deltim = ILOC_NULLVAL;
    Assocs[25].StaInd = 15;    // BUSP
    Assocs[25].arid = 26;
    strcpy(Assocs[25].PhaseHint, "Pn");
    Assocs[25].phaseFixed = 0;
    Assocs[25].ArrivalTime = 1262324314.35;
    Assocs[25].BackAzimuth = ILOC_NULLVAL;
    Assocs[25].Slowness = ILOC_NULLVAL;
    Assocs[25].Timedef = 1;
    Assocs[25].Azimdef = 0;
    Assocs[25].Slowdef = 0;
    Assocs[25].Deltim = ILOC_NULLVAL;
    Assocs[26].StaInd = 16;    // SZP
    Assocs[26].arid = 27;
    strcpy(Assocs[26].PhaseHint, "Pn");
    Assocs[26].phaseFixed = 0;
    Assocs[26].ArrivalTime = 1262324324.76;
    Assocs[26].BackAzimuth = ILOC_NULLVAL;
    Assocs[26].Slowness = ILOC_NULLVAL;
    Assocs[26].Timedef = 1;
    Assocs[26].Azimdef = 0;
    Assocs[26].Slowdef = 0;
    Assocs[26].Deltim = ILOC_NULLVAL;
    Assocs[27].StaInd = 17;    // DAV
    Assocs[27].arid = 28;
    strcpy(Assocs[27].PhaseHint, "Pn");
    Assocs[27].phaseFixed = 0;
    Assocs[27].ArrivalTime = 1262324327.54;
    Assocs[27].BackAzimuth = ILOC_NULLVAL;
    Assocs[27].Slowness = ILOC_NULLVAL;
    Assocs[27].Timedef = 1;
    Assocs[27].Azimdef = 0;
    Assocs[27].Slowdef = 0;
    Assocs[27].Deltim = ILOC_NULLVAL;
    Assocs[28].StaInd = 17;    // DAV
    Assocs[28].arid = 29;
    strcpy(Assocs[28].PhaseHint, "Sn");
    Assocs[28].phaseFixed = 0;
    Assocs[28].ArrivalTime = 1262324400.06;
    Assocs[28].BackAzimuth = ILOC_NULLVAL;
    Assocs[28].Slowness = ILOC_NULLVAL;
    Assocs[28].Timedef = 1;
    Assocs[28].Azimdef = 0;
    Assocs[28].Slowdef = 0;
    Assocs[28].Deltim = ILOC_NULLVAL;
    Assocs[29].StaInd = 18;    // YULB
    Assocs[29].arid = 30;
    strcpy(Assocs[29].PhaseHint, "Pn");
    Assocs[29].phaseFixed = 0;
    Assocs[29].ArrivalTime = 1262324382.62;
    Assocs[29].BackAzimuth = ILOC_NULLVAL;
    Assocs[29].Slowness = ILOC_NULLVAL;
    Assocs[29].Timedef = 1;
    Assocs[29].Azimdef = 0;
    Assocs[29].Slowdef = 0;
    Assocs[29].Deltim = ILOC_NULLVAL;
    Assocs[30].StaInd = 18;    // YULB
    Assocs[30].arid = 31;
    strcpy(Assocs[30].PhaseHint, "Sn");
    Assocs[30].phaseFixed = 0;
    Assocs[30].ArrivalTime = 1262324488.62;
    Assocs[30].BackAzimuth = ILOC_NULLVAL;
    Assocs[30].Slowness = ILOC_NULLVAL;
    Assocs[30].Timedef = 1;
    Assocs[30].Azimdef = 0;
    Assocs[30].Slowdef = 0;
    Assocs[30].Deltim = ILOC_NULLVAL;
    Assocs[31].StaInd = 19;    // JOW
    Assocs[31].arid = 32;
    strcpy(Assocs[31].PhaseHint, "Pn");
    Assocs[31].phaseFixed = 0;
    Assocs[31].ArrivalTime = 1262324420.75;
    Assocs[31].BackAzimuth = ILOC_NULLVAL;
    Assocs[31].Slowness = ILOC_NULLVAL;
    Assocs[31].Timedef = 1;
    Assocs[31].Azimdef = 0;
    Assocs[31].Slowdef = 0;
    Assocs[31].Deltim = ILOC_NULLVAL;
    Assocs[32].StaInd = 19;    // JOW
    Assocs[32].arid = 33;
    strcpy(Assocs[32].PhaseHint, "Pn");
    Assocs[32].phaseFixed = 0;
    Assocs[32].ArrivalTime = 1262324420.5;
    Assocs[32].BackAzimuth = 167.9;
    Assocs[32].Slowness = 18.00;
    Assocs[32].Timedef = 1;
    Assocs[32].Azimdef = 1;
    Assocs[32].Slowdef = 1;
    Assocs[32].Deltim = ILOC_NULLVAL;
    Assocs[33].StaInd = 19;    // JOW
    Assocs[33].arid = 34;
    strcpy(Assocs[33].PhaseHint, "LR");
    Assocs[33].phaseFixed = 0;
    Assocs[33].ArrivalTime = 1262324657.873;
    Assocs[33].BackAzimuth = 0.7;
    Assocs[33].Slowness = 32.10;
    Assocs[33].Timedef = 1;
    Assocs[33].Azimdef = 1;
    Assocs[33].Slowdef = 1;
    Assocs[33].Deltim = ILOC_NULLVAL;
    Assocs[34].StaInd = 20;    // MRSI
    Assocs[34].arid = 35;
    strcpy(Assocs[34].PhaseHint, "P");
    Assocs[34].phaseFixed = 0;
    Assocs[34].ArrivalTime = 1262324434.1;
    Assocs[34].BackAzimuth = ILOC_NULLVAL;
    Assocs[34].Slowness = ILOC_NULLVAL;
    Assocs[34].Timedef = 1;
    Assocs[34].Azimdef = 0;
    Assocs[34].Slowdef = 0;
    Assocs[34].Deltim = ILOC_NULLVAL;
    Assocs[35].StaInd = 21;    // APSI
    Assocs[35].arid = 36;
    strcpy(Assocs[35].PhaseHint, "P");
    Assocs[35].phaseFixed = 0;
    Assocs[35].ArrivalTime = 1262324454.1;
    Assocs[35].BackAzimuth = ILOC_NULLVAL;
    Assocs[35].Slowness = ILOC_NULLVAL;
    Assocs[35].Timedef = 1;
    Assocs[35].Azimdef = 0;
    Assocs[35].Slowdef = 0;
    Assocs[35].Deltim = ILOC_NULLVAL;
    Assocs[36].StaInd = 22;    // SWI
    Assocs[36].arid = 37;
    strcpy(Assocs[36].PhaseHint, "P");
    Assocs[36].phaseFixed = 0;
    Assocs[36].ArrivalTime = 1262324453.4;
    Assocs[36].BackAzimuth = ILOC_NULLVAL;
    Assocs[36].Slowness = ILOC_NULLVAL;
    Assocs[36].Timedef = 1;
    Assocs[36].Azimdef = 0;
    Assocs[36].Slowdef = 0;
    Assocs[36].Deltim = ILOC_NULLVAL;
    Assocs[37].StaInd = 23;    // QIZ
    Assocs[37].arid = 38;
    strcpy(Assocs[37].PhaseHint, "Pn");
    Assocs[37].phaseFixed = 0;
    Assocs[37].ArrivalTime = 1262324456.2;
    Assocs[37].BackAzimuth = ILOC_NULLVAL;
    Assocs[37].Slowness = ILOC_NULLVAL;
    Assocs[37].Timedef = 1;
    Assocs[37].Azimdef = 0;
    Assocs[37].Slowdef = 0;
    Assocs[37].Deltim = ILOC_NULLVAL;
    Assocs[38].StaInd = 23;    // QIZ
    Assocs[38].arid = 39;
    strcpy(Assocs[38].PhaseHint, "Sn");
    Assocs[38].phaseFixed = 0;
    Assocs[38].ArrivalTime = 1262324635.6;
    Assocs[38].BackAzimuth = ILOC_NULLVAL;
    Assocs[38].Slowness = ILOC_NULLVAL;
    Assocs[38].Timedef = 1;
    Assocs[38].Azimdef = 0;
    Assocs[38].Slowdef = 0;
    Assocs[38].Deltim = ILOC_NULLVAL;
    Assocs[39].StaInd = 24;    // NLAI
    Assocs[39].arid = 40;
    strcpy(Assocs[39].PhaseHint, "P");
    Assocs[39].phaseFixed = 0;
    Assocs[39].ArrivalTime = 1262324470.5;
    Assocs[39].BackAzimuth = ILOC_NULLVAL;
    Assocs[39].Slowness = ILOC_NULLVAL;
    Assocs[39].Timedef = 1;
    Assocs[39].Azimdef = 0;
    Assocs[39].Slowdef = 0;
    Assocs[39].Deltim = ILOC_NULLVAL;
    Assocs[40].StaInd = 25;    // MSAI
    Assocs[40].arid = 41;
    strcpy(Assocs[40].PhaseHint, "P");
    Assocs[40].phaseFixed = 0;
    Assocs[40].ArrivalTime = 1262324480.2;
    Assocs[40].BackAzimuth = ILOC_NULLVAL;
    Assocs[40].Slowness = ILOC_NULLVAL;
    Assocs[40].Timedef = 1;
    Assocs[40].Azimdef = 0;
    Assocs[40].Slowdef = 0;
    Assocs[40].Deltim = ILOC_NULLVAL;
    Assocs[41].StaInd = 26;    // SBUM
    Assocs[41].arid = 42;
    strcpy(Assocs[41].PhaseHint, "P");
    Assocs[41].phaseFixed = 0;
    Assocs[41].ArrivalTime = 1262324475.7;
    Assocs[41].BackAzimuth = ILOC_NULLVAL;
    Assocs[41].Slowness = ILOC_NULLVAL;
    Assocs[41].Timedef = 1;
    Assocs[41].Azimdef = 0;
    Assocs[41].Slowdef = 0;
    Assocs[41].Deltim = ILOC_NULLVAL;
    Assocs[42].StaInd = 27;    // TTSI
    Assocs[42].arid = 43;
    strcpy(Assocs[42].PhaseHint, "P");
    Assocs[42].phaseFixed = 0;
    Assocs[42].ArrivalTime = 1262324480.4;
    Assocs[42].BackAzimuth = ILOC_NULLVAL;
    Assocs[42].Slowness = ILOC_NULLVAL;
    Assocs[42].Timedef = 1;
    Assocs[42].Azimdef = 0;
    Assocs[42].Slowdef = 0;
    Assocs[42].Deltim = ILOC_NULLVAL;
    Assocs[43].StaInd = 28;    // SPSI
    Assocs[43].arid = 44;
    strcpy(Assocs[43].PhaseHint, "Pn");
    Assocs[43].phaseFixed = 0;
    Assocs[43].ArrivalTime = 1262324491;
    Assocs[43].BackAzimuth = ILOC_NULLVAL;
    Assocs[43].Slowness = ILOC_NULLVAL;
    Assocs[43].Timedef = 1;
    Assocs[43].Azimdef = 0;
    Assocs[43].Slowdef = 0;
    Assocs[43].Deltim = ILOC_NULLVAL;
    Assocs[44].StaInd = 29;    // BNSI
    Assocs[44].arid = 45;
    strcpy(Assocs[44].PhaseHint, "Pn");
    Assocs[44].phaseFixed = 0;
    Assocs[44].ArrivalTime = 1262324495.5;
    Assocs[44].BackAzimuth = ILOC_NULLVAL;
    Assocs[44].Slowness = ILOC_NULLVAL;
    Assocs[44].Timedef = 1;
    Assocs[44].Azimdef = 0;
    Assocs[44].Slowdef = 0;
    Assocs[44].Deltim = ILOC_NULLVAL;
    Assocs[45].StaInd = 30;    // NJ2
    Assocs[45].arid = 46;
    strcpy(Assocs[45].PhaseHint, "Pn");
    Assocs[45].phaseFixed = 0;
    Assocs[45].ArrivalTime = 1262324497.5;
    Assocs[45].BackAzimuth = ILOC_NULLVAL;
    Assocs[45].Slowness = ILOC_NULLVAL;
    Assocs[45].Timedef = 1;
    Assocs[45].Azimdef = 0;
    Assocs[45].Slowdef = 0;
    Assocs[45].Deltim = ILOC_NULLVAL;
    Assocs[46].StaInd = 31;    // JNU
    Assocs[46].arid = 47;
    strcpy(Assocs[46].PhaseHint, "P");
    Assocs[46].phaseFixed = 0;
    Assocs[46].ArrivalTime = 1262324503.58;
    Assocs[46].BackAzimuth = ILOC_NULLVAL;
    Assocs[46].Slowness = ILOC_NULLVAL;
    Assocs[46].Timedef = 1;
    Assocs[46].Azimdef = 0;
    Assocs[46].Slowdef = 0;
    Assocs[46].Deltim = ILOC_NULLVAL;
    Assocs[47].StaInd = 31;    // JNU
    Assocs[47].arid = 48;
    strcpy(Assocs[47].PhaseHint, "P");
    Assocs[47].phaseFixed = 0;
    Assocs[47].ArrivalTime = 1262324502.275;
    Assocs[47].BackAzimuth = 11.3;
    Assocs[47].Slowness = 10.30;
    Assocs[47].Timedef = 1;
    Assocs[47].Azimdef = 1;
    Assocs[47].Slowdef = 1;
    Assocs[47].Deltim = ILOC_NULLVAL;
    Assocs[48].StaInd = 31;    // JNU
    Assocs[48].arid = 49;
    strcpy(Assocs[48].PhaseHint, "LR");
    Assocs[48].phaseFixed = 0;
    Assocs[48].ArrivalTime = 1262324938.651;
    Assocs[48].BackAzimuth = 41.9;
    Assocs[48].Slowness = 35.50;
    Assocs[48].Timedef = 1;
    Assocs[48].Azimdef = 1;
    Assocs[48].Slowdef = 1;
    Assocs[48].Deltim = ILOC_NULLVAL;
    Assocs[49].StaInd = 32;    // BBKI
    Assocs[49].arid = 50;
    strcpy(Assocs[49].PhaseHint, "P");
    Assocs[49].phaseFixed = 0;
    Assocs[49].ArrivalTime = 1262324505.7;
    Assocs[49].BackAzimuth = ILOC_NULLVAL;
    Assocs[49].Slowness = ILOC_NULLVAL;
    Assocs[49].Timedef = 1;
    Assocs[49].Azimdef = 0;
    Assocs[49].Slowdef = 0;
    Assocs[49].Deltim = ILOC_NULLVAL;
    Assocs[50].StaInd = 33;    // GYA
    Assocs[50].arid = 51;
    strcpy(Assocs[50].PhaseHint, "P");
    Assocs[50].phaseFixed = 0;
    Assocs[50].ArrivalTime = 1262324527.9;
    Assocs[50].BackAzimuth = ILOC_NULLVAL;
    Assocs[50].Slowness = ILOC_NULLVAL;
    Assocs[50].Timedef = 1;
    Assocs[50].Azimdef = 0;
    Assocs[50].Slowdef = 0;
    Assocs[50].Deltim = ILOC_NULLVAL;
    Assocs[51].StaInd = 33;    // GYA
    Assocs[51].arid = 52;
    strcpy(Assocs[51].PhaseHint, "PnPn");
    Assocs[51].phaseFixed = 0;
    Assocs[51].ArrivalTime = 1262324551.8;
    Assocs[51].BackAzimuth = ILOC_NULLVAL;
    Assocs[51].Slowness = ILOC_NULLVAL;
    Assocs[51].Timedef = 1;
    Assocs[51].Azimdef = 0;
    Assocs[51].Slowdef = 0;
    Assocs[51].Deltim = ILOC_NULLVAL;
    Assocs[52].StaInd = 33;    // GYA
    Assocs[52].arid = 53;
    strcpy(Assocs[52].PhaseHint, "S");
    Assocs[52].phaseFixed = 0;
    Assocs[52].ArrivalTime = 1262324762.4;
    Assocs[52].BackAzimuth = ILOC_NULLVAL;
    Assocs[52].Slowness = ILOC_NULLVAL;
    Assocs[52].Timedef = 1;
    Assocs[52].Azimdef = 0;
    Assocs[52].Slowdef = 0;
    Assocs[52].Deltim = ILOC_NULLVAL;
    Assocs[53].StaInd = 33;    // GYA
    Assocs[53].arid = 54;
    strcpy(Assocs[53].PhaseHint, "SnSn");
    Assocs[53].phaseFixed = 0;
    Assocs[53].ArrivalTime = 1262324799;
    Assocs[53].BackAzimuth = ILOC_NULLVAL;
    Assocs[53].Slowness = ILOC_NULLVAL;
    Assocs[53].Timedef = 1;
    Assocs[53].Azimdef = 0;
    Assocs[53].Slowdef = 0;
    Assocs[53].Deltim = ILOC_NULLVAL;
    Assocs[54].StaInd = 33;    // GYA
    Assocs[54].arid = 55;
    strcpy(Assocs[54].PhaseHint, "ScP");
    Assocs[54].phaseFixed = 0;
    Assocs[54].ArrivalTime = 1262324982.2;
    Assocs[54].BackAzimuth = ILOC_NULLVAL;
    Assocs[54].Slowness = ILOC_NULLVAL;
    Assocs[54].Timedef = 1;
    Assocs[54].Azimdef = 0;
    Assocs[54].Slowdef = 0;
    Assocs[54].Deltim = ILOC_NULLVAL;
    Assocs[55].StaInd = 34;    // SLVN
    Assocs[55].arid = 56;
    strcpy(Assocs[55].PhaseHint, "P");
    Assocs[55].phaseFixed = 0;
    Assocs[55].ArrivalTime = 1262324526.37;
    Assocs[55].BackAzimuth = ILOC_NULLVAL;
    Assocs[55].Slowness = ILOC_NULLVAL;
    Assocs[55].Timedef = 1;
    Assocs[55].Azimdef = 0;
    Assocs[55].Slowdef = 0;
    Assocs[55].Deltim = ILOC_NULLVAL;
    Assocs[56].StaInd = 35;    // ENH
    Assocs[56].arid = 57;
    strcpy(Assocs[56].PhaseHint, "P");
    Assocs[56].phaseFixed = 0;
    Assocs[56].ArrivalTime = 1262324527.7;
    Assocs[56].BackAzimuth = ILOC_NULLVAL;
    Assocs[56].Slowness = ILOC_NULLVAL;
    Assocs[56].Timedef = 1;
    Assocs[56].Azimdef = 0;
    Assocs[56].Slowdef = 0;
    Assocs[56].Deltim = ILOC_NULLVAL;
    Assocs[57].StaInd = 36;    // KSAR
    Assocs[57].arid = 58;
    strcpy(Assocs[57].PhaseHint, "P");
    Assocs[57].phaseFixed = 0;
    Assocs[57].ArrivalTime = 1262324544.8;
    Assocs[57].BackAzimuth = ILOC_NULLVAL;
    Assocs[57].Slowness = ILOC_NULLVAL;
    Assocs[57].Timedef = 1;
    Assocs[57].Azimdef = 0;
    Assocs[57].Slowdef = 0;
    Assocs[57].Deltim = ILOC_NULLVAL;
    Assocs[58].StaInd = 37;    // KSRS
    Assocs[58].arid = 59;
    strcpy(Assocs[58].PhaseHint, "P");
    Assocs[58].phaseFixed = 0;
    Assocs[58].ArrivalTime = 1262324544.8;
    Assocs[58].BackAzimuth = 181.8;
    Assocs[58].Slowness = 9.90;
    Assocs[58].Timedef = 1;
    Assocs[58].Azimdef = 1;
    Assocs[58].Slowdef = 1;
    Assocs[58].Deltim = ILOC_NULLVAL;
    Assocs[59].StaInd = 37;    // KSRS
    Assocs[59].arid = 60;
    strcpy(Assocs[59].PhaseHint, "LR");
    Assocs[59].phaseFixed = 0;
    Assocs[59].ArrivalTime = 1262325033.816;
    Assocs[59].BackAzimuth = 200.0;
    Assocs[59].Slowness = 33.80;
    Assocs[59].Timedef = 1;
    Assocs[59].Azimdef = 1;
    Assocs[59].Slowdef = 1;
    Assocs[59].Deltim = ILOC_NULLVAL;
    Assocs[60].StaInd = 38;    // KMI
    Assocs[60].arid = 61;
    strcpy(Assocs[60].PhaseHint, "P");
    Assocs[60].phaseFixed = 0;
    Assocs[60].ArrivalTime = 1262324550.9;
    Assocs[60].BackAzimuth = ILOC_NULLVAL;
    Assocs[60].Slowness = ILOC_NULLVAL;
    Assocs[60].Timedef = 1;
    Assocs[60].Azimdef = 0;
    Assocs[60].Slowdef = 0;
    Assocs[60].Deltim = ILOC_NULLVAL;
    Assocs[61].StaInd = 39;    // JAGI
    Assocs[61].arid = 62;
    strcpy(Assocs[61].PhaseHint, "P");
    Assocs[61].phaseFixed = 0;
    Assocs[61].ArrivalTime = 1262324553.97;
    Assocs[61].BackAzimuth = ILOC_NULLVAL;
    Assocs[61].Slowness = ILOC_NULLVAL;
    Assocs[61].Timedef = 1;
    Assocs[61].Azimdef = 0;
    Assocs[61].Slowdef = 0;
    Assocs[61].Deltim = ILOC_NULLVAL;
    Assocs[62].StaInd = 40;    // XAN
    Assocs[62].arid = 63;
    strcpy(Assocs[62].PhaseHint, "P");
    Assocs[62].phaseFixed = 0;
    Assocs[62].ArrivalTime = 1262324557.5;
    Assocs[62].BackAzimuth = ILOC_NULLVAL;
    Assocs[62].Slowness = ILOC_NULLVAL;
    Assocs[62].Timedef = 1;
    Assocs[62].Azimdef = 0;
    Assocs[62].Slowdef = 0;
    Assocs[62].Deltim = ILOC_NULLVAL;
    Assocs[63].StaInd = 40;    // XAN
    Assocs[63].arid = 64;
    strcpy(Assocs[63].PhaseHint, "pP");
    Assocs[63].phaseFixed = 0;
    Assocs[63].ArrivalTime = 1262324560.4;
    Assocs[63].BackAzimuth = ILOC_NULLVAL;
    Assocs[63].Slowness = ILOC_NULLVAL;
    Assocs[63].Timedef = 1;
    Assocs[63].Azimdef = 0;
    Assocs[63].Slowdef = 0;
    Assocs[63].Deltim = ILOC_NULLVAL;
    Assocs[64].StaInd = 41;    // MJAR
    Assocs[64].arid = 65;
    strcpy(Assocs[64].PhaseHint, "LR");
    Assocs[64].phaseFixed = 0;
    Assocs[64].ArrivalTime = 1262325053.654;
    Assocs[64].BackAzimuth = 205.0;
    Assocs[64].Slowness = 32.40;
    Assocs[64].Timedef = 1;
    Assocs[64].Azimdef = 1;
    Assocs[64].Slowdef = 1;
    Assocs[64].Deltim = ILOC_NULLVAL;
    Assocs[65].StaInd = 42;    // CM31
    Assocs[65].arid = 66;
    strcpy(Assocs[65].PhaseHint, "P");
    Assocs[65].phaseFixed = 0;
    Assocs[65].ArrivalTime = 1262324564.93;
    Assocs[65].BackAzimuth = ILOC_NULLVAL;
    Assocs[65].Slowness = ILOC_NULLVAL;
    Assocs[65].Timedef = 1;
    Assocs[65].Azimdef = 0;
    Assocs[65].Slowdef = 0;
    Assocs[65].Deltim = ILOC_NULLVAL;
    Assocs[66].StaInd = 43;    // CMAR
    Assocs[66].arid = 67;
    strcpy(Assocs[66].PhaseHint, "P");
    Assocs[66].phaseFixed = 0;
    Assocs[66].ArrivalTime = 1262324564.85;
    Assocs[66].BackAzimuth = ILOC_NULLVAL;
    Assocs[66].Slowness = ILOC_NULLVAL;
    Assocs[66].Timedef = 1;
    Assocs[66].Azimdef = 0;
    Assocs[66].Slowdef = 0;
    Assocs[66].Deltim = ILOC_NULLVAL;
    Assocs[67].StaInd = 43;    // CMAR
    Assocs[67].arid = 68;
    strcpy(Assocs[67].PhaseHint, "P");
    Assocs[67].phaseFixed = 0;
    Assocs[67].ArrivalTime = 1262324774.65;
    Assocs[67].BackAzimuth = ILOC_NULLVAL;
    Assocs[67].Slowness = ILOC_NULLVAL;
    Assocs[67].Timedef = 1;
    Assocs[67].Azimdef = 0;
    Assocs[67].Slowdef = 0;
    Assocs[67].Deltim = ILOC_NULLVAL;
    Assocs[68].StaInd = 43;    // CMAR
    Assocs[68].arid = 69;
    strcpy(Assocs[68].PhaseHint, "PcP");
    Assocs[68].phaseFixed = 0;
    Assocs[68].ArrivalTime = 1262324564.85;
    Assocs[68].BackAzimuth = 96.7;
    Assocs[68].Slowness = 7.90;
    Assocs[68].Timedef = 1;
    Assocs[68].Azimdef = 1;
    Assocs[68].Slowdef = 1;
    Assocs[68].Deltim = ILOC_NULLVAL;
    Assocs[69].StaInd = 43;    // CMAR
    Assocs[69].arid = 70;
    strcpy(Assocs[69].PhaseHint, "PcP");
    Assocs[69].phaseFixed = 0;
    Assocs[69].ArrivalTime = 1262324774.65;
    Assocs[69].BackAzimuth = 38.3;
    Assocs[69].Slowness = 1.10;
    Assocs[69].Timedef = 1;
    Assocs[69].Azimdef = 1;
    Assocs[69].Slowdef = 1;
    Assocs[69].Deltim = ILOC_NULLVAL;
    Assocs[70].StaInd = 43;    // CMAR
    Assocs[70].arid = 71;
    strcpy(Assocs[70].PhaseHint, "LR");
    Assocs[70].phaseFixed = 0;
    Assocs[70].ArrivalTime = 1262325219.414;
    Assocs[70].BackAzimuth = 108.5;
    Assocs[70].Slowness = 38.10;
    Assocs[70].Timedef = 1;
    Assocs[70].Azimdef = 1;
    Assocs[70].Slowdef = 1;
    Assocs[70].Deltim = ILOC_NULLVAL;
    Assocs[71].StaInd = 44;    // CD2
    Assocs[71].arid = 72;
    strcpy(Assocs[71].PhaseHint, "P");
    Assocs[71].phaseFixed = 0;
    Assocs[71].ArrivalTime = 1262324569.4;
    Assocs[71].BackAzimuth = ILOC_NULLVAL;
    Assocs[71].Slowness = ILOC_NULLVAL;
    Assocs[71].Timedef = 1;
    Assocs[71].Azimdef = 0;
    Assocs[71].Slowdef = 0;
    Assocs[71].Deltim = ILOC_NULLVAL;
    Assocs[72].StaInd = 44;    // CD2
    Assocs[72].arid = 73;
    strcpy(Assocs[72].PhaseHint, "pP");
    Assocs[72].phaseFixed = 0;
    Assocs[72].ArrivalTime = 1262324573.8;
    Assocs[72].BackAzimuth = ILOC_NULLVAL;
    Assocs[72].Slowness = ILOC_NULLVAL;
    Assocs[72].Timedef = 1;
    Assocs[72].Azimdef = 0;
    Assocs[72].Slowdef = 0;
    Assocs[72].Deltim = ILOC_NULLVAL;
    Assocs[73].StaInd = 44;    // CD2
    Assocs[73].arid = 74;
    strcpy(Assocs[73].PhaseHint, "sP");
    Assocs[73].phaseFixed = 0;
    Assocs[73].ArrivalTime = 1262324576.6;
    Assocs[73].BackAzimuth = ILOC_NULLVAL;
    Assocs[73].Slowness = ILOC_NULLVAL;
    Assocs[73].Timedef = 1;
    Assocs[73].Azimdef = 0;
    Assocs[73].Slowdef = 0;
    Assocs[73].Deltim = ILOC_NULLVAL;
    Assocs[74].StaInd = 44;    // CD2
    Assocs[74].arid = 75;
    strcpy(Assocs[74].PhaseHint, "PnPn");
    Assocs[74].phaseFixed = 0;
    Assocs[74].ArrivalTime = 1262324612.7;
    Assocs[74].BackAzimuth = ILOC_NULLVAL;
    Assocs[74].Slowness = ILOC_NULLVAL;
    Assocs[74].Timedef = 1;
    Assocs[74].Azimdef = 0;
    Assocs[74].Slowdef = 0;
    Assocs[74].Deltim = ILOC_NULLVAL;
    Assocs[75].StaInd = 44;    // CD2
    Assocs[75].arid = 76;
    strcpy(Assocs[75].PhaseHint, "S");
    Assocs[75].phaseFixed = 0;
    Assocs[75].ArrivalTime = 1262324841.1;
    Assocs[75].BackAzimuth = ILOC_NULLVAL;
    Assocs[75].Slowness = ILOC_NULLVAL;
    Assocs[75].Timedef = 1;
    Assocs[75].Azimdef = 0;
    Assocs[75].Slowdef = 0;
    Assocs[75].Deltim = ILOC_NULLVAL;
    Assocs[76].StaInd = 44;    // CD2
    Assocs[76].arid = 77;
    strcpy(Assocs[76].PhaseHint, "SnSn");
    Assocs[76].phaseFixed = 0;
    Assocs[76].ArrivalTime = 1262324910.5;
    Assocs[76].BackAzimuth = ILOC_NULLVAL;
    Assocs[76].Slowness = ILOC_NULLVAL;
    Assocs[76].Timedef = 1;
    Assocs[76].Azimdef = 0;
    Assocs[76].Slowdef = 0;
    Assocs[76].Deltim = ILOC_NULLVAL;
    Assocs[77].StaInd = 45;    // BJI
    Assocs[77].arid = 78;
    strcpy(Assocs[77].PhaseHint, "P");
    Assocs[77].phaseFixed = 0;
    Assocs[77].ArrivalTime = 1262324578.6;
    Assocs[77].BackAzimuth = ILOC_NULLVAL;
    Assocs[77].Slowness = ILOC_NULLVAL;
    Assocs[77].Timedef = 1;
    Assocs[77].Azimdef = 0;
    Assocs[77].Slowdef = 0;
    Assocs[77].Deltim = ILOC_NULLVAL;
    Assocs[78].StaInd = 46;    // PSI
    Assocs[78].arid = 79;
    strcpy(Assocs[78].PhaseHint, "LR");
    Assocs[78].phaseFixed = 0;
    Assocs[78].ArrivalTime = 1262325250.685;
    Assocs[78].BackAzimuth = 141.0;
    Assocs[78].Slowness = 35.90;
    Assocs[78].Timedef = 1;
    Assocs[78].Azimdef = 1;
    Assocs[78].Slowdef = 1;
    Assocs[78].Deltim = ILOC_NULLVAL;
    Assocs[79].StaInd = 47;    // LZH
    Assocs[79].arid = 80;
    strcpy(Assocs[79].PhaseHint, "P");
    Assocs[79].phaseFixed = 0;
    Assocs[79].ArrivalTime = 1262324597.5;
    Assocs[79].BackAzimuth = ILOC_NULLVAL;
    Assocs[79].Slowness = ILOC_NULLVAL;
    Assocs[79].Timedef = 1;
    Assocs[79].Azimdef = 0;
    Assocs[79].Slowdef = 0;
    Assocs[79].Deltim = ILOC_NULLVAL;
    Assocs[80].StaInd = 47;    // LZH
    Assocs[80].arid = 81;
    strcpy(Assocs[80].PhaseHint, "pP");
    Assocs[80].phaseFixed = 0;
    Assocs[80].ArrivalTime = 1262324602.1;
    Assocs[80].BackAzimuth = ILOC_NULLVAL;
    Assocs[80].Slowness = ILOC_NULLVAL;
    Assocs[80].Timedef = 1;
    Assocs[80].Azimdef = 0;
    Assocs[80].Slowdef = 0;
    Assocs[80].Deltim = ILOC_NULLVAL;
    Assocs[81].StaInd = 47;    // LZH
    Assocs[81].arid = 82;
    strcpy(Assocs[81].PhaseHint, "sP");
    Assocs[81].phaseFixed = 0;
    Assocs[81].ArrivalTime = 1262324604.8;
    Assocs[81].BackAzimuth = ILOC_NULLVAL;
    Assocs[81].Slowness = ILOC_NULLVAL;
    Assocs[81].Timedef = 1;
    Assocs[81].Azimdef = 0;
    Assocs[81].Slowdef = 0;
    Assocs[81].Deltim = ILOC_NULLVAL;
    Assocs[82].StaInd = 48;    // HHC
    Assocs[82].arid = 83;
    strcpy(Assocs[82].PhaseHint, "P");
    Assocs[82].phaseFixed = 0;
    Assocs[82].ArrivalTime = 1262324599.5;
    Assocs[82].BackAzimuth = ILOC_NULLVAL;
    Assocs[82].Slowness = ILOC_NULLVAL;
    Assocs[82].Timedef = 1;
    Assocs[82].Azimdef = 0;
    Assocs[82].Slowdef = 0;
    Assocs[82].Deltim = ILOC_NULLVAL;
    Assocs[83].StaInd = 48;    // HHC
    Assocs[83].arid = 84;
    strcpy(Assocs[83].PhaseHint, "pP");
    Assocs[83].phaseFixed = 0;
    Assocs[83].ArrivalTime = 1262324603.8;
    Assocs[83].BackAzimuth = ILOC_NULLVAL;
    Assocs[83].Slowness = ILOC_NULLVAL;
    Assocs[83].Timedef = 1;
    Assocs[83].Azimdef = 0;
    Assocs[83].Slowdef = 0;
    Assocs[83].Deltim = ILOC_NULLVAL;
    Assocs[84].StaInd = 49;    // CN2
    Assocs[84].arid = 85;
    strcpy(Assocs[84].PhaseHint, "P");
    Assocs[84].phaseFixed = 0;
    Assocs[84].ArrivalTime = 1262324601.6;
    Assocs[84].BackAzimuth = ILOC_NULLVAL;
    Assocs[84].Slowness = ILOC_NULLVAL;
    Assocs[84].Timedef = 1;
    Assocs[84].Azimdef = 0;
    Assocs[84].Slowdef = 0;
    Assocs[84].Deltim = ILOC_NULLVAL;
    Assocs[85].StaInd = 49;    // CN2
    Assocs[85].arid = 86;
    strcpy(Assocs[85].PhaseHint, "sP");
    Assocs[85].phaseFixed = 0;
    Assocs[85].ArrivalTime = 1262324610.4;
    Assocs[85].BackAzimuth = ILOC_NULLVAL;
    Assocs[85].Slowness = ILOC_NULLVAL;
    Assocs[85].Timedef = 1;
    Assocs[85].Azimdef = 0;
    Assocs[85].Slowdef = 0;
    Assocs[85].Deltim = ILOC_NULLVAL;
    Assocs[86].StaInd = 49;    // CN2
    Assocs[86].arid = 87;
    strcpy(Assocs[86].PhaseHint, "S");
    Assocs[86].phaseFixed = 0;
    Assocs[86].ArrivalTime = 1262324897.4;
    Assocs[86].BackAzimuth = ILOC_NULLVAL;
    Assocs[86].Slowness = ILOC_NULLVAL;
    Assocs[86].Timedef = 1;
    Assocs[86].Azimdef = 0;
    Assocs[86].Slowdef = 0;
    Assocs[86].Deltim = ILOC_NULLVAL;
    Assocs[87].StaInd = 50;    // USRK
    Assocs[87].arid = 88;
    strcpy(Assocs[87].PhaseHint, "P");
    Assocs[87].phaseFixed = 0;
    Assocs[87].ArrivalTime = 1262324608.775;
    Assocs[87].BackAzimuth = 185.1;
    Assocs[87].Slowness = 9.60;
    Assocs[87].Timedef = 1;
    Assocs[87].Azimdef = 1;
    Assocs[87].Slowdef = 1;
    Assocs[87].Deltim = ILOC_NULLVAL;
    Assocs[88].StaInd = 51;    // MDJ
    Assocs[88].arid = 89;
    strcpy(Assocs[88].PhaseHint, "P");
    Assocs[88].phaseFixed = 0;
    Assocs[88].ArrivalTime = 1262324609.6;
    Assocs[88].BackAzimuth = ILOC_NULLVAL;
    Assocs[88].Slowness = ILOC_NULLVAL;
    Assocs[88].Timedef = 1;
    Assocs[88].Azimdef = 0;
    Assocs[88].Slowdef = 0;
    Assocs[88].Deltim = ILOC_NULLVAL;
    Assocs[89].StaInd = 51;    // MDJ
    Assocs[89].arid = 90;
    strcpy(Assocs[89].PhaseHint, "pP");
    Assocs[89].phaseFixed = 0;
    Assocs[89].ArrivalTime = 1262324612.6;
    Assocs[89].BackAzimuth = ILOC_NULLVAL;
    Assocs[89].Slowness = ILOC_NULLVAL;
    Assocs[89].Timedef = 1;
    Assocs[89].Azimdef = 0;
    Assocs[89].Slowdef = 0;
    Assocs[89].Deltim = ILOC_NULLVAL;
    Assocs[90].StaInd = 51;    // MDJ
    Assocs[90].arid = 91;
    strcpy(Assocs[90].PhaseHint, "sP");
    Assocs[90].phaseFixed = 0;
    Assocs[90].ArrivalTime = 1262324614;
    Assocs[90].BackAzimuth = ILOC_NULLVAL;
    Assocs[90].Slowness = ILOC_NULLVAL;
    Assocs[90].Timedef = 1;
    Assocs[90].Azimdef = 0;
    Assocs[90].Slowdef = 0;
    Assocs[90].Deltim = ILOC_NULLVAL;
    Assocs[91].StaInd = 51;    // MDJ
    Assocs[91].arid = 92;
    strcpy(Assocs[91].PhaseHint, "S");
    Assocs[91].phaseFixed = 0;
    Assocs[91].ArrivalTime = 1262324915.3;
    Assocs[91].BackAzimuth = ILOC_NULLVAL;
    Assocs[91].Slowness = ILOC_NULLVAL;
    Assocs[91].Timedef = 1;
    Assocs[91].Azimdef = 0;
    Assocs[91].Slowdef = 0;
    Assocs[91].Deltim = ILOC_NULLVAL;
    Assocs[92].StaInd = 51;    // MDJ
    Assocs[92].arid = 93;
    strcpy(Assocs[92].PhaseHint, "sS");
    Assocs[92].phaseFixed = 0;
    Assocs[92].ArrivalTime = 1262324920.5;
    Assocs[92].BackAzimuth = ILOC_NULLVAL;
    Assocs[92].Slowness = ILOC_NULLVAL;
    Assocs[92].Timedef = 1;
    Assocs[92].Azimdef = 0;
    Assocs[92].Slowdef = 0;
    Assocs[92].Deltim = ILOC_NULLVAL;
    Assocs[93].StaInd = 51;    // MDJ
    Assocs[93].arid = 94;
    strcpy(Assocs[93].PhaseHint, "ScP");
    Assocs[93].phaseFixed = 0;
    Assocs[93].ArrivalTime = 1262325008;
    Assocs[93].BackAzimuth = ILOC_NULLVAL;
    Assocs[93].Slowness = ILOC_NULLVAL;
    Assocs[93].Timedef = 1;
    Assocs[93].Azimdef = 0;
    Assocs[93].Slowdef = 0;
    Assocs[93].Deltim = ILOC_NULLVAL;
    Assocs[94].StaInd = 51;    // MDJ
    Assocs[94].arid = 95;
    strcpy(Assocs[94].PhaseHint, "PcS");
    Assocs[94].phaseFixed = 0;
    Assocs[94].ArrivalTime = 1262325009.3;
    Assocs[94].BackAzimuth = ILOC_NULLVAL;
    Assocs[94].Slowness = ILOC_NULLVAL;
    Assocs[94].Timedef = 1;
    Assocs[94].Azimdef = 0;
    Assocs[94].Slowdef = 0;
    Assocs[94].Deltim = ILOC_NULLVAL;
    Assocs[95].StaInd = 52;    // COEN
    Assocs[95].arid = 96;
    strcpy(Assocs[95].PhaseHint, "P");
    Assocs[95].phaseFixed = 0;
    Assocs[95].ArrivalTime = 1262324622.55;
    Assocs[95].BackAzimuth = ILOC_NULLVAL;
    Assocs[95].Slowness = ILOC_NULLVAL;
    Assocs[95].Timedef = 1;
    Assocs[95].Azimdef = 0;
    Assocs[95].Slowdef = 0;
    Assocs[95].Deltim = ILOC_NULLVAL;
    Assocs[96].StaInd = 53;    // SHL
    Assocs[96].arid = 97;
    strcpy(Assocs[96].PhaseHint, "P");
    Assocs[96].phaseFixed = 0;
    Assocs[96].ArrivalTime = 1262324635;
    Assocs[96].BackAzimuth = ILOC_NULLVAL;
    Assocs[96].Slowness = ILOC_NULLVAL;
    Assocs[96].Timedef = 1;
    Assocs[96].Azimdef = 0;
    Assocs[96].Slowdef = 0;
    Assocs[96].Deltim = ILOC_NULLVAL;
    Assocs[97].StaInd = 54;    // GTA
    Assocs[97].arid = 98;
    strcpy(Assocs[97].PhaseHint, "P");
    Assocs[97].phaseFixed = 0;
    Assocs[97].ArrivalTime = 1262324638.2;
    Assocs[97].BackAzimuth = ILOC_NULLVAL;
    Assocs[97].Slowness = ILOC_NULLVAL;
    Assocs[97].Timedef = 1;
    Assocs[97].Azimdef = 0;
    Assocs[97].Slowdef = 0;
    Assocs[97].Deltim = ILOC_NULLVAL;
    Assocs[98].StaInd = 54;    // GTA
    Assocs[98].arid = 99;
    strcpy(Assocs[98].PhaseHint, "pP");
    Assocs[98].phaseFixed = 0;
    Assocs[98].ArrivalTime = 1262324642.5;
    Assocs[98].BackAzimuth = ILOC_NULLVAL;
    Assocs[98].Slowness = ILOC_NULLVAL;
    Assocs[98].Timedef = 1;
    Assocs[98].Azimdef = 0;
    Assocs[98].Slowdef = 0;
    Assocs[98].Deltim = ILOC_NULLVAL;
    Assocs[99].StaInd = 54;    // GTA
    Assocs[99].arid = 100;
    strcpy(Assocs[99].PhaseHint, "sP");
    Assocs[99].phaseFixed = 0;
    Assocs[99].ArrivalTime = 1262324646;
    Assocs[99].BackAzimuth = ILOC_NULLVAL;
    Assocs[99].Slowness = ILOC_NULLVAL;
    Assocs[99].Timedef = 1;
    Assocs[99].Azimdef = 0;
    Assocs[99].Slowdef = 0;
    Assocs[99].Deltim = ILOC_NULLVAL;
    Assocs[100].StaInd = 55;    // WRAB
    Assocs[100].arid = 101;
    strcpy(Assocs[100].PhaseHint, "P");
    Assocs[100].phaseFixed = 0;
    Assocs[100].ArrivalTime = 1262324638.7;
    Assocs[100].BackAzimuth = ILOC_NULLVAL;
    Assocs[100].Slowness = ILOC_NULLVAL;
    Assocs[100].Timedef = 1;
    Assocs[100].Azimdef = 0;
    Assocs[100].Slowdef = 0;
    Assocs[100].Deltim = ILOC_NULLVAL;
    Assocs[101].StaInd = 56;    // WRA
    Assocs[101].arid = 102;
    strcpy(Assocs[101].PhaseHint, "P");
    Assocs[101].phaseFixed = 0;
    Assocs[101].ArrivalTime = 1262324639.05;
    Assocs[101].BackAzimuth = ILOC_NULLVAL;
    Assocs[101].Slowness = ILOC_NULLVAL;
    Assocs[101].Timedef = 1;
    Assocs[101].Azimdef = 0;
    Assocs[101].Slowdef = 0;
    Assocs[101].Deltim = ILOC_NULLVAL;
    Assocs[102].StaInd = 56;    // WRA
    Assocs[102].arid = 103;
    strcpy(Assocs[102].PhaseHint, "P");
    Assocs[102].phaseFixed = 0;
    Assocs[102].ArrivalTime = 1262324795.13;
    Assocs[102].BackAzimuth = ILOC_NULLVAL;
    Assocs[102].Slowness = ILOC_NULLVAL;
    Assocs[102].Timedef = 1;
    Assocs[102].Azimdef = 0;
    Assocs[102].Slowdef = 0;
    Assocs[102].Deltim = ILOC_NULLVAL;
    Assocs[103].StaInd = 56;    // WRA
    Assocs[103].arid = 104;
    strcpy(Assocs[103].PhaseHint, "PcP");
    Assocs[103].phaseFixed = 0;
    Assocs[103].ArrivalTime = 1262324639.05;
    Assocs[103].BackAzimuth = 347.3;
    Assocs[103].Slowness = 9.10;
    Assocs[103].Timedef = 1;
    Assocs[103].Azimdef = 1;
    Assocs[103].Slowdef = 1;
    Assocs[103].Deltim = ILOC_NULLVAL;
    Assocs[104].StaInd = 56;    // WRA
    Assocs[104].arid = 105;
    strcpy(Assocs[104].PhaseHint, "PcP");
    Assocs[104].phaseFixed = 0;
    Assocs[104].ArrivalTime = 1262324795.125;
    Assocs[104].BackAzimuth = 341.9;
    Assocs[104].Slowness = 3.40;
    Assocs[104].Timedef = 1;
    Assocs[104].Azimdef = 1;
    Assocs[104].Slowdef = 1;
    Assocs[104].Deltim = ILOC_NULLVAL;
    Assocs[105].StaInd = 57;    // HABR
    Assocs[105].arid = 106;
    strcpy(Assocs[105].PhaseHint, "P");
    Assocs[105].phaseFixed = 0;
    Assocs[105].ArrivalTime = 1262324655.8;
    Assocs[105].BackAzimuth = ILOC_NULLVAL;
    Assocs[105].Slowness = ILOC_NULLVAL;
    Assocs[105].Timedef = 1;
    Assocs[105].Azimdef = 0;
    Assocs[105].Slowdef = 0;
    Assocs[105].Deltim = ILOC_NULLVAL;
    Assocs[106].StaInd = 57;    // HABR
    Assocs[106].arid = 107;
    strcpy(Assocs[106].PhaseHint, "S");
    Assocs[106].phaseFixed = 0;
    Assocs[106].ArrivalTime = 1262324979.8;
    Assocs[106].BackAzimuth = ILOC_NULLVAL;
    Assocs[106].Slowness = ILOC_NULLVAL;
    Assocs[106].Timedef = 1;
    Assocs[106].Azimdef = 0;
    Assocs[106].Slowdef = 0;
    Assocs[106].Deltim = ILOC_NULLVAL;
    Assocs[107].StaInd = 57;    // HABR
    Assocs[107].arid = 108;
    strcpy(Assocs[107].PhaseHint, "sS");
    Assocs[107].phaseFixed = 0;
    Assocs[107].ArrivalTime = 1262324996.6;
    Assocs[107].BackAzimuth = ILOC_NULLVAL;
    Assocs[107].Slowness = ILOC_NULLVAL;
    Assocs[107].Timedef = 1;
    Assocs[107].Azimdef = 0;
    Assocs[107].Slowdef = 0;
    Assocs[107].Deltim = ILOC_NULLVAL;
    Assocs[108].StaInd = 57;    // HABR
    Assocs[108].arid = 109;
    strcpy(Assocs[108].PhaseHint, "SnSn");
    Assocs[108].phaseFixed = 0;
    Assocs[108].ArrivalTime = 1262325123.9;
    Assocs[108].BackAzimuth = ILOC_NULLVAL;
    Assocs[108].Slowness = ILOC_NULLVAL;
    Assocs[108].Timedef = 1;
    Assocs[108].Azimdef = 0;
    Assocs[108].Slowdef = 0;
    Assocs[108].Deltim = ILOC_NULLVAL;
    Assocs[109].StaInd = 58;    // LSA
    Assocs[109].arid = 110;
    strcpy(Assocs[109].PhaseHint, "P");
    Assocs[109].phaseFixed = 0;
    Assocs[109].ArrivalTime = 1262324649.8;
    Assocs[109].BackAzimuth = ILOC_NULLVAL;
    Assocs[109].Slowness = ILOC_NULLVAL;
    Assocs[109].Timedef = 1;
    Assocs[109].Azimdef = 0;
    Assocs[109].Slowdef = 0;
    Assocs[109].Deltim = ILOC_NULLVAL;
    Assocs[110].StaInd = 58;    // LSA
    Assocs[110].arid = 111;
    strcpy(Assocs[110].PhaseHint, "P");
    Assocs[110].phaseFixed = 0;
    Assocs[110].ArrivalTime = 1262324649.79;
    Assocs[110].BackAzimuth = ILOC_NULLVAL;
    Assocs[110].Slowness = ILOC_NULLVAL;
    Assocs[110].Timedef = 1;
    Assocs[110].Azimdef = 0;
    Assocs[110].Slowdef = 0;
    Assocs[110].Deltim = ILOC_NULLVAL;
    Assocs[111].StaInd = 59;    // ULN
    Assocs[111].arid = 112;
    strcpy(Assocs[111].PhaseHint, "P");
    Assocs[111].phaseFixed = 0;
    Assocs[111].ArrivalTime = 1262324664.83;
    Assocs[111].BackAzimuth = ILOC_NULLVAL;
    Assocs[111].Slowness = ILOC_NULLVAL;
    Assocs[111].Timedef = 1;
    Assocs[111].Azimdef = 0;
    Assocs[111].Slowdef = 0;
    Assocs[111].Deltim = ILOC_NULLVAL;
    Assocs[112].StaInd = 60;    // SONM
    Assocs[112].arid = 113;
    strcpy(Assocs[112].PhaseHint, "P");
    Assocs[112].phaseFixed = 0;
    Assocs[112].ArrivalTime = 1262324666.68;
    Assocs[112].BackAzimuth = ILOC_NULLVAL;
    Assocs[112].Slowness = ILOC_NULLVAL;
    Assocs[112].Timedef = 1;
    Assocs[112].Azimdef = 0;
    Assocs[112].Slowdef = 0;
    Assocs[112].Deltim = ILOC_NULLVAL;
    Assocs[113].StaInd = 60;    // SONM
    Assocs[113].arid = 114;
    strcpy(Assocs[113].PhaseHint, "P");
    Assocs[113].phaseFixed = 0;
    Assocs[113].ArrivalTime = 1262324804.06;
    Assocs[113].BackAzimuth = ILOC_NULLVAL;
    Assocs[113].Slowness = ILOC_NULLVAL;
    Assocs[113].Timedef = 1;
    Assocs[113].Azimdef = 0;
    Assocs[113].Slowdef = 0;
    Assocs[113].Deltim = ILOC_NULLVAL;
    Assocs[114].StaInd = 60;    // SONM
    Assocs[114].arid = 115;
    strcpy(Assocs[114].PhaseHint, "PcP");
    Assocs[114].phaseFixed = 0;
    Assocs[114].ArrivalTime = 1262324666.68;
    Assocs[114].BackAzimuth = 155.0;
    Assocs[114].Slowness = 9.70;
    Assocs[114].Timedef = 1;
    Assocs[114].Azimdef = 1;
    Assocs[114].Slowdef = 1;
    Assocs[114].Deltim = ILOC_NULLVAL;
    Assocs[115].StaInd = 60;    // SONM
    Assocs[115].arid = 116;
    strcpy(Assocs[115].PhaseHint, "PcP");
    Assocs[115].phaseFixed = 0;
    Assocs[115].ArrivalTime = 1262324804.06;
    Assocs[115].BackAzimuth = 149.4;
    Assocs[115].Slowness = 2.40;
    Assocs[115].Timedef = 1;
    Assocs[115].Azimdef = 1;
    Assocs[115].Slowdef = 1;
    Assocs[115].Deltim = ILOC_NULLVAL;
    Assocs[116].StaInd = 60;    // SONM
    Assocs[116].arid = 117;
    strcpy(Assocs[116].PhaseHint, "LR");
    Assocs[116].phaseFixed = 0;
    Assocs[116].ArrivalTime = 1262325631.635;
    Assocs[116].BackAzimuth = 163.1;
    Assocs[116].Slowness = 37.30;
    Assocs[116].Timedef = 1;
    Assocs[116].Azimdef = 1;
    Assocs[116].Slowdef = 1;
    Assocs[116].Deltim = ILOC_NULLVAL;
    Assocs[117].StaInd = 61;    // ODAN
    Assocs[117].arid = 118;
    strcpy(Assocs[117].PhaseHint, "P");
    Assocs[117].phaseFixed = 0;
    Assocs[117].ArrivalTime = 1262324670.3;
    Assocs[117].BackAzimuth = ILOC_NULLVAL;
    Assocs[117].Slowness = ILOC_NULLVAL;
    Assocs[117].Timedef = 1;
    Assocs[117].Azimdef = 0;
    Assocs[117].Slowdef = 0;
    Assocs[117].Deltim = ILOC_NULLVAL;
    Assocs[118].StaInd = 62;    // AS31
    Assocs[118].arid = 119;
    strcpy(Assocs[118].PhaseHint, "P");
    Assocs[118].phaseFixed = 0;
    Assocs[118].ArrivalTime = 1262324670.04;
    Assocs[118].BackAzimuth = ILOC_NULLVAL;
    Assocs[118].Slowness = ILOC_NULLVAL;
    Assocs[118].Timedef = 1;
    Assocs[118].Azimdef = 0;
    Assocs[118].Slowdef = 0;
    Assocs[118].Deltim = ILOC_NULLVAL;
    Assocs[119].StaInd = 63;    // ASAR
    Assocs[119].arid = 120;
    strcpy(Assocs[119].PhaseHint, "P");
    Assocs[119].phaseFixed = 0;
    Assocs[119].ArrivalTime = 1262324670;
    Assocs[119].BackAzimuth = ILOC_NULLVAL;
    Assocs[119].Slowness = ILOC_NULLVAL;
    Assocs[119].Timedef = 1;
    Assocs[119].Azimdef = 0;
    Assocs[119].Slowdef = 0;
    Assocs[119].Deltim = ILOC_NULLVAL;
    Assocs[120].StaInd = 63;    // ASAR
    Assocs[120].arid = 121;
    strcpy(Assocs[120].PhaseHint, "P");
    Assocs[120].phaseFixed = 0;
    Assocs[120].ArrivalTime = 1262324804.65;
    Assocs[120].BackAzimuth = ILOC_NULLVAL;
    Assocs[120].Slowness = ILOC_NULLVAL;
    Assocs[120].Timedef = 1;
    Assocs[120].Azimdef = 0;
    Assocs[120].Slowdef = 0;
    Assocs[120].Deltim = ILOC_NULLVAL;
    Assocs[121].StaInd = 63;    // ASAR
    Assocs[121].arid = 122;
    strcpy(Assocs[121].PhaseHint, "PcP");
    Assocs[121].phaseFixed = 0;
    Assocs[121].ArrivalTime = 1262324670;
    Assocs[121].BackAzimuth = 345.2;
    Assocs[121].Slowness = 6.90;
    Assocs[121].Timedef = 1;
    Assocs[121].Azimdef = 1;
    Assocs[121].Slowdef = 1;
    Assocs[121].Deltim = ILOC_NULLVAL;
    Assocs[122].StaInd = 63;    // ASAR
    Assocs[122].arid = 123;
    strcpy(Assocs[122].PhaseHint, "PcP");
    Assocs[122].phaseFixed = 0;
    Assocs[122].ArrivalTime = 1262324804.65;
    Assocs[122].BackAzimuth = 343.6;
    Assocs[122].Slowness = 2.30;
    Assocs[122].Timedef = 1;
    Assocs[122].Azimdef = 1;
    Assocs[122].Slowdef = 1;
    Assocs[122].Deltim = ILOC_NULLVAL;
    Assocs[123].StaInd = 64;    // RAMN
    Assocs[123].arid = 124;
    strcpy(Assocs[123].PhaseHint, "P");
    Assocs[123].phaseFixed = 0;
    Assocs[123].ArrivalTime = 1262324676.7;
    Assocs[123].BackAzimuth = ILOC_NULLVAL;
    Assocs[123].Slowness = ILOC_NULLVAL;
    Assocs[123].Timedef = 1;
    Assocs[123].Azimdef = 0;
    Assocs[123].Slowdef = 0;
    Assocs[123].Deltim = ILOC_NULLVAL;
    Assocs[124].StaInd = 65;    // JIRN
    Assocs[124].arid = 125;
    strcpy(Assocs[124].PhaseHint, "P");
    Assocs[124].phaseFixed = 0;
    Assocs[124].ArrivalTime = 1262324680.9;
    Assocs[124].BackAzimuth = ILOC_NULLVAL;
    Assocs[124].Slowness = ILOC_NULLVAL;
    Assocs[124].Timedef = 1;
    Assocs[124].Azimdef = 0;
    Assocs[124].Slowdef = 0;
    Assocs[124].Deltim = ILOC_NULLVAL;
    Assocs[125].StaInd = 66;    // CTA
    Assocs[125].arid = 126;
    strcpy(Assocs[125].PhaseHint, "P");
    Assocs[125].phaseFixed = 0;
    Assocs[125].ArrivalTime = 1262324681.4;
    Assocs[125].BackAzimuth = ILOC_NULLVAL;
    Assocs[125].Slowness = ILOC_NULLVAL;
    Assocs[125].Timedef = 1;
    Assocs[125].Azimdef = 0;
    Assocs[125].Slowdef = 0;
    Assocs[125].Deltim = ILOC_NULLVAL;
    Assocs[126].StaInd = 66;    // CTA
    Assocs[126].arid = 127;
    strcpy(Assocs[126].PhaseHint, "P");
    Assocs[126].phaseFixed = 0;
    Assocs[126].ArrivalTime = 1262324681.4;
    Assocs[126].BackAzimuth = 334.9;
    Assocs[126].Slowness = 8.40;
    Assocs[126].Timedef = 1;
    Assocs[126].Azimdef = 1;
    Assocs[126].Slowdef = 1;
    Assocs[126].Deltim = ILOC_NULLVAL;
    Assocs[127].StaInd = 67;    // GUN
    Assocs[127].arid = 128;
    strcpy(Assocs[127].PhaseHint, "P");
    Assocs[127].phaseFixed = 0;
    Assocs[127].ArrivalTime = 1262324683.2;
    Assocs[127].BackAzimuth = ILOC_NULLVAL;
    Assocs[127].Slowness = ILOC_NULLVAL;
    Assocs[127].Timedef = 1;
    Assocs[127].Azimdef = 0;
    Assocs[127].Slowdef = 0;
    Assocs[127].Deltim = ILOC_NULLVAL;
    Assocs[128].StaInd = 68;    // PKI
    Assocs[128].arid = 129;
    strcpy(Assocs[128].PhaseHint, "P");
    Assocs[128].phaseFixed = 0;
    Assocs[128].ArrivalTime = 1262324686.5;
    Assocs[128].BackAzimuth = ILOC_NULLVAL;
    Assocs[128].Slowness = ILOC_NULLVAL;
    Assocs[128].Timedef = 1;
    Assocs[128].Azimdef = 0;
    Assocs[128].Slowdef = 0;
    Assocs[128].Deltim = ILOC_NULLVAL;
    Assocs[129].StaInd = 69;    // PKIN
    Assocs[129].arid = 130;
    strcpy(Assocs[129].PhaseHint, "P");
    Assocs[129].phaseFixed = 0;
    Assocs[129].ArrivalTime = 1262324686.9;
    Assocs[129].BackAzimuth = ILOC_NULLVAL;
    Assocs[129].Slowness = ILOC_NULLVAL;
    Assocs[129].Timedef = 1;
    Assocs[129].Azimdef = 0;
    Assocs[129].Slowdef = 0;
    Assocs[129].Deltim = ILOC_NULLVAL;
    Assocs[130].StaInd = 70;    // KKN
    Assocs[130].arid = 131;
    strcpy(Assocs[130].PhaseHint, "P");
    Assocs[130].phaseFixed = 0;
    Assocs[130].ArrivalTime = 1262324686.7;
    Assocs[130].BackAzimuth = ILOC_NULLVAL;
    Assocs[130].Slowness = ILOC_NULLVAL;
    Assocs[130].Timedef = 1;
    Assocs[130].Azimdef = 0;
    Assocs[130].Slowdef = 0;
    Assocs[130].Deltim = ILOC_NULLVAL;
    Assocs[131].StaInd = 71;    // DMN
    Assocs[131].arid = 132;
    strcpy(Assocs[131].PhaseHint, "P");
    Assocs[131].phaseFixed = 0;
    Assocs[131].ArrivalTime = 1262324688.4;
    Assocs[131].BackAzimuth = ILOC_NULLVAL;
    Assocs[131].Slowness = ILOC_NULLVAL;
    Assocs[131].Timedef = 1;
    Assocs[131].Azimdef = 0;
    Assocs[131].Slowdef = 0;
    Assocs[131].Deltim = ILOC_NULLVAL;
    Assocs[132].StaInd = 72;    // GKN
    Assocs[132].arid = 133;
    strcpy(Assocs[132].PhaseHint, "P");
    Assocs[132].phaseFixed = 0;
    Assocs[132].ArrivalTime = 1262324692;
    Assocs[132].BackAzimuth = ILOC_NULLVAL;
    Assocs[132].Slowness = ILOC_NULLVAL;
    Assocs[132].Timedef = 1;
    Assocs[132].Azimdef = 0;
    Assocs[132].Slowdef = 0;
    Assocs[132].Deltim = ILOC_NULLVAL;
    Assocs[133].StaInd = 73;    // ZAK
    Assocs[133].arid = 134;
    strcpy(Assocs[133].PhaseHint, "P");
    Assocs[133].phaseFixed = 0;
    Assocs[133].ArrivalTime = 1262324691.1;
    Assocs[133].BackAzimuth = ILOC_NULLVAL;
    Assocs[133].Slowness = ILOC_NULLVAL;
    Assocs[133].Timedef = 1;
    Assocs[133].Azimdef = 0;
    Assocs[133].Slowdef = 0;
    Assocs[133].Deltim = ILOC_NULLVAL;
    Assocs[134].StaInd = 74;    // HNR
    Assocs[134].arid = 135;
    strcpy(Assocs[134].PhaseHint, "LR");
    Assocs[134].phaseFixed = 0;
    Assocs[134].ArrivalTime = 1262325503.043;
    Assocs[134].BackAzimuth = 133.8;
    Assocs[134].Slowness = 30.90;
    Assocs[134].Timedef = 1;
    Assocs[134].Azimdef = 1;
    Assocs[134].Slowdef = 1;
    Assocs[134].Deltim = ILOC_NULLVAL;
    Assocs[135].StaInd = 75;    // DANN
    Assocs[135].arid = 136;
    strcpy(Assocs[135].PhaseHint, "P");
    Assocs[135].phaseFixed = 0;
    Assocs[135].ArrivalTime = 1262324698.3;
    Assocs[135].BackAzimuth = ILOC_NULLVAL;
    Assocs[135].Slowness = ILOC_NULLVAL;
    Assocs[135].Timedef = 1;
    Assocs[135].Azimdef = 0;
    Assocs[135].Slowdef = 0;
    Assocs[135].Deltim = ILOC_NULLVAL;
    Assocs[136].StaInd = 76;    // KOLN
    Assocs[136].arid = 137;
    strcpy(Assocs[136].PhaseHint, "P");
    Assocs[136].phaseFixed = 0;
    Assocs[136].ArrivalTime = 1262324698.9;
    Assocs[136].BackAzimuth = ILOC_NULLVAL;
    Assocs[136].Slowness = ILOC_NULLVAL;
    Assocs[136].Timedef = 1;
    Assocs[136].Azimdef = 0;
    Assocs[136].Slowdef = 0;
    Assocs[136].Deltim = ILOC_NULLVAL;
    Assocs[137].StaInd = 77;    // TLY
    Assocs[137].arid = 138;
    strcpy(Assocs[137].PhaseHint, "P");
    Assocs[137].phaseFixed = 0;
    Assocs[137].ArrivalTime = 1262324701.9;
    Assocs[137].BackAzimuth = ILOC_NULLVAL;
    Assocs[137].Slowness = ILOC_NULLVAL;
    Assocs[137].Timedef = 1;
    Assocs[137].Azimdef = 0;
    Assocs[137].Slowdef = 0;
    Assocs[137].Deltim = ILOC_NULLVAL;
    Assocs[138].StaInd = 78;    // WMQ
    Assocs[138].arid = 139;
    strcpy(Assocs[138].PhaseHint, "P");
    Assocs[138].phaseFixed = 0;
    Assocs[138].ArrivalTime = 1262324719.2;
    Assocs[138].BackAzimuth = ILOC_NULLVAL;
    Assocs[138].Slowness = ILOC_NULLVAL;
    Assocs[138].Timedef = 1;
    Assocs[138].Azimdef = 0;
    Assocs[138].Slowdef = 0;
    Assocs[138].Deltim = ILOC_NULLVAL;
    Assocs[139].StaInd = 78;    // WMQ
    Assocs[139].arid = 140;
    strcpy(Assocs[139].PhaseHint, "PcP");
    Assocs[139].phaseFixed = 0;
    Assocs[139].ArrivalTime = 1262324824.8;
    Assocs[139].BackAzimuth = ILOC_NULLVAL;
    Assocs[139].Slowness = ILOC_NULLVAL;
    Assocs[139].Timedef = 1;
    Assocs[139].Azimdef = 0;
    Assocs[139].Slowdef = 0;
    Assocs[139].Deltim = ILOC_NULLVAL;
    Assocs[140].StaInd = 78;    // WMQ
    Assocs[140].arid = 141;
    strcpy(Assocs[140].PhaseHint, "ScP");
    Assocs[140].phaseFixed = 0;
    Assocs[140].ArrivalTime = 1262325056.1;
    Assocs[140].BackAzimuth = ILOC_NULLVAL;
    Assocs[140].Slowness = ILOC_NULLVAL;
    Assocs[140].Timedef = 1;
    Assocs[140].Azimdef = 0;
    Assocs[140].Slowdef = 0;
    Assocs[140].Deltim = ILOC_NULLVAL;
    Assocs[141].StaInd = 78;    // WMQ
    Assocs[141].arid = 142;
    strcpy(Assocs[141].PhaseHint, "S");
    Assocs[141].phaseFixed = 0;
    Assocs[141].ArrivalTime = 1262325112;
    Assocs[141].BackAzimuth = ILOC_NULLVAL;
    Assocs[141].Slowness = ILOC_NULLVAL;
    Assocs[141].Timedef = 1;
    Assocs[141].Azimdef = 0;
    Assocs[141].Slowdef = 0;
    Assocs[141].Deltim = ILOC_NULLVAL;
    Assocs[142].StaInd = 79;    // FORT
    Assocs[142].arid = 143;
    strcpy(Assocs[142].PhaseHint, "P");
    Assocs[142].phaseFixed = 0;
    Assocs[142].ArrivalTime = 1262324721.88;
    Assocs[142].BackAzimuth = ILOC_NULLVAL;
    Assocs[142].Slowness = ILOC_NULLVAL;
    Assocs[142].Timedef = 1;
    Assocs[142].Azimdef = 0;
    Assocs[142].Slowdef = 0;
    Assocs[142].Deltim = ILOC_NULLVAL;
    Assocs[143].StaInd = 80;    // HVS
    Assocs[143].arid = 144;
    strcpy(Assocs[143].PhaseHint, "P");
    Assocs[143].phaseFixed = 0;
    Assocs[143].ArrivalTime = 1262324730.8;
    Assocs[143].BackAzimuth = ILOC_NULLVAL;
    Assocs[143].Slowness = ILOC_NULLVAL;
    Assocs[143].Timedef = 1;
    Assocs[143].Azimdef = 0;
    Assocs[143].Slowdef = 0;
    Assocs[143].Deltim = ILOC_NULLVAL;
    Assocs[144].StaInd = 81;    // EIDS
    Assocs[144].arid = 145;
    strcpy(Assocs[144].PhaseHint, "P");
    Assocs[144].phaseFixed = 0;
    Assocs[144].ArrivalTime = 1262324736.49;
    Assocs[144].BackAzimuth = ILOC_NULLVAL;
    Assocs[144].Slowness = ILOC_NULLVAL;
    Assocs[144].Timedef = 1;
    Assocs[144].Azimdef = 0;
    Assocs[144].Slowdef = 0;
    Assocs[144].Deltim = ILOC_NULLVAL;
    Assocs[145].StaInd = 82;    // PETK
    Assocs[145].arid = 146;
    strcpy(Assocs[145].PhaseHint, "P");
    Assocs[145].phaseFixed = 0;
    Assocs[145].ArrivalTime = 1262324742.4;
    Assocs[145].BackAzimuth = ILOC_NULLVAL;
    Assocs[145].Slowness = ILOC_NULLVAL;
    Assocs[145].Timedef = 1;
    Assocs[145].Azimdef = 0;
    Assocs[145].Slowdef = 0;
    Assocs[145].Deltim = ILOC_NULLVAL;
    Assocs[146].StaInd = 82;    // PETK
    Assocs[146].arid = 147;
    strcpy(Assocs[146].PhaseHint, "P");
    Assocs[146].phaseFixed = 0;
    Assocs[146].ArrivalTime = 1262324742.35;
    Assocs[146].BackAzimuth = ILOC_NULLVAL;
    Assocs[146].Slowness = ILOC_NULLVAL;
    Assocs[146].Timedef = 1;
    Assocs[146].Azimdef = 0;
    Assocs[146].Slowdef = 0;
    Assocs[146].Deltim = ILOC_NULLVAL;
    Assocs[147].StaInd = 82;    // PETK
    Assocs[147].arid = 148;
    strcpy(Assocs[147].PhaseHint, "P");
    Assocs[147].phaseFixed = 0;
    Assocs[147].ArrivalTime = 1262324742.355;
    Assocs[147].BackAzimuth = 186.3;
    Assocs[147].Slowness = 6.50;
    Assocs[147].Timedef = 1;
    Assocs[147].Azimdef = 1;
    Assocs[147].Slowdef = 1;
    Assocs[147].Deltim = ILOC_NULLVAL;
    Assocs[148].StaInd = 82;    // PETK
    Assocs[148].arid = 149;
    strcpy(Assocs[148].PhaseHint, "LR");
    Assocs[148].phaseFixed = 0;
    Assocs[148].ArrivalTime = 1262325734.503;
    Assocs[148].BackAzimuth = 20.1;
    Assocs[148].Slowness = 32.10;
    Assocs[148].Timedef = 1;
    Assocs[148].Azimdef = 1;
    Assocs[148].Slowdef = 1;
    Assocs[148].Deltim = ILOC_NULLVAL;
    Assocs[149].StaInd = 83;    // NWAO
    Assocs[149].arid = 150;
    strcpy(Assocs[149].PhaseHint, "P");
    Assocs[149].phaseFixed = 0;
    Assocs[149].ArrivalTime = 1262324742.55;
    Assocs[149].BackAzimuth = ILOC_NULLVAL;
    Assocs[149].Slowness = ILOC_NULLVAL;
    Assocs[149].Timedef = 1;
    Assocs[149].Azimdef = 0;
    Assocs[149].Slowdef = 0;
    Assocs[149].Deltim = ILOC_NULLVAL;
    Assocs[150].StaInd = 83;    // NWAO
    Assocs[150].arid = 151;
    strcpy(Assocs[150].PhaseHint, "P");
    Assocs[150].phaseFixed = 0;
    Assocs[150].ArrivalTime = 1262324742.55;
    Assocs[150].BackAzimuth = 38.7;
    Assocs[150].Slowness = 6.50;
    Assocs[150].Timedef = 1;
    Assocs[150].Azimdef = 1;
    Assocs[150].Slowdef = 1;
    Assocs[150].Deltim = ILOC_NULLVAL;
    Assocs[151].StaInd = 84;    // BBOO
    Assocs[151].arid = 152;
    strcpy(Assocs[151].PhaseHint, "P");
    Assocs[151].phaseFixed = 0;
    Assocs[151].ArrivalTime = 1262324745.68;
    Assocs[151].BackAzimuth = ILOC_NULLVAL;
    Assocs[151].Slowness = ILOC_NULLVAL;
    Assocs[151].Timedef = 1;
    Assocs[151].Azimdef = 0;
    Assocs[151].Slowdef = 0;
    Assocs[151].Deltim = ILOC_NULLVAL;
    Assocs[152].StaInd = 85;    // STKA
    Assocs[152].arid = 153;
    strcpy(Assocs[152].PhaseHint, "P");
    Assocs[152].phaseFixed = 0;
    Assocs[152].ArrivalTime = 1262324750.38;
    Assocs[152].BackAzimuth = ILOC_NULLVAL;
    Assocs[152].Slowness = ILOC_NULLVAL;
    Assocs[152].Timedef = 1;
    Assocs[152].Azimdef = 0;
    Assocs[152].Slowdef = 0;
    Assocs[152].Deltim = ILOC_NULLVAL;
    Assocs[153].StaInd = 85;    // STKA
    Assocs[153].arid = 154;
    strcpy(Assocs[153].PhaseHint, "P");
    Assocs[153].phaseFixed = 0;
    Assocs[153].ArrivalTime = 1262324750.55;
    Assocs[153].BackAzimuth = 334.7;
    Assocs[153].Slowness = 7.20;
    Assocs[153].Timedef = 1;
    Assocs[153].Azimdef = 1;
    Assocs[153].Slowdef = 1;
    Assocs[153].Deltim = ILOC_NULLVAL;
    Assocs[154].StaInd = 86;    // MK31
    Assocs[154].arid = 155;
    strcpy(Assocs[154].PhaseHint, "P");
    Assocs[154].phaseFixed = 0;
    Assocs[154].ArrivalTime = 1262324758.2;
    Assocs[154].BackAzimuth = ILOC_NULLVAL;
    Assocs[154].Slowness = ILOC_NULLVAL;
    Assocs[154].Timedef = 1;
    Assocs[154].Azimdef = 0;
    Assocs[154].Slowdef = 0;
    Assocs[154].Deltim = ILOC_NULLVAL;
    Assocs[155].StaInd = 86;    // MK31
    Assocs[155].arid = 156;
    strcpy(Assocs[155].PhaseHint, "P");
    Assocs[155].phaseFixed = 0;
    Assocs[155].ArrivalTime = 1262324758.15;
    Assocs[155].BackAzimuth = ILOC_NULLVAL;
    Assocs[155].Slowness = ILOC_NULLVAL;
    Assocs[155].Timedef = 1;
    Assocs[155].Azimdef = 0;
    Assocs[155].Slowdef = 0;
    Assocs[155].Deltim = ILOC_NULLVAL;
    Assocs[156].StaInd = 87;    // MKAR
    Assocs[156].arid = 157;
    strcpy(Assocs[156].PhaseHint, "P");
    Assocs[156].phaseFixed = 0;
    Assocs[156].ArrivalTime = 1262324758.43;
    Assocs[156].BackAzimuth = ILOC_NULLVAL;
    Assocs[156].Slowness = ILOC_NULLVAL;
    Assocs[156].Timedef = 1;
    Assocs[156].Azimdef = 0;
    Assocs[156].Slowdef = 0;
    Assocs[156].Deltim = ILOC_NULLVAL;
    Assocs[157].StaInd = 87;    // MKAR
    Assocs[157].arid = 158;
    strcpy(Assocs[157].PhaseHint, "P");
    Assocs[157].phaseFixed = 0;
    Assocs[157].ArrivalTime = 1262324758.425;
    Assocs[157].BackAzimuth = 112.1;
    Assocs[157].Slowness = 8.50;
    Assocs[157].Timedef = 1;
    Assocs[157].Azimdef = 1;
    Assocs[157].Slowdef = 1;
    Assocs[157].Deltim = ILOC_NULLVAL;
    Assocs[158].StaInd = 88;    // KSH
    Assocs[158].arid = 159;
    strcpy(Assocs[158].PhaseHint, "P");
    Assocs[158].phaseFixed = 0;
    Assocs[158].ArrivalTime = 1262324772.4;
    Assocs[158].BackAzimuth = ILOC_NULLVAL;
    Assocs[158].Slowness = ILOC_NULLVAL;
    Assocs[158].Timedef = 1;
    Assocs[158].Azimdef = 0;
    Assocs[158].Slowdef = 0;
    Assocs[158].Deltim = ILOC_NULLVAL;
    Assocs[159].StaInd = 88;    // KSH
    Assocs[159].arid = 160;
    strcpy(Assocs[159].PhaseHint, "sP");
    Assocs[159].phaseFixed = 0;
    Assocs[159].ArrivalTime = 1262324779.3;
    Assocs[159].BackAzimuth = ILOC_NULLVAL;
    Assocs[159].Slowness = ILOC_NULLVAL;
    Assocs[159].Timedef = 1;
    Assocs[159].Azimdef = 0;
    Assocs[159].Slowdef = 0;
    Assocs[159].Deltim = ILOC_NULLVAL;
    Assocs[160].StaInd = 88;    // KSH
    Assocs[160].arid = 161;
    strcpy(Assocs[160].PhaseHint, "PcP");
    Assocs[160].phaseFixed = 0;
    Assocs[160].ArrivalTime = 1262324849.1;
    Assocs[160].BackAzimuth = ILOC_NULLVAL;
    Assocs[160].Slowness = ILOC_NULLVAL;
    Assocs[160].Timedef = 1;
    Assocs[160].Azimdef = 0;
    Assocs[160].Slowdef = 0;
    Assocs[160].Deltim = ILOC_NULLVAL;
    Assocs[161].StaInd = 88;    // KSH
    Assocs[161].arid = 162;
    strcpy(Assocs[161].PhaseHint, "PP");
    Assocs[161].phaseFixed = 0;
    Assocs[161].ArrivalTime = 1262324888.7;
    Assocs[161].BackAzimuth = ILOC_NULLVAL;
    Assocs[161].Slowness = ILOC_NULLVAL;
    Assocs[161].Timedef = 1;
    Assocs[161].Azimdef = 0;
    Assocs[161].Slowdef = 0;
    Assocs[161].Deltim = ILOC_NULLVAL;
    Assocs[162].StaInd = 88;    // KSH
    Assocs[162].arid = 163;
    strcpy(Assocs[162].PhaseHint, "S");
    Assocs[162].phaseFixed = 0;
    Assocs[162].ArrivalTime = 1262325206.1;
    Assocs[162].BackAzimuth = ILOC_NULLVAL;
    Assocs[162].Slowness = ILOC_NULLVAL;
    Assocs[162].Timedef = 1;
    Assocs[162].Azimdef = 0;
    Assocs[162].Slowdef = 0;
    Assocs[162].Deltim = ILOC_NULLVAL;
    Assocs[163].StaInd = 88;    // KSH
    Assocs[163].arid = 164;
    strcpy(Assocs[163].PhaseHint, "SS");
    Assocs[163].phaseFixed = 0;
    Assocs[163].ArrivalTime = 1262325418.3;
    Assocs[163].BackAzimuth = ILOC_NULLVAL;
    Assocs[163].Slowness = ILOC_NULLVAL;
    Assocs[163].Timedef = 1;
    Assocs[163].Azimdef = 0;
    Assocs[163].Slowdef = 0;
    Assocs[163].Deltim = ILOC_NULLVAL;
    Assocs[164].StaInd = 89;    // ULHL
    Assocs[164].arid = 165;
    strcpy(Assocs[164].PhaseHint, "P");
    Assocs[164].phaseFixed = 0;
    Assocs[164].ArrivalTime = 1262324775.8;
    Assocs[164].BackAzimuth = ILOC_NULLVAL;
    Assocs[164].Slowness = ILOC_NULLVAL;
    Assocs[164].Timedef = 1;
    Assocs[164].Azimdef = 0;
    Assocs[164].Slowdef = 0;
    Assocs[164].Deltim = ILOC_NULLVAL;
    Assocs[165].StaInd = 89;    // ULHL
    Assocs[165].arid = 166;
    strcpy(Assocs[165].PhaseHint, "P");
    Assocs[165].phaseFixed = 0;
    Assocs[165].ArrivalTime = 1262324775.84;
    Assocs[165].BackAzimuth = ILOC_NULLVAL;
    Assocs[165].Slowness = ILOC_NULLVAL;
    Assocs[165].Timedef = 1;
    Assocs[165].Azimdef = 0;
    Assocs[165].Slowdef = 0;
    Assocs[165].Deltim = ILOC_NULLVAL;
    Assocs[166].StaInd = 90;    // ZAA0
    Assocs[166].arid = 167;
    strcpy(Assocs[166].PhaseHint, "P");
    Assocs[166].phaseFixed = 0;
    Assocs[166].ArrivalTime = 1262324775.18;
    Assocs[166].BackAzimuth = ILOC_NULLVAL;
    Assocs[166].Slowness = ILOC_NULLVAL;
    Assocs[166].Timedef = 1;
    Assocs[166].Azimdef = 0;
    Assocs[166].Slowdef = 0;
    Assocs[166].Deltim = ILOC_NULLVAL;
    Assocs[167].StaInd = 91;    // ZALV
    Assocs[167].arid = 168;
    strcpy(Assocs[167].PhaseHint, "P");
    Assocs[167].phaseFixed = 0;
    Assocs[167].ArrivalTime = 1262324775.4;
    Assocs[167].BackAzimuth = ILOC_NULLVAL;
    Assocs[167].Slowness = ILOC_NULLVAL;
    Assocs[167].Timedef = 1;
    Assocs[167].Azimdef = 0;
    Assocs[167].Slowdef = 0;
    Assocs[167].Deltim = ILOC_NULLVAL;
    Assocs[168].StaInd = 91;    // ZALV
    Assocs[168].arid = 169;
    strcpy(Assocs[168].PhaseHint, "P");
    Assocs[168].phaseFixed = 0;
    Assocs[168].ArrivalTime = 1262324775.43;
    Assocs[168].BackAzimuth = ILOC_NULLVAL;
    Assocs[168].Slowness = ILOC_NULLVAL;
    Assocs[168].Timedef = 1;
    Assocs[168].Azimdef = 0;
    Assocs[168].Slowdef = 0;
    Assocs[168].Deltim = ILOC_NULLVAL;
    Assocs[169].StaInd = 91;    // ZALV
    Assocs[169].arid = 170;
    strcpy(Assocs[169].PhaseHint, "P");
    Assocs[169].phaseFixed = 0;
    Assocs[169].ArrivalTime = 1262324775.425;
    Assocs[169].BackAzimuth = 117.0;
    Assocs[169].Slowness = 8.50;
    Assocs[169].Timedef = 1;
    Assocs[169].Azimdef = 1;
    Assocs[169].Slowdef = 1;
    Assocs[169].Deltim = ILOC_NULLVAL;
    Assocs[170].StaInd = 92;    // KZA
    Assocs[170].arid = 171;
    strcpy(Assocs[170].PhaseHint, "P");
    Assocs[170].phaseFixed = 0;
    Assocs[170].ArrivalTime = 1262324781.7;
    Assocs[170].BackAzimuth = ILOC_NULLVAL;
    Assocs[170].Slowness = ILOC_NULLVAL;
    Assocs[170].Timedef = 1;
    Assocs[170].Azimdef = 0;
    Assocs[170].Slowdef = 0;
    Assocs[170].Deltim = ILOC_NULLVAL;
    Assocs[171].StaInd = 92;    // KZA
    Assocs[171].arid = 172;
    strcpy(Assocs[171].PhaseHint, "P");
    Assocs[171].phaseFixed = 0;
    Assocs[171].ArrivalTime = 1262324781.79;
    Assocs[171].BackAzimuth = ILOC_NULLVAL;
    Assocs[171].Slowness = ILOC_NULLVAL;
    Assocs[171].Timedef = 1;
    Assocs[171].Azimdef = 0;
    Assocs[171].Slowdef = 0;
    Assocs[171].Deltim = ILOC_NULLVAL;
    Assocs[172].StaInd = 93;    // TKM2
    Assocs[172].arid = 173;
    strcpy(Assocs[172].PhaseHint, "P");
    Assocs[172].phaseFixed = 0;
    Assocs[172].ArrivalTime = 1262324779.9;
    Assocs[172].BackAzimuth = ILOC_NULLVAL;
    Assocs[172].Slowness = ILOC_NULLVAL;
    Assocs[172].Timedef = 1;
    Assocs[172].Azimdef = 0;
    Assocs[172].Slowdef = 0;
    Assocs[172].Deltim = ILOC_NULLVAL;
    Assocs[173].StaInd = 93;    // TKM2
    Assocs[173].arid = 174;
    strcpy(Assocs[173].PhaseHint, "P");
    Assocs[173].phaseFixed = 0;
    Assocs[173].ArrivalTime = 1262324781.03;
    Assocs[173].BackAzimuth = ILOC_NULLVAL;
    Assocs[173].Slowness = ILOC_NULLVAL;
    Assocs[173].Timedef = 1;
    Assocs[173].Azimdef = 0;
    Assocs[173].Slowdef = 0;
    Assocs[173].Deltim = ILOC_NULLVAL;
    Assocs[174].StaInd = 94;    // UCH
    Assocs[174].arid = 175;
    strcpy(Assocs[174].PhaseHint, "P");
    Assocs[174].phaseFixed = 0;
    Assocs[174].ArrivalTime = 1262324785.6;
    Assocs[174].BackAzimuth = ILOC_NULLVAL;
    Assocs[174].Slowness = ILOC_NULLVAL;
    Assocs[174].Timedef = 1;
    Assocs[174].Azimdef = 0;
    Assocs[174].Slowdef = 0;
    Assocs[174].Deltim = ILOC_NULLVAL;
    Assocs[175].StaInd = 94;    // UCH
    Assocs[175].arid = 176;
    strcpy(Assocs[175].PhaseHint, "P");
    Assocs[175].phaseFixed = 0;
    Assocs[175].ArrivalTime = 1262324785.61;
    Assocs[175].BackAzimuth = ILOC_NULLVAL;
    Assocs[175].Slowness = ILOC_NULLVAL;
    Assocs[175].Timedef = 1;
    Assocs[175].Azimdef = 0;
    Assocs[175].Slowdef = 0;
    Assocs[175].Deltim = ILOC_NULLVAL;
    Assocs[176].StaInd = 95;    // CHMS
    Assocs[176].arid = 177;
    strcpy(Assocs[176].PhaseHint, "P");
    Assocs[176].phaseFixed = 0;
    Assocs[176].ArrivalTime = 1262324784.8;
    Assocs[176].BackAzimuth = ILOC_NULLVAL;
    Assocs[176].Slowness = ILOC_NULLVAL;
    Assocs[176].Timedef = 1;
    Assocs[176].Azimdef = 0;
    Assocs[176].Slowdef = 0;
    Assocs[176].Deltim = ILOC_NULLVAL;
    Assocs[177].StaInd = 95;    // CHMS
    Assocs[177].arid = 178;
    strcpy(Assocs[177].PhaseHint, "P");
    Assocs[177].phaseFixed = 0;
    Assocs[177].ArrivalTime = 1262324784.89;
    Assocs[177].BackAzimuth = ILOC_NULLVAL;
    Assocs[177].Slowness = ILOC_NULLVAL;
    Assocs[177].Timedef = 1;
    Assocs[177].Azimdef = 0;
    Assocs[177].Slowdef = 0;
    Assocs[177].Deltim = ILOC_NULLVAL;
    Assocs[178].StaInd = 96;    // USP
    Assocs[178].arid = 179;
    strcpy(Assocs[178].PhaseHint, "P");
    Assocs[178].phaseFixed = 0;
    Assocs[178].ArrivalTime = 1262324786.6;
    Assocs[178].BackAzimuth = ILOC_NULLVAL;
    Assocs[178].Slowness = ILOC_NULLVAL;
    Assocs[178].Timedef = 1;
    Assocs[178].Azimdef = 0;
    Assocs[178].Slowdef = 0;
    Assocs[178].Deltim = ILOC_NULLVAL;
    Assocs[179].StaInd = 96;    // USP
    Assocs[179].arid = 180;
    strcpy(Assocs[179].PhaseHint, "P");
    Assocs[179].phaseFixed = 0;
    Assocs[179].ArrivalTime = 1262324786.6;
    Assocs[179].BackAzimuth = ILOC_NULLVAL;
    Assocs[179].Slowness = ILOC_NULLVAL;
    Assocs[179].Timedef = 1;
    Assocs[179].Azimdef = 0;
    Assocs[179].Slowdef = 0;
    Assocs[179].Deltim = ILOC_NULLVAL;
    Assocs[180].StaInd = 97;    // SEY
    Assocs[180].arid = 181;
    strcpy(Assocs[180].PhaseHint, "P");
    Assocs[180].phaseFixed = 0;
    Assocs[180].ArrivalTime = 1262324786.6;
    Assocs[180].BackAzimuth = 178.9;
    Assocs[180].Slowness = 23.40;
    Assocs[180].Timedef = 1;
    Assocs[180].Azimdef = 1;
    Assocs[180].Slowdef = 1;
    Assocs[180].Deltim = ILOC_NULLVAL;
    Assocs[181].StaInd = 97;    // SEY
    Assocs[181].arid = 182;
    strcpy(Assocs[181].PhaseHint, "P");
    Assocs[181].phaseFixed = 0;
    Assocs[181].ArrivalTime = 1262324786.5;
    Assocs[181].BackAzimuth = ILOC_NULLVAL;
    Assocs[181].Slowness = ILOC_NULLVAL;
    Assocs[181].Timedef = 1;
    Assocs[181].Azimdef = 0;
    Assocs[181].Slowdef = 0;
    Assocs[181].Deltim = ILOC_NULLVAL;
    Assocs[182].StaInd = 97;    // SEY
    Assocs[182].arid = 183;
    strcpy(Assocs[182].PhaseHint, "P");
    Assocs[182].phaseFixed = 0;
    Assocs[182].ArrivalTime = 1262324786.57;
    Assocs[182].BackAzimuth = ILOC_NULLVAL;
    Assocs[182].Slowness = ILOC_NULLVAL;
    Assocs[182].Timedef = 1;
    Assocs[182].Azimdef = 0;
    Assocs[182].Slowdef = 0;
    Assocs[182].Deltim = ILOC_NULLVAL;
    Assocs[183].StaInd = 98;    // AML
    Assocs[183].arid = 184;
    strcpy(Assocs[183].PhaseHint, "P");
    Assocs[183].phaseFixed = 0;
    Assocs[183].ArrivalTime = 1262324789.2;
    Assocs[183].BackAzimuth = ILOC_NULLVAL;
    Assocs[183].Slowness = ILOC_NULLVAL;
    Assocs[183].Timedef = 1;
    Assocs[183].Azimdef = 0;
    Assocs[183].Slowdef = 0;
    Assocs[183].Deltim = ILOC_NULLVAL;
    Assocs[184].StaInd = 98;    // AML
    Assocs[184].arid = 185;
    strcpy(Assocs[184].PhaseHint, "P");
    Assocs[184].phaseFixed = 0;
    Assocs[184].ArrivalTime = 1262324789.19;
    Assocs[184].BackAzimuth = ILOC_NULLVAL;
    Assocs[184].Slowness = ILOC_NULLVAL;
    Assocs[184].Timedef = 1;
    Assocs[184].Azimdef = 0;
    Assocs[184].Slowdef = 0;
    Assocs[184].Deltim = ILOC_NULLVAL;
    Assocs[185].StaInd = 99;    // EKS2
    Assocs[185].arid = 186;
    strcpy(Assocs[185].PhaseHint, "P");
    Assocs[185].phaseFixed = 0;
    Assocs[185].ArrivalTime = 1262324789.3;
    Assocs[185].BackAzimuth = ILOC_NULLVAL;
    Assocs[185].Slowness = ILOC_NULLVAL;
    Assocs[185].Timedef = 1;
    Assocs[185].Azimdef = 0;
    Assocs[185].Slowdef = 0;
    Assocs[185].Deltim = ILOC_NULLVAL;
    Assocs[186].StaInd = 99;    // EKS2
    Assocs[186].arid = 187;
    strcpy(Assocs[186].PhaseHint, "P");
    Assocs[186].phaseFixed = 0;
    Assocs[186].ArrivalTime = 1262324789.93;
    Assocs[186].BackAzimuth = ILOC_NULLVAL;
    Assocs[186].Slowness = ILOC_NULLVAL;
    Assocs[186].Timedef = 1;
    Assocs[186].Azimdef = 0;
    Assocs[186].Slowdef = 0;
    Assocs[186].Deltim = ILOC_NULLVAL;
    Assocs[187].StaInd = 100;    // KK31
    Assocs[187].arid = 188;
    strcpy(Assocs[187].PhaseHint, "P");
    Assocs[187].phaseFixed = 0;
    Assocs[187].ArrivalTime = 1262324806.4;
    Assocs[187].BackAzimuth = ILOC_NULLVAL;
    Assocs[187].Slowness = ILOC_NULLVAL;
    Assocs[187].Timedef = 1;
    Assocs[187].Azimdef = 0;
    Assocs[187].Slowdef = 0;
    Assocs[187].Deltim = ILOC_NULLVAL;
    Assocs[188].StaInd = 101;    // KKAR
    Assocs[188].arid = 189;
    strcpy(Assocs[188].PhaseHint, "P");
    Assocs[188].phaseFixed = 0;
    Assocs[188].ArrivalTime = 1262324806.3;
    Assocs[188].BackAzimuth = ILOC_NULLVAL;
    Assocs[188].Slowness = ILOC_NULLVAL;
    Assocs[188].Timedef = 1;
    Assocs[188].Azimdef = 0;
    Assocs[188].Slowdef = 0;
    Assocs[188].Deltim = ILOC_NULLVAL;
    Assocs[189].StaInd = 101;    // KKAR
    Assocs[189].arid = 190;
    strcpy(Assocs[189].PhaseHint, "P");
    Assocs[189].phaseFixed = 0;
    Assocs[189].ArrivalTime = 1262324806.32;
    Assocs[189].BackAzimuth = ILOC_NULLVAL;
    Assocs[189].Slowness = ILOC_NULLVAL;
    Assocs[189].Timedef = 1;
    Assocs[189].Azimdef = 0;
    Assocs[189].Slowdef = 0;
    Assocs[189].Deltim = ILOC_NULLVAL;
    Assocs[190].StaInd = 102;    // TIXI
    Assocs[190].arid = 191;
    strcpy(Assocs[190].PhaseHint, "P");
    Assocs[190].phaseFixed = 0;
    Assocs[190].ArrivalTime = 1262324822.3;
    Assocs[190].BackAzimuth = ILOC_NULLVAL;
    Assocs[190].Slowness = ILOC_NULLVAL;
    Assocs[190].Timedef = 1;
    Assocs[190].Azimdef = 0;
    Assocs[190].Slowdef = 0;
    Assocs[190].Deltim = ILOC_NULLVAL;
    Assocs[191].StaInd = 103;    // BVA0
    Assocs[191].arid = 192;
    strcpy(Assocs[191].PhaseHint, "P");
    Assocs[191].phaseFixed = 0;
    Assocs[191].ArrivalTime = 1262324828;
    Assocs[191].BackAzimuth = ILOC_NULLVAL;
    Assocs[191].Slowness = ILOC_NULLVAL;
    Assocs[191].Timedef = 1;
    Assocs[191].Azimdef = 0;
    Assocs[191].Slowdef = 0;
    Assocs[191].Deltim = ILOC_NULLVAL;
    Assocs[192].StaInd = 104;    // BRVK
    Assocs[192].arid = 193;
    strcpy(Assocs[192].PhaseHint, "P");
    Assocs[192].phaseFixed = 0;
    Assocs[192].ArrivalTime = 1262324827.5;
    Assocs[192].BackAzimuth = ILOC_NULLVAL;
    Assocs[192].Slowness = ILOC_NULLVAL;
    Assocs[192].Timedef = 1;
    Assocs[192].Azimdef = 0;
    Assocs[192].Slowdef = 0;
    Assocs[192].Deltim = ILOC_NULLVAL;
    Assocs[193].StaInd = 105;    // ABKAR
    Assocs[193].arid = 194;
    strcpy(Assocs[193].PhaseHint, "P");
    Assocs[193].phaseFixed = 0;
    Assocs[193].ArrivalTime = 1262324863.77;
    Assocs[193].BackAzimuth = ILOC_NULLVAL;
    Assocs[193].Slowness = ILOC_NULLVAL;
    Assocs[193].Timedef = 1;
    Assocs[193].Azimdef = 0;
    Assocs[193].Slowdef = 0;
    Assocs[193].Deltim = ILOC_NULLVAL;
    Assocs[194].StaInd = 106;    // SVE
    Assocs[194].arid = 195;
    strcpy(Assocs[194].PhaseHint, "P");
    Assocs[194].phaseFixed = 0;
    Assocs[194].ArrivalTime = 1262324869.4;
    Assocs[194].BackAzimuth = ILOC_NULLVAL;
    Assocs[194].Slowness = ILOC_NULLVAL;
    Assocs[194].Timedef = 1;
    Assocs[194].Azimdef = 0;
    Assocs[194].Slowdef = 0;
    Assocs[194].Deltim = ILOC_NULLVAL;
    Assocs[195].StaInd = 107;    // AKTO
    Assocs[195].arid = 196;
    strcpy(Assocs[195].PhaseHint, "P");
    Assocs[195].phaseFixed = 0;
    Assocs[195].ArrivalTime = 1262324873.2;
    Assocs[195].BackAzimuth = ILOC_NULLVAL;
    Assocs[195].Slowness = ILOC_NULLVAL;
    Assocs[195].Timedef = 1;
    Assocs[195].Azimdef = 0;
    Assocs[195].Slowdef = 0;
    Assocs[195].Deltim = ILOC_NULLVAL;
    Assocs[196].StaInd = 107;    // AKTO
    Assocs[196].arid = 197;
    strcpy(Assocs[196].PhaseHint, "P");
    Assocs[196].phaseFixed = 0;
    Assocs[196].ArrivalTime = 1262324873.15;
    Assocs[196].BackAzimuth = ILOC_NULLVAL;
    Assocs[196].Slowness = ILOC_NULLVAL;
    Assocs[196].Timedef = 1;
    Assocs[196].Azimdef = 0;
    Assocs[196].Slowdef = 0;
    Assocs[196].Deltim = ILOC_NULLVAL;
    Assocs[197].StaInd = 107;    // AKTO
    Assocs[197].arid = 198;
    strcpy(Assocs[197].PhaseHint, "P");
    Assocs[197].phaseFixed = 0;
    Assocs[197].ArrivalTime = 1262324873.15;
    Assocs[197].BackAzimuth = 83.8;
    Assocs[197].Slowness = 5.80;
    Assocs[197].Timedef = 1;
    Assocs[197].Azimdef = 1;
    Assocs[197].Slowdef = 1;
    Assocs[197].Deltim = ILOC_NULLVAL;
    Assocs[198].StaInd = 108;    // ARU
    Assocs[198].arid = 199;
    strcpy(Assocs[198].PhaseHint, "P");
    Assocs[198].phaseFixed = 0;
    Assocs[198].ArrivalTime = 1262324877.5;
    Assocs[198].BackAzimuth = ILOC_NULLVAL;
    Assocs[198].Slowness = ILOC_NULLVAL;
    Assocs[198].Timedef = 1;
    Assocs[198].Azimdef = 0;
    Assocs[198].Slowdef = 0;
    Assocs[198].Deltim = ILOC_NULLVAL;
    Assocs[199].StaInd = 108;    // ARU
    Assocs[199].arid = 200;
    strcpy(Assocs[199].PhaseHint, "S");
    Assocs[199].phaseFixed = 0;
    Assocs[199].ArrivalTime = 1262325405.1;
    Assocs[199].BackAzimuth = ILOC_NULLVAL;
    Assocs[199].Slowness = ILOC_NULLVAL;
    Assocs[199].Timedef = 1;
    Assocs[199].Azimdef = 0;
    Assocs[199].Slowdef = 0;
    Assocs[199].Deltim = ILOC_NULLVAL;
    Assocs[200].StaInd = 109;    // OHAK
    Assocs[200].arid = 201;
    strcpy(Assocs[200].PhaseHint, "P");
    Assocs[200].phaseFixed = 0;
    Assocs[200].ArrivalTime = 1262324926.21;
    Assocs[200].BackAzimuth = ILOC_NULLVAL;
    Assocs[200].Slowness = ILOC_NULLVAL;
    Assocs[200].Timedef = 1;
    Assocs[200].Azimdef = 0;
    Assocs[200].Slowdef = 0;
    Assocs[200].Deltim = ILOC_NULLVAL;
    Assocs[201].StaInd = 110;    // ZEI
    Assocs[201].arid = 202;
    strcpy(Assocs[201].PhaseHint, "P");
    Assocs[201].phaseFixed = 0;
    Assocs[201].ArrivalTime = 1262324931.2;
    Assocs[201].BackAzimuth = ILOC_NULLVAL;
    Assocs[201].Slowness = ILOC_NULLVAL;
    Assocs[201].Timedef = 1;
    Assocs[201].Azimdef = 0;
    Assocs[201].Slowdef = 0;
    Assocs[201].Deltim = ILOC_NULLVAL;
    Assocs[202].StaInd = 111;    // COLD
    Assocs[202].arid = 203;
    strcpy(Assocs[202].PhaseHint, "P");
    Assocs[202].phaseFixed = 0;
    Assocs[202].ArrivalTime = 1262324935.24;
    Assocs[202].BackAzimuth = ILOC_NULLVAL;
    Assocs[202].Slowness = ILOC_NULLVAL;
    Assocs[202].Timedef = 1;
    Assocs[202].Azimdef = 0;
    Assocs[202].Slowdef = 0;
    Assocs[202].Deltim = ILOC_NULLVAL;
    Assocs[203].StaInd = 112;    // KBZ
    Assocs[203].arid = 204;
    strcpy(Assocs[203].PhaseHint, "P");
    Assocs[203].phaseFixed = 0;
    Assocs[203].ArrivalTime = 1262324935.3;
    Assocs[203].BackAzimuth = ILOC_NULLVAL;
    Assocs[203].Slowness = ILOC_NULLVAL;
    Assocs[203].Timedef = 1;
    Assocs[203].Azimdef = 0;
    Assocs[203].Slowdef = 0;
    Assocs[203].Deltim = ILOC_NULLVAL;
    Assocs[204].StaInd = 112;    // KBZ
    Assocs[204].arid = 205;
    strcpy(Assocs[204].PhaseHint, "P");
    Assocs[204].phaseFixed = 0;
    Assocs[204].ArrivalTime = 1262324935.25;
    Assocs[204].BackAzimuth = ILOC_NULLVAL;
    Assocs[204].Slowness = ILOC_NULLVAL;
    Assocs[204].Timedef = 1;
    Assocs[204].Azimdef = 0;
    Assocs[204].Slowdef = 0;
    Assocs[204].Deltim = ILOC_NULLVAL;
    Assocs[205].StaInd = 112;    // KBZ
    Assocs[205].arid = 206;
    strcpy(Assocs[205].PhaseHint, "P");
    Assocs[205].phaseFixed = 0;
    Assocs[205].ArrivalTime = 1262324935.25;
    Assocs[205].BackAzimuth = 77.9;
    Assocs[205].Slowness = 14.50;
    Assocs[205].Timedef = 1;
    Assocs[205].Azimdef = 1;
    Assocs[205].Slowdef = 1;
    Assocs[205].Deltim = ILOC_NULLVAL;
    Assocs[206].StaInd = 113;    // KIV
    Assocs[206].arid = 207;
    strcpy(Assocs[206].PhaseHint, "P");
    Assocs[206].phaseFixed = 0;
    Assocs[206].ArrivalTime = 1262324937.2;
    Assocs[206].BackAzimuth = ILOC_NULLVAL;
    Assocs[206].Slowness = ILOC_NULLVAL;
    Assocs[206].Timedef = 1;
    Assocs[206].Azimdef = 0;
    Assocs[206].Slowdef = 0;
    Assocs[206].Deltim = ILOC_NULLVAL;
    Assocs[207].StaInd = 113;    // KIV
    Assocs[207].arid = 208;
    strcpy(Assocs[207].PhaseHint, "SS");
    Assocs[207].phaseFixed = 0;
    Assocs[207].ArrivalTime = 1262325807.8;
    Assocs[207].BackAzimuth = ILOC_NULLVAL;
    Assocs[207].Slowness = ILOC_NULLVAL;
    Assocs[207].Timedef = 1;
    Assocs[207].Azimdef = 0;
    Assocs[207].Slowdef = 0;
    Assocs[207].Deltim = ILOC_NULLVAL;
    Assocs[208].StaInd = 114;    // MCK
    Assocs[208].arid = 209;
    strcpy(Assocs[208].PhaseHint, "P");
    Assocs[208].phaseFixed = 0;
    Assocs[208].ArrivalTime = 1262324937.6;
    Assocs[208].BackAzimuth = ILOC_NULLVAL;
    Assocs[208].Slowness = ILOC_NULLVAL;
    Assocs[208].Timedef = 1;
    Assocs[208].Azimdef = 0;
    Assocs[208].Slowdef = 0;
    Assocs[208].Deltim = ILOC_NULLVAL;
    Assocs[209].StaInd = 115;    // COLA
    Assocs[209].arid = 210;
    strcpy(Assocs[209].PhaseHint, "P");
    Assocs[209].phaseFixed = 0;
    Assocs[209].ArrivalTime = 1262324940.8;
    Assocs[209].BackAzimuth = ILOC_NULLVAL;
    Assocs[209].Slowness = ILOC_NULLVAL;
    Assocs[209].Timedef = 1;
    Assocs[209].Azimdef = 0;
    Assocs[209].Slowdef = 0;
    Assocs[209].Deltim = ILOC_NULLVAL;
    Assocs[210].StaInd = 115;    // COLA
    Assocs[210].arid = 211;
    strcpy(Assocs[210].PhaseHint, "P");
    Assocs[210].phaseFixed = 0;
    Assocs[210].ArrivalTime = 1262324940.82;
    Assocs[210].BackAzimuth = ILOC_NULLVAL;
    Assocs[210].Slowness = ILOC_NULLVAL;
    Assocs[210].Timedef = 1;
    Assocs[210].Azimdef = 0;
    Assocs[210].Slowdef = 0;
    Assocs[210].Deltim = ILOC_NULLVAL;
    Assocs[211].StaInd = 116;    // ILAR
    Assocs[211].arid = 212;
    strcpy(Assocs[211].PhaseHint, "P");
    Assocs[211].phaseFixed = 0;
    Assocs[211].ArrivalTime = 1262324942.5;
    Assocs[211].BackAzimuth = ILOC_NULLVAL;
    Assocs[211].Slowness = ILOC_NULLVAL;
    Assocs[211].Timedef = 1;
    Assocs[211].Azimdef = 0;
    Assocs[211].Slowdef = 0;
    Assocs[211].Deltim = ILOC_NULLVAL;
    Assocs[212].StaInd = 116;    // ILAR
    Assocs[212].arid = 213;
    strcpy(Assocs[212].PhaseHint, "P");
    Assocs[212].phaseFixed = 0;
    Assocs[212].ArrivalTime = 1262324942.5;
    Assocs[212].BackAzimuth = 246.2;
    Assocs[212].Slowness = 6.90;
    Assocs[212].Timedef = 1;
    Assocs[212].Azimdef = 1;
    Assocs[212].Slowdef = 1;
    Assocs[212].Deltim = ILOC_NULLVAL;
    Assocs[213].StaInd = 116;    // ILAR
    Assocs[213].arid = 214;
    strcpy(Assocs[213].PhaseHint, "LR");
    Assocs[213].phaseFixed = 0;
    Assocs[213].ArrivalTime = 1262326975.315;
    Assocs[213].BackAzimuth = 279.6;
    Assocs[213].Slowness = 35.80;
    Assocs[213].Timedef = 1;
    Assocs[213].Azimdef = 1;
    Assocs[213].Slowdef = 1;
    Assocs[213].Deltim = ILOC_NULLVAL;
    Assocs[214].StaInd = 117;    // DOT
    Assocs[214].arid = 215;
    strcpy(Assocs[214].PhaseHint, "P");
    Assocs[214].phaseFixed = 0;
    Assocs[214].ArrivalTime = 1262324950.56;
    Assocs[214].BackAzimuth = ILOC_NULLVAL;
    Assocs[214].Slowness = ILOC_NULLVAL;
    Assocs[214].Timedef = 1;
    Assocs[214].Azimdef = 0;
    Assocs[214].Slowdef = 0;
    Assocs[214].Deltim = ILOC_NULLVAL;
    Assocs[215].StaInd = 118;    // OBN
    Assocs[215].arid = 216;
    strcpy(Assocs[215].PhaseHint, "P");
    Assocs[215].phaseFixed = 0;
    Assocs[215].ArrivalTime = 1262324950.5;
    Assocs[215].BackAzimuth = ILOC_NULLVAL;
    Assocs[215].Slowness = ILOC_NULLVAL;
    Assocs[215].Timedef = 1;
    Assocs[215].Azimdef = 0;
    Assocs[215].Slowdef = 0;
    Assocs[215].Deltim = ILOC_NULLVAL;
    Assocs[216].StaInd = 118;    // OBN
    Assocs[216].arid = 217;
    strcpy(Assocs[216].PhaseHint, "PcP");
    Assocs[216].phaseFixed = 0;
    Assocs[216].ArrivalTime = 1262324962.8;
    Assocs[216].BackAzimuth = ILOC_NULLVAL;
    Assocs[216].Slowness = ILOC_NULLVAL;
    Assocs[216].Timedef = 1;
    Assocs[216].Azimdef = 0;
    Assocs[216].Slowdef = 0;
    Assocs[216].Deltim = ILOC_NULLVAL;
    Assocs[217].StaInd = 118;    // OBN
    Assocs[217].arid = 218;
    strcpy(Assocs[217].PhaseHint, "S");
    Assocs[217].phaseFixed = 0;
    Assocs[217].ArrivalTime = 1262325545.9;
    Assocs[217].BackAzimuth = ILOC_NULLVAL;
    Assocs[217].Slowness = ILOC_NULLVAL;
    Assocs[217].Timedef = 1;
    Assocs[217].Azimdef = 0;
    Assocs[217].Slowdef = 0;
    Assocs[217].Deltim = ILOC_NULLVAL;
    Assocs[218].StaInd = 119;    // EGAK
    Assocs[218].arid = 219;
    strcpy(Assocs[218].PhaseHint, "P");
    Assocs[218].phaseFixed = 0;
    Assocs[218].ArrivalTime = 1262324956.65;
    Assocs[218].BackAzimuth = ILOC_NULLVAL;
    Assocs[218].Slowness = ILOC_NULLVAL;
    Assocs[218].Timedef = 1;
    Assocs[218].Azimdef = 0;
    Assocs[218].Slowdef = 0;
    Assocs[218].Deltim = ILOC_NULLVAL;
    Assocs[219].StaInd = 120;    // JOF
    Assocs[219].arid = 220;
    strcpy(Assocs[219].PhaseHint, "P");
    Assocs[219].phaseFixed = 0;
    Assocs[219].ArrivalTime = 1262324960;
    Assocs[219].BackAzimuth = ILOC_NULLVAL;
    Assocs[219].Slowness = ILOC_NULLVAL;
    Assocs[219].Timedef = 1;
    Assocs[219].Azimdef = 0;
    Assocs[219].Slowdef = 0;
    Assocs[219].Deltim = ILOC_NULLVAL;
    Assocs[220].StaInd = 121;    // ARCES
    Assocs[220].arid = 221;
    strcpy(Assocs[220].PhaseHint, "P");
    Assocs[220].phaseFixed = 0;
    Assocs[220].ArrivalTime = 1262324964.63;
    Assocs[220].BackAzimuth = ILOC_NULLVAL;
    Assocs[220].Slowness = ILOC_NULLVAL;
    Assocs[220].Timedef = 1;
    Assocs[220].Azimdef = 0;
    Assocs[220].Slowdef = 0;
    Assocs[220].Deltim = ILOC_NULLVAL;
    Assocs[221].StaInd = 121;    // ARCES
    Assocs[221].arid = 222;
    strcpy(Assocs[221].PhaseHint, "P");
    Assocs[221].phaseFixed = 0;
    Assocs[221].ArrivalTime = 1262324964.626;
    Assocs[221].BackAzimuth = 65.6;
    Assocs[221].Slowness = 5.90;
    Assocs[221].Timedef = 1;
    Assocs[221].Azimdef = 1;
    Assocs[221].Slowdef = 1;
    Assocs[221].Deltim = ILOC_NULLVAL;
    Assocs[222].StaInd = 122;    // INK
    Assocs[222].arid = 223;
    strcpy(Assocs[222].PhaseHint, "P");
    Assocs[222].phaseFixed = 0;
    Assocs[222].ArrivalTime = 1262324968.17;
    Assocs[222].BackAzimuth = ILOC_NULLVAL;
    Assocs[222].Slowness = ILOC_NULLVAL;
    Assocs[222].Timedef = 1;
    Assocs[222].Azimdef = 0;
    Assocs[222].Slowdef = 0;
    Assocs[222].Deltim = ILOC_NULLVAL;
    Assocs[223].StaInd = 122;    // INK
    Assocs[223].arid = 224;
    strcpy(Assocs[223].PhaseHint, "P");
    Assocs[223].phaseFixed = 0;
    Assocs[223].ArrivalTime = 1262324968.025;
    Assocs[223].BackAzimuth = 291.7;
    Assocs[223].Slowness = 6.60;
    Assocs[223].Timedef = 1;
    Assocs[223].Azimdef = 1;
    Assocs[223].Slowdef = 1;
    Assocs[223].Deltim = ILOC_NULLVAL;
    Assocs[224].StaInd = 123;    // KAF
    Assocs[224].arid = 225;
    strcpy(Assocs[224].PhaseHint, "P");
    Assocs[224].phaseFixed = 0;
    Assocs[224].ArrivalTime = 1262324972.2;
    Assocs[224].BackAzimuth = ILOC_NULLVAL;
    Assocs[224].Slowness = ILOC_NULLVAL;
    Assocs[224].Timedef = 1;
    Assocs[224].Azimdef = 0;
    Assocs[224].Slowdef = 0;
    Assocs[224].Deltim = ILOC_NULLVAL;
    Assocs[225].StaInd = 123;    // KAF
    Assocs[225].arid = 226;
    strcpy(Assocs[225].PhaseHint, "P");
    Assocs[225].phaseFixed = 0;
    Assocs[225].ArrivalTime = 1262324972.2;
    Assocs[225].BackAzimuth = ILOC_NULLVAL;
    Assocs[225].Slowness = ILOC_NULLVAL;
    Assocs[225].Timedef = 1;
    Assocs[225].Azimdef = 0;
    Assocs[225].Slowdef = 0;
    Assocs[225].Deltim = ILOC_NULLVAL;
    Assocs[226].StaInd = 124;    // FINES
    Assocs[226].arid = 227;
    strcpy(Assocs[226].PhaseHint, "P");
    Assocs[226].phaseFixed = 0;
    Assocs[226].ArrivalTime = 1262324973.98;
    Assocs[226].BackAzimuth = ILOC_NULLVAL;
    Assocs[226].Slowness = ILOC_NULLVAL;
    Assocs[226].Timedef = 1;
    Assocs[226].Azimdef = 0;
    Assocs[226].Slowdef = 0;
    Assocs[226].Deltim = ILOC_NULLVAL;
    Assocs[227].StaInd = 124;    // FINES
    Assocs[227].arid = 228;
    strcpy(Assocs[227].PhaseHint, "P");
    Assocs[227].phaseFixed = 0;
    Assocs[227].ArrivalTime = 1262324973.975;
    Assocs[227].BackAzimuth = 62.4;
    Assocs[227].Slowness = 4.50;
    Assocs[227].Timedef = 1;
    Assocs[227].Azimdef = 1;
    Assocs[227].Slowdef = 1;
    Assocs[227].Deltim = ILOC_NULLVAL;
    Assocs[228].StaInd = 125;    // BRTR
    Assocs[228].arid = 229;
    strcpy(Assocs[228].PhaseHint, "P");
    Assocs[228].phaseFixed = 0;
    Assocs[228].ArrivalTime = 1262324975.6;
    Assocs[228].BackAzimuth = ILOC_NULLVAL;
    Assocs[228].Slowness = ILOC_NULLVAL;
    Assocs[228].Timedef = 1;
    Assocs[228].Azimdef = 0;
    Assocs[228].Slowdef = 0;
    Assocs[228].Deltim = ILOC_NULLVAL;
    Assocs[229].StaInd = 125;    // BRTR
    Assocs[229].arid = 230;
    strcpy(Assocs[229].PhaseHint, "P");
    Assocs[229].phaseFixed = 0;
    Assocs[229].ArrivalTime = 1262324975.6;
    Assocs[229].BackAzimuth = 95.7;
    Assocs[229].Slowness = 3.70;
    Assocs[229].Timedef = 1;
    Assocs[229].Azimdef = 1;
    Assocs[229].Slowdef = 1;
    Assocs[229].Deltim = ILOC_NULLVAL;
    Assocs[230].StaInd = 126;    // AKASG
    Assocs[230].arid = 231;
    strcpy(Assocs[230].PhaseHint, "P");
    Assocs[230].phaseFixed = 0;
    Assocs[230].ArrivalTime = 1262324978.2;
    Assocs[230].BackAzimuth = ILOC_NULLVAL;
    Assocs[230].Slowness = ILOC_NULLVAL;
    Assocs[230].Timedef = 1;
    Assocs[230].Azimdef = 0;
    Assocs[230].Slowdef = 0;
    Assocs[230].Deltim = ILOC_NULLVAL;
    Assocs[231].StaInd = 126;    // AKASG
    Assocs[231].arid = 232;
    strcpy(Assocs[231].PhaseHint, "P");
    Assocs[231].phaseFixed = 0;
    Assocs[231].ArrivalTime = 1262324978.17;
    Assocs[231].BackAzimuth = ILOC_NULLVAL;
    Assocs[231].Slowness = ILOC_NULLVAL;
    Assocs[231].Timedef = 1;
    Assocs[231].Azimdef = 0;
    Assocs[231].Slowdef = 0;
    Assocs[231].Deltim = ILOC_NULLVAL;
    Assocs[232].StaInd = 126;    // AKASG
    Assocs[232].arid = 233;
    strcpy(Assocs[232].PhaseHint, "P");
    Assocs[232].phaseFixed = 0;
    Assocs[232].ArrivalTime = 1262324978.175;
    Assocs[232].BackAzimuth = 65.4;
    Assocs[232].Slowness = 4.20;
    Assocs[232].Timedef = 1;
    Assocs[232].Azimdef = 1;
    Assocs[232].Slowdef = 1;
    Assocs[232].Deltim = ILOC_NULLVAL;
    Assocs[233].StaInd = 127;    // BUR08
    Assocs[233].arid = 234;
    strcpy(Assocs[233].PhaseHint, "P");
    Assocs[233].phaseFixed = 0;
    Assocs[233].ArrivalTime = 1262324996.4;
    Assocs[233].BackAzimuth = ILOC_NULLVAL;
    Assocs[233].Slowness = ILOC_NULLVAL;
    Assocs[233].Timedef = 1;
    Assocs[233].Azimdef = 0;
    Assocs[233].Slowdef = 0;
    Assocs[233].Deltim = ILOC_NULLVAL;
    Assocs[234].StaInd = 128;    // MLR
    Assocs[234].arid = 235;
    strcpy(Assocs[234].PhaseHint, "P");
    Assocs[234].phaseFixed = 0;
    Assocs[234].ArrivalTime = 1262324997.14;
    Assocs[234].BackAzimuth = ILOC_NULLVAL;
    Assocs[234].Slowness = ILOC_NULLVAL;
    Assocs[234].Timedef = 1;
    Assocs[234].Azimdef = 0;
    Assocs[234].Slowdef = 0;
    Assocs[234].Deltim = ILOC_NULLVAL;
    Assocs[235].StaInd = 128;    // MLR
    Assocs[235].arid = 236;
    strcpy(Assocs[235].PhaseHint, "P");
    Assocs[235].phaseFixed = 0;
    Assocs[235].ArrivalTime = 1262324997.144;
    Assocs[235].BackAzimuth = 74.6;
    Assocs[235].Slowness = 4.60;
    Assocs[235].Timedef = 1;
    Assocs[235].Azimdef = 1;
    Assocs[235].Slowdef = 1;
    Assocs[235].Deltim = ILOC_NULLVAL;
    Assocs[236].StaInd = 129;    // RES
    Assocs[236].arid = 237;
    strcpy(Assocs[236].PhaseHint, "P");
    Assocs[236].phaseFixed = 0;
    Assocs[236].ArrivalTime = 1262325002.57;
    Assocs[236].BackAzimuth = ILOC_NULLVAL;
    Assocs[236].Slowness = ILOC_NULLVAL;
    Assocs[236].Timedef = 1;
    Assocs[236].Azimdef = 0;
    Assocs[236].Slowdef = 0;
    Assocs[236].Deltim = ILOC_NULLVAL;
    Assocs[237].StaInd = 129;    // RES
    Assocs[237].arid = 238;
    strcpy(Assocs[237].PhaseHint, "P");
    Assocs[237].phaseFixed = 0;
    Assocs[237].ArrivalTime = 1262325002.85;
    Assocs[237].BackAzimuth = 305.8;
    Assocs[237].Slowness = 5.80;
    Assocs[237].Timedef = 1;
    Assocs[237].Azimdef = 1;
    Assocs[237].Slowdef = 1;
    Assocs[237].Deltim = ILOC_NULLVAL;
    Assocs[238].StaInd = 130;    // KOLS
    Assocs[238].arid = 239;
    strcpy(Assocs[238].PhaseHint, "P");
    Assocs[238].phaseFixed = 0;
    Assocs[238].ArrivalTime = 1262325000.4;
    Assocs[238].BackAzimuth = ILOC_NULLVAL;
    Assocs[238].Slowness = ILOC_NULLVAL;
    Assocs[238].Timedef = 1;
    Assocs[238].Azimdef = 0;
    Assocs[238].Slowdef = 0;
    Assocs[238].Deltim = ILOC_NULLVAL;
    Assocs[239].StaInd = 131;    // CRVS
    Assocs[239].arid = 240;
    strcpy(Assocs[239].PhaseHint, "P");
    Assocs[239].phaseFixed = 0;
    Assocs[239].ArrivalTime = 1262325005.8;
    Assocs[239].BackAzimuth = ILOC_NULLVAL;
    Assocs[239].Slowness = ILOC_NULLVAL;
    Assocs[239].Timedef = 1;
    Assocs[239].Azimdef = 0;
    Assocs[239].Slowdef = 0;
    Assocs[239].Deltim = ILOC_NULLVAL;
    Assocs[240].StaInd = 132;    // NOA
    Assocs[240].arid = 241;
    strcpy(Assocs[240].PhaseHint, "LR");
    Assocs[240].phaseFixed = 0;
    Assocs[240].ArrivalTime = 1262327421.242;
    Assocs[240].BackAzimuth = 310.0;
    Assocs[240].Slowness = 35.70;
    Assocs[240].Timedef = 1;
    Assocs[240].Azimdef = 1;
    Assocs[240].Slowdef = 1;
    Assocs[240].Deltim = ILOC_NULLVAL;
    Assocs[241].StaInd = 133;    // YKA
    Assocs[241].arid = 242;
    strcpy(Assocs[241].PhaseHint, "P");
    Assocs[241].phaseFixed = 0;
    Assocs[241].ArrivalTime = 1262325015.35;
    Assocs[241].BackAzimuth = ILOC_NULLVAL;
    Assocs[241].Slowness = ILOC_NULLVAL;
    Assocs[241].Timedef = 1;
    Assocs[241].Azimdef = 0;
    Assocs[241].Slowdef = 0;
    Assocs[241].Deltim = ILOC_NULLVAL;
    Assocs[242].StaInd = 133;    // YKA
    Assocs[242].arid = 243;
    strcpy(Assocs[242].PhaseHint, "P");
    Assocs[242].phaseFixed = 0;
    Assocs[242].ArrivalTime = 1262325015.345;
    Assocs[242].BackAzimuth = 301.2;
    Assocs[242].Slowness = 5.00;
    Assocs[242].Timedef = 1;
    Assocs[242].Azimdef = 1;
    Assocs[242].Slowdef = 1;
    Assocs[242].Deltim = ILOC_NULLVAL;
    Assocs[243].StaInd = 134;    // PHWY
    Assocs[243].arid = 244;
    strcpy(Assocs[243].PhaseHint, "Pdif");
    Assocs[243].phaseFixed = 0;
    Assocs[243].ArrivalTime = 1262325089.91;
    Assocs[243].BackAzimuth = ILOC_NULLVAL;
    Assocs[243].Slowness = ILOC_NULLVAL;
    Assocs[243].Timedef = 1;
    Assocs[243].Azimdef = 0;
    Assocs[243].Slowdef = 0;
    Assocs[243].Deltim = ILOC_NULLVAL;
    Assocs[244].StaInd = 135;    // TXAR
    Assocs[244].arid = 245;
    strcpy(Assocs[244].PhaseHint, "PKPdf");
    Assocs[244].phaseFixed = 0;
    Assocs[244].ArrivalTime = 1262325356.11;
    Assocs[244].BackAzimuth = ILOC_NULLVAL;
    Assocs[244].Slowness = ILOC_NULLVAL;
    Assocs[244].Timedef = 1;
    Assocs[244].Azimdef = 0;
    Assocs[244].Slowdef = 0;
    Assocs[244].Deltim = ILOC_NULLVAL;
    Assocs[245].StaInd = 135;    // TXAR
    Assocs[245].arid = 246;
    strcpy(Assocs[245].PhaseHint, "PKPdf");
    Assocs[245].phaseFixed = 0;
    Assocs[245].ArrivalTime = 1262325356.113;
    Assocs[245].BackAzimuth = 317.9;
    Assocs[245].Slowness = 2.60;
    Assocs[245].Timedef = 1;
    Assocs[245].Azimdef = 1;
    Assocs[245].Slowdef = 1;
    Assocs[245].Deltim = ILOC_NULLVAL;
    Assocs[246].StaInd = 136;    // TORD
    Assocs[246].arid = 247;
    strcpy(Assocs[246].PhaseHint, "PKPdf");
    Assocs[246].phaseFixed = 0;
    Assocs[246].ArrivalTime = 1262325359.65;
    Assocs[246].BackAzimuth = ILOC_NULLVAL;
    Assocs[246].Slowness = ILOC_NULLVAL;
    Assocs[246].Timedef = 1;
    Assocs[246].Azimdef = 0;
    Assocs[246].Slowdef = 0;
    Assocs[246].Deltim = ILOC_NULLVAL;
    Assocs[247].StaInd = 136;    // TORD
    Assocs[247].arid = 248;
    strcpy(Assocs[247].PhaseHint, "PKPdf");
    Assocs[247].phaseFixed = 0;
    Assocs[247].ArrivalTime = 1262325359.65;
    Assocs[247].BackAzimuth = 69.4;
    Assocs[247].Slowness = 2.10;
    Assocs[247].Timedef = 1;
    Assocs[247].Azimdef = 1;
    Assocs[247].Slowdef = 1;
    Assocs[247].Deltim = ILOC_NULLVAL;
    Assocs[248].StaInd = 137;    // PLCA
    Assocs[248].arid = 249;
    strcpy(Assocs[248].PhaseHint, "PKPbc");
    Assocs[248].phaseFixed = 0;
    Assocs[248].ArrivalTime = 1262325422.15;
    Assocs[248].BackAzimuth = ILOC_NULLVAL;
    Assocs[248].Slowness = ILOC_NULLVAL;
    Assocs[248].Timedef = 1;
    Assocs[248].Azimdef = 0;
    Assocs[248].Slowdef = 0;
    Assocs[248].Deltim = ILOC_NULLVAL;
    Assocs[249].StaInd = 137;    // PLCA
    Assocs[249].arid = 250;
    strcpy(Assocs[249].PhaseHint, "PKPbc");
    Assocs[249].phaseFixed = 0;
    Assocs[249].ArrivalTime = 1262325421.975;
    Assocs[249].BackAzimuth = 215.0;
    Assocs[249].Slowness = 5.80;
    Assocs[249].Timedef = 1;
    Assocs[249].Azimdef = 1;
    Assocs[249].Slowdef = 1;
    Assocs[249].Deltim = ILOC_NULLVAL;
    Assocs[250].StaInd = 138;    // SDV
    Assocs[250].arid = 251;
    strcpy(Assocs[250].PhaseHint, "PKPdf");
    Assocs[250].phaseFixed = 0;
    Assocs[250].ArrivalTime = 1262325421.18;
    Assocs[250].BackAzimuth = ILOC_NULLVAL;
    Assocs[250].Slowness = ILOC_NULLVAL;
    Assocs[250].Timedef = 1;
    Assocs[250].Azimdef = 0;
    Assocs[250].Slowdef = 0;
    Assocs[250].Deltim = ILOC_NULLVAL;
    Assocs[251].StaInd = 138;    // SDV
    Assocs[251].arid = 252;
    strcpy(Assocs[251].PhaseHint, "PKPbc");
    Assocs[251].phaseFixed = 0;
    Assocs[251].ArrivalTime = 1262325428.92;
    Assocs[251].BackAzimuth = ILOC_NULLVAL;
    Assocs[251].Slowness = ILOC_NULLVAL;
    Assocs[251].Timedef = 1;
    Assocs[251].Azimdef = 0;
    Assocs[251].Slowdef = 0;
    Assocs[251].Deltim = ILOC_NULLVAL;
    Assocs[252].StaInd = 139;    // LPAZ
    Assocs[252].arid = 253;
    strcpy(Assocs[252].PhaseHint, "PKPdf");
    Assocs[252].phaseFixed = 0;
    Assocs[252].ArrivalTime = 1262325439.15;
    Assocs[252].BackAzimuth = ILOC_NULLVAL;
    Assocs[252].Slowness = ILOC_NULLVAL;
    Assocs[252].Timedef = 1;
    Assocs[252].Azimdef = 0;
    Assocs[252].Slowdef = 0;
    Assocs[252].Deltim = ILOC_NULLVAL;
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

