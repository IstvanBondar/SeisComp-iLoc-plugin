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
    strcpy(iLocConfig.TTmodel, "iasp91");
    strcpy(iLocConfig.RSTTmodel, "/Users/istvanbondar/iLoc4.0/auxdata/RSTTmodels/pdu202009Du.geotess");
    iLocConfig.UseRSTTPnSn = 1;
    iLocConfig.UseRSTTPgLg = 1;
    strcpy(iLocConfig.LocalVmodel, "");
    iLocConfig.UseLocalTT = 0;
    iLocConfig.MaxLocalTTDelta = 3.;
    iLocConfig.UseRSTT = 1;
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
    iLocConfig.NAinitialSample = 1000;
    iLocConfig.NAnextSample = 100;
    iLocConfig.NAcells = 25;
/*
 *
 * set Hypocenter structure for test purposes
 *    slow converging Ukrainian ammunition explosion
 *
 */
    Hypocenter.isManMade = 0;
    Hypocenter.numSta = 29;
    Hypocenter.numPhase = 52;
    Hypocenter.Time = 1164129172.22;
    Hypocenter.Lat = 43.5334;
    Hypocenter.Lon = 45.9757;
    Hypocenter.Depth = 129.8;
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
    StaLocs[0].StaLat = 43.72330; StaLocs[0].StaLon = 44.73210; StaLocs[0].StaElevation = 141.0;
    StaLocs[1].StaLat = 43.07000; StaLocs[1].StaLon = 46.63000; StaLocs[1].StaElevation = 420.0;
    StaLocs[2].StaLat = 43.02000; StaLocs[2].StaLon = 46.83000; StaLocs[2].StaElevation = 900.0;
    StaLocs[3].StaLat = 43.06850; StaLocs[3].StaLon = 44.81170; StaLocs[3].StaElevation = 671.0;
    StaLocs[4].StaLat = 42.66389; StaLocs[4].StaLon = 46.22222; StaLocs[4].StaElevation = 870.0;
    StaLocs[5].StaLat = 43.04650; StaLocs[5].StaLon = 44.67730; StaLocs[5].StaElevation = 684.0;
    StaLocs[6].StaLat = 42.82694; StaLocs[6].StaLon = 46.90694; StaLocs[6].StaElevation = 1150.0;
    StaLocs[7].StaLat = 43.75150; StaLocs[7].StaLon = 44.28200; StaLocs[7].StaElevation = 136.0;
    StaLocs[8].StaLat = 42.71389; StaLocs[8].StaLon = 46.79444; StaLocs[8].StaElevation = 650.0;
    StaLocs[9].StaLat = 42.82500; StaLocs[9].StaLon = 47.10833; StaLocs[9].StaElevation = 480.0;
    StaLocs[10].StaLat = 44.80000; StaLocs[10].StaLon = 37.43300; StaLocs[10].StaElevation = 35.0;
    StaLocs[11].StaLat = 42.96100; StaLocs[11].StaLon = 47.50500; StaLocs[11].StaElevation = 42.0;
    StaLocs[12].StaLat = 42.82700; StaLocs[12].StaLon = 44.29700; StaLocs[12].StaElevation = 1271.0;
    StaLocs[13].StaLat = 43.08610; StaLocs[13].StaLon = 44.06770; StaLocs[13].StaElevation = 621.0;
    StaLocs[14].StaLat = 42.78810; StaLocs[14].StaLon = 43.90140; StaLocs[14].StaElevation = 1926.0;
    StaLocs[15].StaLat = 43.80000; StaLocs[15].StaLon = 43.40972; StaLocs[15].StaElevation = 665.0;
    StaLocs[16].StaLat = 42.09025; StaLocs[16].StaLon = 44.70189; StaLocs[16].StaElevation = 946.0;
    StaLocs[17].StaLat = 42.89940; StaLocs[17].StaLon = 43.58070; StaLocs[17].StaElevation = 1907.0;
    StaLocs[18].StaLat = 44.03333; StaLocs[18].StaLon = 43.05833; StaLocs[18].StaElevation = 544.0;
    StaLocs[19].StaLat = 42.59030; StaLocs[19].StaLon = 43.45250; StaLocs[19].StaElevation = 836.0;
    StaLocs[20].StaLat = 41.69392; StaLocs[20].StaLon = 44.79253; StaLocs[20].StaElevation = 515.0;
    StaLocs[21].StaLat = 41.45072; StaLocs[21].StaLon = 45.37317; StaLocs[21].StaElevation = 690.0;
    StaLocs[22].StaLat = 43.95531; StaLocs[22].StaLon = 42.68631; StaLocs[22].StaElevation = 1054.0;
    StaLocs[23].StaLat = 41.38150; StaLocs[23].StaLon = 44.41530; StaLocs[23].StaElevation = 740.0;
    StaLocs[24].StaLat = 40.14950; StaLocs[24].StaLon = 44.74139; StaLocs[24].StaElevation = 1583.0;
    StaLocs[25].StaLat = 44.80000; StaLocs[25].StaLon = 37.43300; StaLocs[25].StaElevation = 35.0;
    StaLocs[26].StaLat = 50.43480; StaLocs[26].StaLon = 58.01640; StaLocs[26].StaElevation = 379.0;
    StaLocs[27].StaLat = 49.25560; StaLocs[27].StaLon = 59.94310; StaLocs[27].StaElevation = 243.0;
    StaLocs[28].StaLat = 55.11380; StaLocs[28].StaLon = 36.56870; StaLocs[28].StaElevation = 0.0;
/*
 *  set Assocs structure
 */
    Assocs[0].StaInd = 0;   // TRKR
    Assocs[0].arid = 0;
    strcpy(Assocs[0].PhaseHint, "P");
    Assocs[0].phaseFixed = 0;
    Assocs[0].ArrivalTime = 1164107595.50;
    Assocs[0].BackAzimuth = ILOC_NULLVAL;
    Assocs[0].Slowness = ILOC_NULLVAL;
    Assocs[0].Timedef = 1;
    Assocs[0].Azimdef = 0;
    Assocs[0].Slowdef = 0;
    Assocs[0].Deltim = ILOC_NULLVAL;

    Assocs[1].StaInd = 0;
    Assocs[1].arid = 1;
    strcpy(Assocs[1].PhaseHint, "S");
    Assocs[1].phaseFixed = 0;
    Assocs[1].ArrivalTime = 1164107612.7;
    Assocs[1].BackAzimuth = ILOC_NULLVAL;
    Assocs[1].Slowness = ILOC_NULLVAL;
    Assocs[1].Timedef = 1;
    Assocs[1].Azimdef = 0;
    Assocs[1].Slowdef = 0;
    Assocs[1].Deltim = ILOC_NULLVAL;

    Assocs[2].StaInd = 1;   // DLMR
    Assocs[2].arid = 2;
    strcpy(Assocs[2].PhaseHint, "P");
    Assocs[2].phaseFixed = 0;
    Assocs[2].ArrivalTime = 1164129191.5;
    Assocs[2].BackAzimuth = ILOC_NULLVAL;
    Assocs[2].Slowness = ILOC_NULLVAL;
    Assocs[2].Timedef = 1;
    Assocs[2].Azimdef = 0;
    Assocs[2].Slowdef = 0;
    Assocs[2].Deltim = ILOC_NULLVAL;

    Assocs[3].StaInd = 1;
    Assocs[3].arid = 3;
    strcpy(Assocs[3].PhaseHint, "S");
    Assocs[3].phaseFixed = 1;
    Assocs[3].ArrivalTime = 1164129207.0;
    Assocs[3].BackAzimuth = ILOC_NULLVAL;
    Assocs[3].Slowness = ILOC_NULLVAL;
    Assocs[3].Timedef = 1;
    Assocs[3].Azimdef = 0;
    Assocs[3].Slowdef = 0;
    Assocs[3].Deltim = ILOC_NULLVAL;

    Assocs[4].StaInd = 2;   // DBC
    Assocs[4].arid = 4;
    strcpy(Assocs[4].PhaseHint, "P");
    Assocs[4].phaseFixed = 0;
    Assocs[4].ArrivalTime = 1164129193.7;
    Assocs[4].BackAzimuth = ILOC_NULLVAL;
    Assocs[4].Slowness = ILOC_NULLVAL;
    Assocs[4].Timedef = 1;
    Assocs[4].Azimdef = 0;
    Assocs[4].Slowdef = 0;
    Assocs[4].Deltim = ILOC_NULLVAL;

    Assocs[5].StaInd = 2;
    Assocs[5].arid = 5;
    strcpy(Assocs[5].PhaseHint, "S");
    Assocs[5].phaseFixed = 0;
    Assocs[5].ArrivalTime = 1164129210;
    Assocs[5].BackAzimuth = ILOC_NULLVAL;
    Assocs[5].Slowness = ILOC_NULLVAL;
    Assocs[5].Timedef = 1;
    Assocs[5].Azimdef = 0;
    Assocs[5].Slowdef = 0;
    Assocs[5].Deltim = ILOC_NULLVAL;

    Assocs[6].StaInd = 3;   // SNJR
    Assocs[6].arid = 6;
    strcpy(Assocs[6].PhaseHint, "P");
    Assocs[6].phaseFixed = 0;
    Assocs[6].ArrivalTime = 1164129195.6;
    Assocs[6].BackAzimuth = ILOC_NULLVAL;
    Assocs[6].Slowness = ILOC_NULLVAL;
    Assocs[6].Timedef = 1;
    Assocs[6].Azimdef = 0;
    Assocs[6].Slowdef = 0;
    Assocs[6].Deltim = ILOC_NULLVAL;

    Assocs[7].StaInd = 3;
    Assocs[7].arid = 7;
    strcpy(Assocs[7].PhaseHint, "S");
    Assocs[7].phaseFixed = 0;
    Assocs[7].ArrivalTime = 1164129213;
    Assocs[7].BackAzimuth = ILOC_NULLVAL;
    Assocs[7].Slowness = ILOC_NULLVAL;
    Assocs[7].Timedef = 1;
    Assocs[7].Azimdef = 0;
    Assocs[7].Slowdef = 0;
    Assocs[7].Deltim = ILOC_NULLVAL;

    Assocs[8].StaInd = 4;   // BTLR
    Assocs[8].arid = 8;
    strcpy(Assocs[8].PhaseHint, "P");
    Assocs[8].phaseFixed = 0;
    Assocs[8].ArrivalTime = 1164129193.5;
    Assocs[8].BackAzimuth = ILOC_NULLVAL;
    Assocs[8].Slowness = ILOC_NULLVAL;
    Assocs[8].Timedef = 1;
    Assocs[8].Azimdef = 0;
    Assocs[8].Slowdef = 0;
    Assocs[8].Deltim = ILOC_NULLVAL;

    Assocs[9].StaInd = 4;
    Assocs[9].arid = 9;
    strcpy(Assocs[9].PhaseHint, "S");
    Assocs[9].phaseFixed = 0;
    Assocs[9].ArrivalTime = 1164129209.5;
    Assocs[9].BackAzimuth = ILOC_NULLVAL;
    Assocs[9].Slowness = ILOC_NULLVAL;
    Assocs[9].Timedef = 1;
    Assocs[9].Azimdef = 0;
    Assocs[9].Slowdef = 0;
    Assocs[9].Deltim = ILOC_NULLVAL;

    Assocs[10].StaInd = 5;   // VLKR
    Assocs[10].arid = 10;
    strcpy(Assocs[10].PhaseHint, "P");
    Assocs[10].phaseFixed = 0;
    Assocs[10].ArrivalTime = 1164129196.3;
    Assocs[10].BackAzimuth = ILOC_NULLVAL;
    Assocs[10].Slowness = ILOC_NULLVAL;
    Assocs[10].Timedef = 1;
    Assocs[10].Azimdef = 0;
    Assocs[10].Slowdef = 0;
    Assocs[10].Deltim = ILOC_NULLVAL;

    Assocs[11].StaInd = 5;
    Assocs[11].arid = 11;
    strcpy(Assocs[11].PhaseHint, "S");
    Assocs[11].phaseFixed = 0;
    Assocs[11].ArrivalTime = 1164129214.3;
    Assocs[11].BackAzimuth = ILOC_NULLVAL;
    Assocs[11].Slowness = ILOC_NULLVAL;
    Assocs[11].Timedef = 1;
    Assocs[11].Azimdef = 0;
    Assocs[11].Slowdef = 0;
    Assocs[11].Deltim = ILOC_NULLVAL;

    Assocs[12].StaInd = 6;   // KRNR
    Assocs[12].arid = 12;
    strcpy(Assocs[12].PhaseHint, "P");
    Assocs[12].phaseFixed = 0;
    Assocs[12].ArrivalTime = 1164129194.5;
    Assocs[12].BackAzimuth = ILOC_NULLVAL;
    Assocs[12].Slowness = ILOC_NULLVAL;
    Assocs[12].Timedef = 1;
    Assocs[12].Azimdef = 0;
    Assocs[12].Slowdef = 0;
    Assocs[12].Deltim = ILOC_NULLVAL;

    Assocs[13].StaInd = 6;
    Assocs[13].arid = 13;
    strcpy(Assocs[13].PhaseHint, "S");
    Assocs[13].phaseFixed = 0;
    Assocs[13].ArrivalTime = 1164129212;
    Assocs[13].BackAzimuth = ILOC_NULLVAL;
    Assocs[13].Slowness = ILOC_NULLVAL;
    Assocs[13].Timedef = 1;
    Assocs[13].Azimdef = 0;
    Assocs[13].Slowdef = 0;
    Assocs[13].Deltim = ILOC_NULLVAL;

    Assocs[14].StaInd = 7;   // PRTR
    Assocs[14].arid = 14;
    strcpy(Assocs[14].PhaseHint, "P");
    Assocs[14].phaseFixed = 0;
    Assocs[14].ArrivalTime = 1164129198.2;
    Assocs[14].BackAzimuth = ILOC_NULLVAL;
    Assocs[14].Slowness = ILOC_NULLVAL;
    Assocs[14].Timedef = 1;
    Assocs[14].Azimdef = 0;
    Assocs[14].Slowdef = 0;
    Assocs[14].Deltim = ILOC_NULLVAL;

    Assocs[15].StaInd = 7;
    Assocs[15].arid = 15;
    strcpy(Assocs[15].PhaseHint, "S");
    Assocs[15].phaseFixed = 0;
    Assocs[15].ArrivalTime = 1164129217.5;
    Assocs[15].BackAzimuth = ILOC_NULLVAL;
    Assocs[15].Slowness = ILOC_NULLVAL;
    Assocs[15].Timedef = 1;
    Assocs[15].Azimdef = 0;
    Assocs[15].Slowdef = 0;
    Assocs[15].Deltim = ILOC_NULLVAL;

    Assocs[16].StaInd = 8;   // UNCR
    Assocs[16].arid = 16;
    strcpy(Assocs[16].PhaseHint, "P");
    Assocs[16].phaseFixed = 0;
    Assocs[16].ArrivalTime = 1164129194.5;
    Assocs[16].BackAzimuth = ILOC_NULLVAL;
    Assocs[16].Slowness = ILOC_NULLVAL;
    Assocs[16].Timedef = 1;
    Assocs[16].Azimdef = 0;
    Assocs[16].Slowdef = 0;
    Assocs[16].Deltim = ILOC_NULLVAL;

    Assocs[17].StaInd = 8;
    Assocs[17].arid = 17;
    strcpy(Assocs[17].PhaseHint, "S");
    Assocs[17].phaseFixed = 0;
    Assocs[17].ArrivalTime = 1164129212;
    Assocs[17].BackAzimuth = ILOC_NULLVAL;
    Assocs[17].Slowness = ILOC_NULLVAL;
    Assocs[17].Timedef = 1;
    Assocs[17].Azimdef = 0;
    Assocs[17].Slowdef = 0;
    Assocs[17].Deltim = ILOC_NULLVAL;

    Assocs[18].StaInd = 9;   // BURJ
    Assocs[18].arid = 18;
    strcpy(Assocs[18].PhaseHint, "P");
    Assocs[18].phaseFixed = 0;
    Assocs[18].ArrivalTime = 1164129196;
    Assocs[18].BackAzimuth = ILOC_NULLVAL;
    Assocs[18].Slowness = ILOC_NULLVAL;
    Assocs[18].Timedef = 1;
    Assocs[18].Azimdef = 0;
    Assocs[18].Slowdef = 0;
    Assocs[18].Deltim = ILOC_NULLVAL;

    Assocs[19].StaInd = 9;
    Assocs[19].arid = 19;
    strcpy(Assocs[19].PhaseHint, "S");
    Assocs[19].phaseFixed = 0;
    Assocs[19].ArrivalTime = 1164129215;
    Assocs[19].BackAzimuth = ILOC_NULLVAL;
    Assocs[19].Slowness = ILOC_NULLVAL;
    Assocs[19].Timedef = 1;
    Assocs[19].Azimdef = 0;
    Assocs[19].Slowdef = 0;
    Assocs[19].Deltim = ILOC_NULLVAL;

    Assocs[20].StaInd = 10;   // ARNR
    Assocs[20].arid = 20;
    strcpy(Assocs[20].PhaseHint, "P");
    Assocs[20].phaseFixed = 0;
    Assocs[20].ArrivalTime = 1164129198.5;
    Assocs[20].BackAzimuth = ILOC_NULLVAL;
    Assocs[20].Slowness = ILOC_NULLVAL;
    Assocs[20].Timedef = 1;
    Assocs[20].Azimdef = 0;
    Assocs[20].Slowdef = 0;
    Assocs[20].Deltim = ILOC_NULLVAL;

    Assocs[21].StaInd = 10;
    Assocs[21].arid = 21;
    strcpy(Assocs[21].PhaseHint, "S");
    Assocs[21].phaseFixed = 0;
    Assocs[21].ArrivalTime = 1164129215;
    Assocs[21].BackAzimuth = ILOC_NULLVAL;
    Assocs[21].Slowness = ILOC_NULLVAL;
    Assocs[21].Timedef = 1;
    Assocs[21].Azimdef = 0;
    Assocs[21].Slowdef = 0;
    Assocs[21].Deltim = ILOC_NULLVAL;

    Assocs[22].StaInd = 11;   // MAK
    Assocs[22].arid = 22;
    strcpy(Assocs[22].PhaseHint, "P");
    Assocs[22].phaseFixed = 0;
    Assocs[22].ArrivalTime = 1164129201.2;
    Assocs[22].BackAzimuth = ILOC_NULLVAL;
    Assocs[22].Slowness = ILOC_NULLVAL;
    Assocs[22].Timedef = 1;
    Assocs[22].Azimdef = 0;
    Assocs[22].Slowdef = 0;
    Assocs[22].Deltim = ILOC_NULLVAL;

    Assocs[23].StaInd = 11;
    Assocs[23].arid = 23;
    strcpy(Assocs[23].PhaseHint, "S");
    Assocs[23].phaseFixed = 0;
    Assocs[23].ArrivalTime = 1164129221.9;
    Assocs[23].BackAzimuth = ILOC_NULLVAL;
    Assocs[23].Slowness = ILOC_NULLVAL;
    Assocs[23].Timedef = 1;
    Assocs[23].Azimdef = 0;
    Assocs[23].Slowdef = 0;
    Assocs[23].Deltim = ILOC_NULLVAL;

    Assocs[24].StaInd = 12;   // LACR
    Assocs[24].arid = 24;
    strcpy(Assocs[24].PhaseHint, "P");
    Assocs[24].phaseFixed = 0;
    Assocs[24].ArrivalTime = 1164129198.9;
    Assocs[24].BackAzimuth = ILOC_NULLVAL;
    Assocs[24].Slowness = ILOC_NULLVAL;
    Assocs[24].Timedef = 1;
    Assocs[24].Azimdef = 0;
    Assocs[24].Slowdef = 0;
    Assocs[24].Deltim = ILOC_NULLVAL;

    Assocs[25].StaInd = 12;
    Assocs[25].arid = 25;
    strcpy(Assocs[25].PhaseHint, "S");
    Assocs[25].phaseFixed = 0;
    Assocs[25].ArrivalTime = 1164129217.8;
    Assocs[25].BackAzimuth = ILOC_NULLVAL;
    Assocs[25].Slowness = ILOC_NULLVAL;
    Assocs[25].Timedef = 1;
    Assocs[25].Azimdef = 0;
    Assocs[25].Slowdef = 0;
    Assocs[25].Deltim = ILOC_NULLVAL;

    Assocs[26].StaInd = 13;   // KORR
    Assocs[26].arid = 26;
    strcpy(Assocs[26].PhaseHint, "P");
    Assocs[26].phaseFixed = 0;
    Assocs[26].ArrivalTime = 1164129199.8;
    Assocs[26].BackAzimuth = ILOC_NULLVAL;
    Assocs[26].Slowness = ILOC_NULLVAL;
    Assocs[26].Timedef = 1;
    Assocs[26].Azimdef = 0;
    Assocs[26].Slowdef = 0;
    Assocs[26].Deltim = ILOC_NULLVAL;

    Assocs[27].StaInd = 13;
    Assocs[27].arid = 27;
    strcpy(Assocs[27].PhaseHint, "S");
    Assocs[27].phaseFixed = 0;
    Assocs[27].ArrivalTime = 1164129220.3;
    Assocs[27].BackAzimuth = ILOC_NULLVAL;
    Assocs[27].Slowness = ILOC_NULLVAL;
    Assocs[27].Timedef = 1;
    Assocs[27].Azimdef = 0;
    Assocs[27].Slowdef = 0;
    Assocs[27].Deltim = ILOC_NULLVAL;

    Assocs[28].StaInd = 14;   // ZEI
    Assocs[28].arid = 28;
    strcpy(Assocs[28].PhaseHint, "P");
    Assocs[28].phaseFixed = 0;
    Assocs[28].ArrivalTime = 1164129201.8;
    Assocs[28].BackAzimuth = ILOC_NULLVAL;
    Assocs[28].Slowness = ILOC_NULLVAL;
    Assocs[28].Timedef = 1;
    Assocs[28].Azimdef = 0;
    Assocs[28].Slowdef = 0;
    Assocs[28].Deltim = ILOC_NULLVAL;

    Assocs[29].StaInd = 14;
    Assocs[29].arid = 29;
    strcpy(Assocs[29].PhaseHint, "S");
    Assocs[29].phaseFixed = 0;
    Assocs[29].ArrivalTime = 1164129222.7;
    Assocs[29].BackAzimuth = ILOC_NULLVAL;
    Assocs[29].Slowness = ILOC_NULLVAL;
    Assocs[29].Timedef = 1;
    Assocs[29].Azimdef = 0;
    Assocs[29].Slowdef = 0;
    Assocs[29].Deltim = ILOC_NULLVAL;

    Assocs[30].StaInd = 15;   // KUBR
    Assocs[30].arid = 30;
    strcpy(Assocs[30].PhaseHint, "P");
    Assocs[30].phaseFixed = 0;
    Assocs[30].ArrivalTime = 1164129204.3;
    Assocs[30].BackAzimuth = ILOC_NULLVAL;
    Assocs[30].Slowness = ILOC_NULLVAL;
    Assocs[30].Timedef = 1;
    Assocs[30].Azimdef = 0;
    Assocs[30].Slowdef = 0;
    Assocs[30].Deltim = ILOC_NULLVAL;

    Assocs[31].StaInd = 16;   // DUS
    Assocs[31].arid = 31;
    strcpy(Assocs[31].PhaseHint, "P");
    Assocs[31].phaseFixed = 0;
    Assocs[31].ArrivalTime = 1164129203.781;
    Assocs[31].BackAzimuth = ILOC_NULLVAL;
    Assocs[31].Slowness = ILOC_NULLVAL;
    Assocs[31].Timedef = 1;
    Assocs[31].Azimdef = 0;
    Assocs[31].Slowdef = 0;
    Assocs[31].Deltim = ILOC_NULLVAL;

    Assocs[32].StaInd = 16;
    Assocs[32].arid = 32;
    strcpy(Assocs[32].PhaseHint, "S");
    Assocs[32].phaseFixed = 0;
    Assocs[32].ArrivalTime = 1164129227.2;
    Assocs[32].BackAzimuth = ILOC_NULLVAL;
    Assocs[32].Slowness = ILOC_NULLVAL;
    Assocs[32].Timedef = 1;
    Assocs[32].Azimdef = 0;
    Assocs[32].Slowdef = 0;
    Assocs[32].Deltim = ILOC_NULLVAL;

    Assocs[33].StaInd = 17;   // DIGR
    Assocs[33].arid = 33;
    strcpy(Assocs[33].PhaseHint, "P");
    Assocs[33].phaseFixed = 0;
    Assocs[33].ArrivalTime = 1164129203.7;
    Assocs[33].BackAzimuth = ILOC_NULLVAL;
    Assocs[33].Slowness = ILOC_NULLVAL;
    Assocs[33].Timedef = 1;
    Assocs[33].Azimdef = 0;
    Assocs[33].Slowdef = 0;
    Assocs[33].Deltim = ILOC_NULLVAL;

    Assocs[34].StaInd = 17;
    Assocs[34].arid = 34;
    strcpy(Assocs[34].PhaseHint, "S");
    Assocs[34].phaseFixed = 0;
    Assocs[34].ArrivalTime = 1164129226.1;
    Assocs[34].BackAzimuth = ILOC_NULLVAL;
    Assocs[34].Slowness = ILOC_NULLVAL;
    Assocs[34].Timedef = 1;
    Assocs[34].Azimdef = 0;
    Assocs[34].Slowdef = 0;
    Assocs[34].Deltim = ILOC_NULLVAL;

    Assocs[35].StaInd = 18;   // PYA
    Assocs[35].arid = 35;
    strcpy(Assocs[35].PhaseHint, "P");
    Assocs[35].phaseFixed = 0;
    Assocs[35].ArrivalTime = 1164129207.1;
    Assocs[35].BackAzimuth = ILOC_NULLVAL;
    Assocs[35].Slowness = ILOC_NULLVAL;
    Assocs[35].Timedef = 1;
    Assocs[35].Azimdef = 0;
    Assocs[35].Slowdef = 0;
    Assocs[35].Deltim = ILOC_NULLVAL;

    Assocs[36].StaInd = 19;   // ONI
    Assocs[36].arid = 36;
    strcpy(Assocs[36].PhaseHint, "P");
    Assocs[36].phaseFixed = 0;
    Assocs[36].ArrivalTime = 1164129207;
    Assocs[36].BackAzimuth = ILOC_NULLVAL;
    Assocs[36].Slowness = ILOC_NULLVAL;
    Assocs[36].Timedef = 1;
    Assocs[36].Azimdef = 0;
    Assocs[36].Slowdef = 0;
    Assocs[36].Deltim = ILOC_NULLVAL;

    Assocs[37].StaInd = 19;
    Assocs[37].arid = 37;
    strcpy(Assocs[37].PhaseHint, "S");
    Assocs[37].phaseFixed = 0;
    Assocs[37].ArrivalTime = 1164129232.5;
    Assocs[37].BackAzimuth = ILOC_NULLVAL;
    Assocs[37].Slowness = ILOC_NULLVAL;
    Assocs[37].Timedef = 1;
    Assocs[37].Azimdef = 0;
    Assocs[37].Slowdef = 0;
    Assocs[37].Deltim = ILOC_NULLVAL;

    Assocs[38].StaInd = 20;   // MTA
    Assocs[38].arid = 38;
    strcpy(Assocs[38].PhaseHint, "P");
    Assocs[38].phaseFixed = 0;
    Assocs[38].ArrivalTime = 1164129216.2;
    Assocs[38].BackAzimuth = ILOC_NULLVAL;
    Assocs[38].Slowness = ILOC_NULLVAL;
    Assocs[38].Timedef = 1;
    Assocs[38].Azimdef = 0;
    Assocs[38].Slowdef = 0;
    Assocs[38].Deltim = ILOC_NULLVAL;

    Assocs[39].StaInd = 20;
    Assocs[39].arid = 39;
    strcpy(Assocs[39].PhaseHint, "S");
    Assocs[39].phaseFixed = 0;
    Assocs[39].ArrivalTime = 1164129241.6;
    Assocs[39].BackAzimuth = ILOC_NULLVAL;
    Assocs[39].Slowness = ILOC_NULLVAL;
    Assocs[39].Timedef = 1;
    Assocs[39].Azimdef = 0;
    Assocs[39].Slowdef = 0;
    Assocs[39].Deltim = ILOC_NULLVAL;

    Assocs[40].StaInd = 21;   // DGRG
    Assocs[40].arid = 40;
    strcpy(Assocs[40].PhaseHint, "P");
    Assocs[40].phaseFixed = 0;
    Assocs[40].ArrivalTime = 1164129209.1;
    Assocs[40].BackAzimuth = ILOC_NULLVAL;
    Assocs[40].Slowness = ILOC_NULLVAL;
    Assocs[40].Timedef = 1;
    Assocs[40].Azimdef = 0;
    Assocs[40].Slowdef = 0;
    Assocs[40].Deltim = ILOC_NULLVAL;

    Assocs[41].StaInd = 21;
    Assocs[41].arid = 41;
    strcpy(Assocs[41].PhaseHint, "S");
    Assocs[41].phaseFixed = 0;
    Assocs[41].ArrivalTime = 1164129235.8;
    Assocs[41].BackAzimuth = ILOC_NULLVAL;
    Assocs[41].Slowness = ILOC_NULLVAL;
    Assocs[41].Timedef = 1;
    Assocs[41].Azimdef = 0;
    Assocs[41].Slowdef = 0;
    Assocs[41].Deltim = ILOC_NULLVAL;

    Assocs[42].StaInd = 22;   // KIV
    Assocs[42].arid = 42;
    strcpy(Assocs[42].PhaseHint, "P");
    Assocs[42].phaseFixed = 0;
    Assocs[42].ArrivalTime = 1164129210.3;
    Assocs[42].BackAzimuth = ILOC_NULLVAL;
    Assocs[42].Slowness = ILOC_NULLVAL;
    Assocs[42].Timedef = 1;
    Assocs[42].Azimdef = 0;
    Assocs[42].Slowdef = 0;
    Assocs[42].Deltim = ILOC_NULLVAL;

    Assocs[43].StaInd = 23;   // KZR
    Assocs[43].arid = 43;
    strcpy(Assocs[43].PhaseHint, "P");
    Assocs[43].phaseFixed = 0;
    Assocs[43].ArrivalTime = 1164129211.586;
    Assocs[43].BackAzimuth = ILOC_NULLVAL;
    Assocs[43].Slowness = ILOC_NULLVAL;
    Assocs[43].Timedef = 1;
    Assocs[43].Azimdef = 0;
    Assocs[43].Slowdef = 0;
    Assocs[43].Deltim = ILOC_NULLVAL;

    Assocs[44].StaInd = 23;
    Assocs[44].arid = 44;
    strcpy(Assocs[44].PhaseHint, "S");
    Assocs[44].phaseFixed = 0;
    Assocs[44].ArrivalTime = 1164129241.886;
    Assocs[44].BackAzimuth = ILOC_NULLVAL;
    Assocs[44].Slowness = ILOC_NULLVAL;
    Assocs[44].Timedef = 1;
    Assocs[44].Azimdef = 0;
    Assocs[44].Slowdef = 0;
    Assocs[44].Deltim = ILOC_NULLVAL;

    Assocs[45].StaInd = 24;   // GNI
    Assocs[45].arid = 45;
    strcpy(Assocs[45].PhaseHint, "P");
    Assocs[45].phaseFixed = 0;
    Assocs[45].ArrivalTime = 1164129226.5;
    Assocs[45].BackAzimuth = ILOC_NULLVAL;
    Assocs[45].Slowness = ILOC_NULLVAL;
    Assocs[45].Timedef = 1;
    Assocs[45].Azimdef = 0;
    Assocs[45].Slowdef = 0;
    Assocs[45].Deltim = ILOC_NULLVAL;

    Assocs[46].StaInd = 25;   // ANN
    Assocs[46].arid = 46;
    strcpy(Assocs[46].PhaseHint, "P");
    Assocs[46].phaseFixed = 0;
    Assocs[46].ArrivalTime = 1164129264.3;
    Assocs[46].BackAzimuth = ILOC_NULLVAL;
    Assocs[46].Slowness = ILOC_NULLVAL;
    Assocs[46].Timedef = 1;
    Assocs[46].Azimdef = 0;
    Assocs[46].Slowdef = 0;
    Assocs[46].Deltim = ILOC_NULLVAL;

    Assocs[47].StaInd = 25;
    Assocs[47].arid = 47;
    strcpy(Assocs[47].PhaseHint, "S");
    Assocs[47].phaseFixed = 0;
    Assocs[47].ArrivalTime = 1164129333.4;
    Assocs[47].BackAzimuth = ILOC_NULLVAL;
    Assocs[47].Slowness = ILOC_NULLVAL;
    Assocs[47].Timedef = 1;
    Assocs[47].Azimdef = 0;
    Assocs[47].Slowdef = 0;
    Assocs[47].Deltim = ILOC_NULLVAL;

    Assocs[48].StaInd = 26;   // AKTO
    Assocs[48].arid = 48;
    strcpy(Assocs[48].PhaseHint, "P");
    Assocs[48].phaseFixed = 0;
    Assocs[48].ArrivalTime = 1164129321.484;
    Assocs[48].BackAzimuth = ILOC_NULLVAL;
    Assocs[48].Slowness = ILOC_NULLVAL;
    Assocs[48].Timedef = 1;
    Assocs[48].Azimdef = 0;
    Assocs[48].Slowdef = 0;
    Assocs[48].Deltim = ILOC_NULLVAL;

    Assocs[49].StaInd = 26;
    Assocs[49].arid = 49;
    strcpy(Assocs[49].PhaseHint, "S");
    Assocs[49].phaseFixed = 0;
    Assocs[49].ArrivalTime = 1164129434.831;
    Assocs[49].BackAzimuth = ILOC_NULLVAL;
    Assocs[49].Slowness = ILOC_NULLVAL;
    Assocs[49].Timedef = 1;
    Assocs[49].Azimdef = 0;
    Assocs[49].Slowdef = 0;
    Assocs[49].Deltim = ILOC_NULLVAL;

    Assocs[50].StaInd = 27;   // AB31
    Assocs[50].arid = 50;
    strcpy(Assocs[50].PhaseHint, "P");
    Assocs[50].phaseFixed = 0;
    Assocs[50].ArrivalTime = 1164129328.214;
    Assocs[50].BackAzimuth = ILOC_NULLVAL;
    Assocs[50].Slowness = ILOC_NULLVAL;
    Assocs[50].Timedef = 1;
    Assocs[50].Azimdef = 0;
    Assocs[50].Slowdef = 0;
    Assocs[50].Deltim = ILOC_NULLVAL;

    Assocs[51].StaInd = 27;
    Assocs[51].arid = 51;
    strcpy(Assocs[51].PhaseHint, "S");
    Assocs[51].phaseFixed = 0;
    Assocs[51].ArrivalTime = 1164129445.812;
    Assocs[51].BackAzimuth = ILOC_NULLVAL;
    Assocs[51].Slowness = ILOC_NULLVAL;
    Assocs[51].Timedef = 1;
    Assocs[51].Azimdef = 0;
    Assocs[51].Slowdef = 0;
    Assocs[51].Deltim = ILOC_NULLVAL;

    Assocs[52].StaInd = 28;   // OBN
    Assocs[52].arid = 52;
    strcpy(Assocs[52].PhaseHint, "P");
    Assocs[52].phaseFixed = 0;
    Assocs[52].ArrivalTime = 1164129328.214;
    Assocs[52].BackAzimuth = ILOC_NULLVAL;
    Assocs[52].Slowness = ILOC_NULLVAL;
    Assocs[52].Timedef = 1;
    Assocs[52].Azimdef = 0;
    Assocs[52].Slowdef = 0;
    Assocs[52].Deltim = ILOC_NULLVAL;

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
    retval = iLoc_Locator(&iLocConfig, &PhaseIdInfo, &fe, &DefaultDepth,
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

