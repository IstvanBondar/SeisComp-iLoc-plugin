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
 * sciLocTTimes
 *
 * Istvan Bondar
 * ibondar2014@gmail.com
 *
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
 */

#include "sciLocInterface.h"
static int sortPredictedTT(int numPhase, ILOC_TT *predictedTT);

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
    ILOC_TT *predictedTT = (ILOC_TT *)NULL;
    char phase[ILOC_PHALEN];
    double lat, lon, depth, staLat, staLon, staElev, Delta, Esaz;
    int i, numPhase;
/*
 *
 *  set defaults for iLocConfig for tests
 *
 */

/*
 *  directory of auxiliary data files
 */
    strcpy(iLocConfig.auxdir, "/Users/istvanbondar/sciLocGit/iLocAuxDir");
    iLocConfig.Verbose = 0;
/*
 *  Travel time predictions
 */
    strcpy(iLocConfig.TTmodel, "iasp91");
    iLocConfig.UseRSTT = 0;
    strcpy(iLocConfig.RSTTmodel, "/Users/istvanbondar/sciLocGit/iLocAuxDir/RSTTmodels/pdu202009Du.geotess");
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
 *
    iLocConfig.MinIterations = 4;
    iLocConfig.MaxIterations = 100;
    iLocConfig.MinNdefPhases = 4;
    iLocConfig.SigmaThreshold = 6.;
    iLocConfig.DoCorrelatedErrors = 1;
    iLocConfig.AllowDamping = 1;
    iLocConfig.DoNotRenamePhases = 0;
*
 *  depth resolution
 *
    iLocConfig.MaxLocalDistDeg = 0.2;
    iLocConfig.MinLocalStations = 1;
    iLocConfig.MaxSPDistDeg = 2.;
    iLocConfig.MinSPpairs = 3;
    iLocConfig.MinCorePhases = 3;
    iLocConfig.MinDepthPhases = 3;
    iLocConfig.MaxShallowDepthError = 30.;
    iLocConfig.MaxDeepDepthError = 60.;
*
 *  NA search parameters
 *
    iLocConfig.DoGridSearch = 0;
    iLocConfig.NAsearchRadius = 5.;
    iLocConfig.NAsearchDepth = 300.;
    iLocConfig.NAsearchOT = 30.;
    iLocConfig.NAlpNorm = 1.;
    iLocConfig.NAiterMax = 5;
    iLocConfig.NAinitialSample = 1000;
    iLocConfig.NAnextSample = 100;
    iLocConfig.NAcells = 25;
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
 *  Travel time predictions
 *
 */
    strcpy(phase, "");
    for(;;) {
/*
 *      Infinite loop for source
 */
        printf("Enter source coordinates (finish with negative depth)\n");
        printf("  Source depth (km): ");
        scanf("%le", &depth);
        if (depth < 0) break;
        printf("  Source latitude (deg): ");
        scanf("%le", &lat);
        printf("  Source longitude (deg): ");
        scanf("%le", &lon);
        printf("Source lat=%.3f lon=%.3f depth=%.2f\n", lat, lon, depth);
        for (;;) {
/*
 *          Infinite loop for receiver
 */
            printf("Enter receiver coordinates (finish with negative elevation)\n");
            printf("  Receiver elevation (m): ");
            scanf("%le", &staElev);
            if (staElev < 0) break;
            printf("  Receiver latitude (deg): ");
            scanf("%le", &staLat);
            printf("  Receiver longitude (deg): ");
            scanf("%le", &staLon);
            printf("Receiver lat=%.3f lon=%.3f elevation=%.2f\n", staLat, staLon, staElev);
            predictedTT = iLoc_TravelTimePredictions(&iLocConfig, ec,
                               DefaultDepth.Topo, &TTInfo, TTtables,
                               &LocalTTInfo, LocalTTtables, phase, lat, lon, depth,
                               staLat, staLon, staElev, &Delta, &Esaz, &numPhase);
            if (numPhase) {
                sortPredictedTT(numPhase, predictedTT);
                printf("delta = %6.2f  esaz = %5.1f depth = %6.2f  n = %d\n",
                        Delta, Esaz, depth, numPhase);
                printf("  # phase    TTime[s] dTdD[s/deg] dTdh[s/km] vmodel\n");
                for (i = 0; i < numPhase; i++) {
                    printf("%3d %-6s %10.3f %11.6f %10.6f %s\n",
                          i+1, predictedTT[i].Phase, predictedTT[i].ttime,
                          predictedTT[i].dtdd, predictedTT[i].dtdh,
                          predictedTT[i].Vmodel);
                }
            }
            else {
                printf("no valid travel time has been found!\n");
            }
            iLoc_Free(predictedTT);
        }
    }
/*
 *
 *  Free allocated memory
 *
 */
    iLoc_FreeAuxData(&PhaseIdInfo, &fe, &DefaultDepth, &variogram, &TTInfo,
            TTtables, ec, &LocalTTInfo, LocalTTtables, iLocConfig.UseRSTT);
    return ILOC_SUCCESS;
}

/*
 *
 * SortPredictedTT: sorts predictedTT records by travel time
 *
 */
static int sortPredictedTT(int numPhase, ILOC_TT *predictedTT)
{
    int i, j;
    ILOC_TT temp;
    for (i = 1; i < numPhase; i++) {
        for (j = i - 1; j > -1; j--) {
            if (predictedTT[j].ttime > predictedTT[j+1].ttime) {
                ILOC_SWAP(predictedTT[j], predictedTT[j+1]);
            }
        }
    }
    return 0;
}
