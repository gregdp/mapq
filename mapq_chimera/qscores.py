

# Copyright (c) 2020 Greg Pintilie - gregp@slac.stanford.edu

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.



import numpy
import _multiscale
from CGLutil.AdaptiveTree import AdaptiveTree
import chimera
import FitMap
import os
import Matrix
import VolumeData



chargedIons = { "MG":2, "NA":1, "CL":-1, "CA":2, "ZN":2, "MN":2, "FE":3, "CO":2, "NI":2 }



# returns the min and max density value in a map

def MinMaxD ( dmap ) :

    # dmap - the map

    M = dmap.data.full_matrix()
    maxM = numpy.max(M)
    minM = numpy.min(M)

    maxD = min ( numpy.average(M)+numpy.std(M)*10, maxM )
    minD = max ( numpy.average(M)-numpy.std(M)*1, minM )

    # xray
    #maxD = min ( numpy.average(M)+numpy.std(M)*3.5, maxM )
    #minD = max ( numpy.average(M)-numpy.std(M)*0.77, minM )

    #print "%s - %.2f->%.2f, %.2f->%.2f" % (dmap.name, minD, maxD, minM, maxM )
    #minD = numpy.min(M)
    #minD, maxD = numpy.min(M), numpy.max(M)
    return minD, maxD





# this method calculates CC between radial points placed around the atoms and the map
# - two values are returned - basic CC and CC about the mean - the latter is the Q-scre


def Qscore ( atoms, dmap, sigma, allAtTree = None, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :


    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    #sigma = 1.0

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)
    #print " __%s__ " % (atoms[0].name), pts[0]


    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals should have the reference gaussian...
    g_vals = (numpy.ones ( [len(pts)*numPts,1] ) * refG).astype(numpy.float64, copy=False)
    g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False)

    if mol == None :
        mol = atoms[0].molecule


    # r_avg holds the average values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( pts, mol.openState.xform ).astype(numpy.float64, copy=False)
    d_vals = numpy.repeat ( d_vals, numPts )

    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,len(pts)*numPts] ]

    d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)


    # make smaller atom tree...
    if 1 and allAtTree != None :
        ats_near = []
        for at in atoms :
            anear = allAtTree.searchTree ( at.coord().data(), toRAD*2.0 )
            ats_near.extend ( anear )

        points = _multiscale.get_atom_coordinates ( ats_near, transformed = False )
        if log :
            print " - new search tree: %d pts" % ( len(ats_near) )
        allAtTree = AdaptiveTree ( points.tolist(), ats_near, 1.0)



    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        #outRad2 = outRad * outRad
        pts = []
        for at in atoms :
            #npts = numPts # 8 # int ( npts )
            npts = int (numPts * RAD*RAD / (dRAD*dRAD)) if show else numPts
            #npts = numPts * (RAD*RAD / (dRAD*dRAD))
            #print RAD, dRAD, numPts, " -> ", npts
            for i in range (0, 100) :
                outPts = SpherePts ( at.coord(), RAD, npts+i*2 )
                at_pts, at_pts_i = [None]*len(outPts), 0
                for pt in outPts :
                    vPt = [pt[0], pt[1], pt[2]]
                    apt = numpy.array ( vPt )
                    if allAtTree != None :
                        opointsNear = allAtTree.searchTree ( vPt, outRad )

                        if 1 :
                            foundNearPt = False
                            for npt in opointsNear :
                                v = apt - npt.coord().data()
                                r2 = numpy.sum ( v * v )
                                if r2 < outRad2 :
                                    foundNearPt = True
                                    break
                            if not foundNearPt :
                                at_pts[at_pts_i] = vPt
                                at_pts_i += 1

                        else :
                            if len(opointsNear) == 0 :
                                at_pts[at_pts_i] = vPt
                                at_pts_i += 1
                    else :
                        at_pts[at_pts_i] = vPt
                        at_pts_i += 1
                #if log :
                #    print " - %d, %d pts" % (i, len(at_pts))
                if at_pts_i >= npts or i >= 95 : # or show :
                    pts.extend ( at_pts[0:at_pts_i] )
                    break

        if show :
            AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f %s" % (RAD,atoms[0].name) )

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

            r_avg.append ( [RAD,0,0] )


        else :
            d_vals_n = dmap.interpolated_values ( pts, mol.openState.xform )
            d_vals = numpy.append ( d_vals, d_vals_n )
            avg = numpy.average ( d_vals_n )

            #gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #A, B = GV, 0
            #A, B = GV - minD, minD
            A,B = maxD - minD, minD
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            g_vals = numpy.append ( g_vals, numpy.ones([len(pts),1]) * gv )

            g_vals_avg = numpy.append ( g_vals_avg, gv )
            d_vals_avg = numpy.append ( d_vals_avg, avg )

            r_avg.append ( [RAD,avg,len(pts)] )

            #if log :
            #    print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, len(pts))

        RAD += dRAD
        i+=1

    if log and not fitg :
        min, max = r_avg[0][1], r_avg[0][1]
        for RAD, avg, numPts in r_avg :
            if avg < min : min = avg
            if avg > max : max = avg
        A,B = max-min, min
        #A,B = maxD - minD, minD
        #A,B = GV - minD, minD
        for RAD, avg, numPts in r_avg :
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            #print "%.1f\t%f\t%f\t%d" % (RAD, avg+0.02, gv+0.02, numPts)
            print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, numPts)

    #d_vals = d_vals + 0.02
    #g_vals = g_vals + 0.02

    # this is the CC between averaged radial values - not at robust
    if 0 :
        olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals_avg, g_vals_avg )
        if log :
            print "olap -avg-: %.3f cc: %.3f, Q: %.3f -- %d" % (olap, CC, Qs, len(d_vals_avg))
            #print "%f\t%f\t%f" % (olap, CC, Qs)

    olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals, g_vals )
    # this is the CC between _all_ radial values
    Qs = CCmean
    if log :
        print "olap --N--: %.3f cc: %.3f, ccmean (Q-score): %.3f -- %d" % (olap, CC, Qs, len(d_vals))
        #print "%f\t%f\t%f" % (olap, CC, Qs)


    if fitg :
        if log : print "fitting gaussian : "
        #V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        V, N = [ [x[0],x[1]] for x in r_avg[0:25] ], float(25)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f" % (sdev, A, B, err)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            #print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        return Qs, yds, err

    else :
        return Qs



# this is an older Q-score function which does not try to make sure to use numPts around each atom

def Qscore_ ( atoms, dmap, sigma, allAtTree = None, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    #sigma = 1.0

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)
    #print " __%s__ " % (atoms[0].name), pts[0]


    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals should have the reference gaussian...
    g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False)

    if mol == None :
        mol = atoms[0].molecule


    # r_avg holds the average values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( pts, mol.openState.xform ).astype(numpy.float64, copy=False)
    d_vals = numpy.repeat ( d_vals, numPts )

    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,len(pts)*numPts] ]

    d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)


    # make smaller atom tree...
    if 1 and allAtTree != None :
        ats_near = []
        for at in atoms :
            anear = allAtTree.searchTree ( at.coord().data(), toRAD*2.0 )
            ats_near.extend ( anear )

        points = _multiscale.get_atom_coordinates ( ats_near, transformed = False )
        if log :
            print " - new search tree: %d pts" % ( len(ats_near) )
        allAtTree = AdaptiveTree ( points.tolist(), ats_near, 1.0)



    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        #outRad2 = outRad * outRad
        pts = []
        for at in atoms :

            outPts = SpherePts ( at.coord(), RAD, numPts )
            at_pts, at_pts_i = [None]*len(outPts), 0

            for pt in outPts :
                vPt = [pt[0], pt[1], pt[2]]
                apt = numpy.array ( vPt )
                if allAtTree != None :
                    opointsNear = allAtTree.searchTree ( vPt, outRad )

                    if 1 :
                        foundNearPt = False
                        for npt in opointsNear :
                            v = apt - npt.coord().data()
                            r2 = numpy.sum ( v * v )
                            if r2 < outRad2 :
                                foundNearPt = True
                                break
                        if not foundNearPt :
                            at_pts[at_pts_i] = vPt
                            at_pts_i += 1

                    else :
                        if len(opointsNear) == 0 :
                            at_pts[at_pts_i] = vPt
                            at_pts_i += 1
                else :
                    at_pts[at_pts_i] = vPt
                    at_pts_i += 1

                pts.extend ( at_pts[0:at_pts_i] )

        if show :
            AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f" % RAD )

        if len (pts) < 1 :
            if 0 and log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

            r_avg.append ( [RAD,0,0] )


        else :
            d_vals_n = dmap.interpolated_values ( pts, mol.openState.xform )
            #d_vals = numpy.append ( d_vals, d_vals_n )
            avg = numpy.average ( d_vals_n )

            A,B = maxD - minD, minD
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B

            g_vals_avg = numpy.append ( g_vals_avg, gv )
            d_vals_avg = numpy.append ( d_vals_avg, avg )

            r_avg.append ( [RAD,avg,len(pts)] )


        RAD += dRAD
        i+=1

    if 0 and log :
        min, max = r_avg[0][1], r_avg[0][1]
        for RAD, avg, numPts in r_avg :
            if avg < min : min = avg
            if avg > max : max = avg
        A,B = max-min, min
        A,B = maxD - minD, minD
        #A,B = GV - minD, minD
        for RAD, avg, numPts in r_avg :
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            #print "%.1f\t%f\t%f\t%d" % (RAD, avg+0.02, gv+0.02, numPts)
            print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, numPts)

    #d_vals = d_vals + 0.02
    #g_vals = g_vals + 0.02

    olap, CC, CCm = FitMap.overlap_and_correlation ( d_vals_avg, g_vals_avg )
    Qscore = CCm
    if log :
        print "olap -avg-: %.3f cc: %.3f, ccm (Q-score): %.3f -- %d" % (olap, CC, CCm, len(d_vals_avg))
        #print "%f\t%f\t%f" % (olap, CC, CCm)


    if fitg :
        if log : print "fitting gaussian : "
        #V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        V, N = [ [x[0],x[1]] for x in r_avg[0:15] ], float(15)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f" % (sdev, A, B, err)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        return Qscore, yds, err

    else :
        return Qscore




def QscorePt ( atPt, xfI, dmap, sigma, allAtTree = None, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0 ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    #xfI = chimera.Xform()
    atPtC = chimera.Point ( *atPt )

    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals should have the reference gaussian...
    g_vals = (numpy.ones ( [numPts,1] ) * refG).astype(numpy.float64, copy=False )
    g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False )

    # r_avg holds the average values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( [atPt], xfI ).astype(numpy.float64, copy=False)
    d_vals = numpy.repeat ( d_vals, numPts )

    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,numPts] ]

    d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)


    # make smaller atom tree...
    if 1 and allAtTree != None :
        ats_near = []
        anear = allAtTree.searchTree ( atPt, toRAD*2.0 )
        ats_near.extend ( anear )

        points = _multiscale.get_atom_coordinates ( ats_near, transformed = False )
        if log :
            print " - new search tree: %d pts" % ( len(ats_near) )
        allAtTree = AdaptiveTree ( points.tolist(), ats_near, 1.0)

    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        #outRad2 = outRad * outRad
        pts = []

        for i in range (0, 100) :
            outPts = SpherePts ( atPtC, RAD, numPts+i*2 )
            at_pts, at_pts_i = [None]*len(outPts), 0
            for pt in outPts :
                vPt = [pt[0], pt[1], pt[2]]
                apt = numpy.array ( vPt )
                if allAtTree != None :
                    opointsNear = allAtTree.searchTree ( vPt, outRad )
                    foundNearPt = False
                    for npt in opointsNear :
                        v = apt - npt.coord().data()
                        r2 = numpy.sum ( v * v )
                        if r2 < outRad2 :
                            foundNearPt = True
                            break
                    if not foundNearPt :
                        at_pts[at_pts_i] = vPt
                        at_pts_i += 1

                else :
                    at_pts[at_pts_i] = vPt
                    at_pts_i += 1
            #if log :
            #    print " - %d, %d pts" % (i, len(at_pts))
            if at_pts_i >= numPts or i >= 15 : # or show :
                pts.extend ( at_pts[0:at_pts_i] )
                break

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

            r_avg.append ( [RAD,0,0] )


        else :
            d_vals_n = dmap.interpolated_values ( pts, xfI )
            d_vals = numpy.append ( d_vals, d_vals_n )
            avg = numpy.average ( d_vals_n )

            #gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #A, B = GV, 0
            #A, B = GV - minD, minD
            A,B = maxD - minD, minD
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            g_vals = numpy.append ( g_vals, numpy.ones([len(pts),1]) * gv )

            g_vals_avg = numpy.append ( g_vals_avg, gv )
            d_vals_avg = numpy.append ( d_vals_avg, avg )

            r_avg.append ( [RAD,avg,len(pts)] )

            #if log :
            #    print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, len(pts))

        RAD += dRAD
        i+=1

    if log and not fitg :
        min, max = r_avg[0][1], r_avg[0][1]
        for RAD, avg, numPts in r_avg :
            if avg < min : min = avg
            if avg > max : max = avg
        A,B = max-min, min
        #A,B = maxD - minD, minD
        #A,B = GV - minD, minD
        for RAD, avg, numPts in r_avg :
            gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
            #print "%.1f\t%f\t%f\t%d" % (RAD, avg+0.02, gv+0.02, numPts)
            print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, numPts)

    #d_vals = d_vals + 0.02
    #g_vals = g_vals + 0.02

    #if log :
    #    olap, CC, CCm = FitMap.overlap_and_correlation ( d_vals_avg, g_vals_avg )
    #    print "olap -avg-: %.3f cc: %.3f, ccm: %.3f -- %d" % (olap, CC, CCm, len(d_vals_avg))
    #    #print "%f\t%f\t%f" % (olap, CC, CCm)


    olap, CC, CCm = FitMap.overlap_and_correlation ( d_vals, g_vals )
    qscore = CCm
    if log :
        print "olap --N--: %.3f cc: %.3f, ccm: %.3f -- %d" % (olap, CC, CCm, len(d_vals))
        #print "%f\t%f\t%f" % (olap, CC, CCm)

    if fitg :
        if log : print "fitting gaussian : "
        #V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        V, N = [ [x[0],x[1]] for x in r_avg[0:25] ], float(25)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f" % (sdev, A, B, err)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            #print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        return qscore, yds, err

    else :
        return qscore



def RadAts ( atoms, dmap, allAtTree = None, show=0, log=0, numPts=20, toRAD=2.0, dRAD=0.1 ) :

    if len(atoms) == 0 :
        #print " - no RAD atoms?"
        return None

    #pts = []
    #for at in atoms :
    #    p = at.coord()
    #    pts.append ( [p[0], p[1], p[2]] )

    from _multiscale import get_atom_coordinates
    pts = get_atom_coordinates(atoms, transformed = False)

    RD_, X, Y = [], [], []
    d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
    avg = numpy.average ( d_vals )

    RD_.append ( [0,avg] ); X.append (0); Y.append (avg)


    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        pts = []
        for at in atoms :
            npts = (numPts * RAD*RAD / (dRAD*dRAD)) if show else numPts
            npts = int ( npts )
            #print RAD, dRAD, numPts, " -> ", npts
            outPts = SpherePts ( at.coord(), RAD, npts )
            for pt in outPts :
                ppt = [pt[0], pt[1], pt[2]]
                if allAtTree != None :
                    vPt = numpy.array ( ppt )
                    opointsNear = allAtTree.searchTree ( ppt, outRad )
                    if 1 :
                        clash = False
                        for p in opointsNear :
                            v = vPt - p.coord().data()
                            sqSum = numpy.sum ( v * v )
                            if sqSum < outRad2 :
                                clash = True
                                break
                        if clash == False :
                            pts.append ( ppt )

                    else :
                        if len(opointsNear) == 0 :
                            pts.append ( ppt )
                else :
                    pts.append ( ppt )

        if show :
            AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f" % RAD )

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - %d.%s - " % (RAD, atoms[0].residue.id.position, atoms[0].residue.type),
                print "SC" if atoms[0].isSC else "BB"

        else :
            d_vals = dmap.interpolated_values ( pts, atoms[0].molecule.openState.xform )
            avg = numpy.average ( d_vals )
            RD_.append ( [RAD,avg] );
            if log :
                print RAD, avg, len(pts)
                X.append (RAD); Y.append (avg)

        RAD += dRAD

    #minSd = opt0 ( RD_, 0.1 )
    #if minSd != None :
    #    if show :
    #        print " SD0: %.1f" % minSd

    sdev = toRAD
    slope = 0

    if RD_[0][1] <=  RD_[-1][1] :
        sdev = 10.0

    else :

        #for i in range ( len(RD_) ) :
        #    RD_[i][1] = RD_[i][1] - RD_[-1][1]
        #    if log :
        #        Y[i] = Y[i] - Y[-1]


        #import time
        #start = time.time()
        sdev, A, B = optSGD ( RD_, 9000, 0.2 )
        sdev, A, B = optSGD ( RD_, 9000, 0.02, sdev, A, B )
        sdev, A, B = optSGD ( RD_, 9000, 0.002, sdev, A, B )
        #end = time.time()
        #if log : print " sgd - sdev: %.4f, A %.4f, B %.4f -- %f" % (sdev, A, B, (end - start))
        sdev = sdev
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

        #start = time.time()
        #sdev, A, B = optGN ( RD_, 0.0001 )
        #print " gn - sdev: %.4f, A %.4f, B %.4f -- %f" % (sdev, A, B, (end - start))
        #end = time.time()

        if 1 :
            if 0 and sdev != None :

                if log :
                    print " gn1 - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

            else :
                sdev, A, B = optSGD ( RD_, 10000, 0.01 )

                if log :
                    print " sgd - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)

                sdev2, A2, B2 = optGN ( RD_, 0.0001, sdev, A, B )
                if sdev2 != None :
                    sdev, A, B = sdev2, A2, B2
                    if log :
                        print " gn2 - sdev: %.4f, A %.4f, B %.4f" % (sdev, A, B)
                #else :
                #    return 10.0


        if log :
            r = numpy.polyfit ( X, Y, 1, rcond=None, full=False, w=None, cov=False)
            print " sdev: %.4f, A %.4f, B %.4f // slope: %.4f y %.4f" % (sdev, A, B, r[0], r[1])

            #A, B = 0.26+0.08, -0.08
            lastX = 0
            for i in range ( len(RD_) ) :
                x, y = RD_[i]
                gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
                gvRef = A * numpy.exp ( -0.5 * numpy.power(x/0.5,2) ) + B
                lv = x * r[0] + r[1]
                print "%.1f\t%f\t%f\t%f" % (x, y, gv, gvRef)
                lastX = x

            if 1 :
                x = lastX + dRAD
                #while x < min(4 * sdev,50.0) :
                while x < min(10.0,50.0) :
                    gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
                    gvRef = A * numpy.exp ( -0.5 * numpy.power(x/0.5,2) ) + B
                    lv = x * r[0] + r[1]
                    print "%.1f\t\t%f\t%f" % (x, gv, gvRef)
                    x += dRAD


    #return abs(sdev), abs(slope)
    return abs(sdev)


def TimeLeftStr ( atI, totI, totSec ) :

    leftTime = ""
    leftSec = 0.0
    iPerSec = float(atI) / totSec
    if iPerSec > 0 :
        leftSec = float ( totI - atI ) / iPerSec
        leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
        leftSec = leftSec - leftHour * 60.0 * 60.0
        leftMin = numpy.floor ( leftSec / 60.0 )
        leftSec = leftSec - leftMin * 60.0
        leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)
        return leftTime
    return ""


def optGN ( V, err, S=None, A=None, B=None ) :

    y0 = V[0][1]
    yN = V[-1][1]

    if S == None :
        S = 0.5
        A = y0+yN
        B = yN

    an = numpy.array ( [A,B,S] )
    #print " _ -- A %.3f B %.3f s %.3f" % (A, B, S)

    reg = 1.0
    badMatCount = 0

    for i in range ( 1000 ) :

        J = numpy.zeros ( [len(V),3] )
        e = numpy.zeros ( [len(V),1] )

        err0 = 0
        j = 0
        for x,y in V :
            expv = numpy.exp ( -0.5 * numpy.power(x/S,2) )
            v = A * expv + B
            yd = v - y
            err0 += yd * yd
            #print "%.2f,%.2f/%.2f(%.2f)" % (x, y, v, yd),

            dA = expv
            dB = 1
            dS = A*x*x*numpy.power(S,-3) * expv
            J[j,:] = [dA, dB, dS]
            e[j,0] = yd
            j += 1

        Jt = numpy.transpose(J)

        try :
            J_ = numpy.dot ( numpy.linalg.inv ( numpy.dot(Jt,J) ), Jt )
        except :
            #print " - bad matrix?"
            #print numpy.dot(Jt,J)
            badMatCount += 1

            if badMatCount > 3 :
                return None, None, None

            from numpy import random as R
            an = numpy.array ( [R.random()*(y0+yN),R.random()*yN,R.random()*10.0] )
            A,B,S = an[0], an[1], an[2]
            #print " ? -- A %.3f B %.3f s %.3f" % (A, B, S)
            reg = 1.0

            continue

        ad = numpy.dot ( J_, e )
        ann = an - ( ad[:,0] * reg )
        A,B,S = ann[0], ann[1], ann[2]

        err1 = err3 ( V, S, A, B )
        #if err1 > err0 :
        #    reg = reg * 0.1
        #    if reg < err :
        #        break
        #else :
        an = ann
        #print " %d -- A %.3f B %.3f s %.3f - err %.3f, reg %.5f" % (i, A, B, S, err1, reg)

        if abs(err0 - err1) < err :
            #print " - done"
            break

        i += 1

    return S,A,B



def optSGD ( V, N, err, S=None, A=None, B=None ) :

    if S == None :
        y0 = V[0][1]
        yN = V[-1][1]
        S = 0.5
        A = y0+yN
        B = yN

    from numpy import random

    lastE = err3 ( V, S, A, B )
    #while True :
    for i in range(N) :

        S_ = S + random.normal ( 0, err ) # mean, sigma
        A_ = A + random.normal ( 0, err ) # mean, sigma
        B_ = B + random.normal ( 0, err ) # mean, sigma

        e = err3 ( V, S_, A_, B_ )
        #print "%d %.2f %f %f %.4f" % (i, sdAt, e, numpy.log(e), dd)
        if e < lastE :
            S, A, B = S_, A_, B_
            lastE = e

    return S,A,B


def err3 ( XYz, sd, A, B ) :

    y0 = XYz[0][1]
    err = 0
    #for x,y in XYz[1:] :
    for x,y in XYz :
        yd = y - A * numpy.exp ( -0.5 * numpy.power(x/sd,2) ) - B
        err += yd * yd
    #err /= float(len(XYz))
    return err



def err ( XYz, sd ) :

    y0 = XYz[0][1]
    err = 0
    for x,y in XYz[1:] :
        yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
        err += yd * yd
    #err /= float(len(XYz))
    return err


def opt0 ( RD_, dStep ) :

    sd = 0.1
    y0 = RD_[0][1]
    minSd, minErr, N = None, 1e99, float ( len(RD_)-1 )
    while sd < 10.0 :

        err = 0
        for x,y in RD_[1:] :
            yd = y - y0 * numpy.exp ( -0.5 * numpy.power(x/sd,2) )
            err += yd * yd
        err /= N

        #print err

        if err < minErr :
            minErr = err
            minSd = sd

        sd += dStep


def opt ( V, maxErr ) :

    dd = 1.0
    sdAt = 0.1
    lastE = err ( V, sdAt )
    #while True :
    for i in range(10000) :
        sdAt += dd
        e = err ( V, sdAt )
        #print "%d %.2f %f %f %.4f" % (i, sdAt, e, numpy.log(e), dd)
        if e >= lastE :
            dd *= -0.75
            if abs(dd) < maxErr :
                return sdAt
        lastE = e
    return sdAt








def Calc ( chimeraPath, numProc, res=3.0, bfactorF=-1, sigma=0.6 ) :

    print "Calc Q scores"
    print " - chimera path: ", chimeraPath
    print " - num processors: ", numProc
    print " - resolution: ", res
    print " - sigma: ", sigma
    if bfactorF > 0 :
        print " - b-factor: ", bfactorF


    from VolumeViewer import Volume
    vols = chimera.openModels.list(modelTypes = [Volume])
    if len(vols) == 0 :
        print " - no volumes loaded"
        return

    dmap = vols[0]
    print " - volume: %s" % dmap.name

    from chimera import Molecule

    mols = chimera.openModels.list(modelTypes = [Molecule])
    if len(mols) == 0 :
        print " - no molecules loaded"
        return

    for mi, mol in enumerate (mols) :

        print ""
        print "Model %d/%d: %s" % (mi+1, len(mols), mol.name)
        SetBBAts ( mol )

        ats = [at for at in mol.atoms if not at.element.name == "H"]
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
        #allAtTree = None

        if numProc == 1 :
            CalcQ ( mol, None, dmap, sigma, allAtTree=allAtTree )
        else :
            CalcQp ( mol, None, dmap, sigma, allAtTree=allAtTree, numProc=numProc )

        SaveQStats ( mol, "All", dmap, sigma, res )

        if bfactorF > 0 :
            minb, maxb = 1.0e9, 0.0
            for at in mol.atoms :
                at.bfactor = bfactorF * (1.0 - at.Q)
                #at.occupancy = 1.0 # max(0,at.Q)
                #dval = self.cur_dmap.interpolated_values ( [ at.coord()  ], self.cur_mol.openState.xform ).astype(numpy.float64, copy=False)[0]
                #at.occupancy = (dval - minD) / (maxD - minD)
                minb = min ( minb, at.bfactor )
                maxb = max ( maxb, at.bfactor )

            molPath = os.path.splitext(mol.openedAs[0])[0]
            nname = molPath + "_B%.0f.pdb" % bfactorF
            print "Saving pdb with B'-factors, f=%.0f:" % bfactorF
            print "  -> ", nname
            print "  - bfactor = %.0f*(1-Qscore), range %.2f to %.2f" % (bfactorF, minb, maxb)
            #print "  -  occupancies set to 1"
            print ""
            chimera.PDBio().writePDBfile ( [mol], nname )




# this is the function that the MP version executes once Chimera is opened
# with partial model and map

def CalcQForOpenModelsRess () :

    from VolumeViewer import Volume
    dmap = chimera.openModels.list(modelTypes = [Volume])[0]
    print " - dmap: %s" % dmap.name

    #minD, maxD = MinMaxD ( dmap )
    #print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mol = chimera.openModels.list(modelTypes = [Molecule])[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )


    #rids = {}
    #for r in mol.residues :
    #    rids["%d.%s" % (r.id.position,r.id.chainId)] = r

    atids = {}
    for r in mol.residues :
        for at in r.atoms :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            atids["%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,altLoc)] = at


    ats = [at for at in mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
    allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)


    fin = open ( os.path.splitext ( dmap.data.path )[0] + ".txt" )
    fout = open ( os.path.splitext ( dmap.data.path )[0] + "_out.txt", "w" )
    foutn = os.path.splitext ( dmap.data.path )[0] + "_stat.txt"

    sig_at = []

    for l in fin :
        #print l,
        sigma, minD, maxD, atIdStr = l.split()
        if not atIdStr in atids :
            print " - atid not found: ", atIdStr
        at = atids[atIdStr.strip()]
        sigma = float(sigma)
        minD, maxD = float(minD), float(maxD)
        sig_at.append ( [sigma, minD, maxD, at, atIdStr] )

    fs = open ( foutn, "w" ); fs.write ( "%d/%d" % (0,len(sig_at) ) ); fs.close()

    import time
    start = time.time()

    i = 0
    for sigma, minD, maxD, at, atId in sig_at :
        #print "%d.%s.%s" % (r.id.position,r.id.chainId,at.name),

        sig = sigma

        # Q-scores for ions and water are using sigma of 0.4
        #if at.name.upper() in chargedIons or at.residue.type.upper() == "HOH" :
        #    sig = 0.4

        qs = Qscore ( [at], dmap, sig, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
        fout.write ( "%s %f\n" % (atId, qs) )

        if i%10 == 0 :

            end = time.time()
            totSec = end - start

            leftTime = ""
            leftSec = 0.0
            iPerSec = float(i) / totSec
            if iPerSec > 0 :
                leftSec = float ( len(sig_at) - i ) / iPerSec
                leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                leftSec = leftSec - leftHour * 60.0 * 60.0
                leftMin = numpy.floor ( leftSec / 60.0 )
                leftSec = leftSec - leftMin * 60.0
                leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)


            fs = open ( foutn, "w" ); fs.write ( "%d/%d - %s" % (i+1,len(sig_at),leftTime) ); fs.close()

        i += 1


    fin.close()
    fout.close()

    fs = open ( foutn, "w" ); fs.write ( "done" ); fs.close()




def CalcQp ( mol, cid, dmap, sigma, allAtTree=None, useOld=False, log=False, numProc=None ) :


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    nname = molPath + "__Q__" + mapName + ".pdb"

    if useOld :
        SetBBAts ( mol )
        if QsFromFile ( mol, nname ) :
            Qavg = QStats1 ( mol, cid )
            return Qavg


    #numProc = 2

    if numProc == None :
        import multiprocessing
        numProc = multiprocessing.cpu_count() / 2

    M = dmap.data.full_matrix()
    minD, maxD = numpy.min(M), numpy.max(M)

    print "Q Scores - p - %d" % numProc
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    import time
    start = time.time()

    SetBBAts ( mol )

    ress = []
    atoms = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            if 1 or r.isNA :
                for at in r.atoms :
                    if 0 or not at.element.name == "H" :
                        atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    import subprocess
    import sys
    mapPath = os.path.split ( dmap.data.path )[0]
    mapBase = os.path.splitext ( dmap.data.path )[0]

    print "cmd:",
    #print sys.argv
    for arg in sys.argv :
        print arg,
    print ""

    # '/Users/greg/_mol/Chimera.app/Contents/Resources/share/__main__.py'
    chiPath = os.path.split ( sys.argv[0] )[0]
    #mapQPPath = os.path.join ( chiPath, 'Segger' )
    #mapQPPath = os.path.join ( chiPath, 'mapqp.py' )
    #print " -- path to mapQ script:", mapQPPath

    # for Mac
    chiPath, share = os.path.split ( chiPath )
    #print chiPath, share
    chiPath2, resOrChim = os.path.split ( chiPath )
    #print chiPath, resOrChim
    if "Chimera" in resOrChim :
        print " -- on unix"
        chiPath = os.path.join ( chiPath, 'bin' )
        chiPath = os.path.join ( chiPath, 'chimera' )
    else :
        print " -- on mac"
        #chiPath2, contents = os.path.split ( chiPath2 )
        #print chiPath2, contents
        chiPath = os.path.join ( chiPath2, 'MacOS' )
        chiPath = os.path.join ( chiPath, 'chimera' )

    print " -- path to Chimera:", chiPath

    dir_path = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.split(dir_path)[0]
    print " -- working dir:", inDir
    #mapQPPath = os.path.join ( inDir, 'Segger' )
    mapQPPath = os.path.join ( dir_path, 'mapqp.py' )
    print " -- path to mapQ script:", mapQPPath

    mapBase = mapBase + "_Q-score-mp"

    n = len(atoms)
    g = [atoms[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]
    procs = []
    for mi, atoms1 in enumerate(g) :

        ress1 = atoms1[0].residue
        ressN = atoms1[-1].residue
        print " - %d/%d, %d-%d" % (mi+1, numProc, ress1.id.position, ressN.id.position)

        fout = open ( mapBase + "_%d.txt" % mi, "w" )
        for at in atoms1 :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            fout.write ( "%.3f %f %f %d.%s.%s.%s\n" % (sigma, minD, maxD, r.id.position,r.id.chainId,at.name,altLoc) )
        fout.close()

        nmap_path = mapBase + "_%d.mrc" % mi
        #print " -> ", nmap_path
        nmap = MaskMapResize ( atoms1, 4.0, dmap, nmap_path )
        #nmap.write_file ( nmap_path , "mrc" )

        args = [chiPath, '--nogui', '--silent', '--nostatus', mol.openedAs[0], nmap_path, mapQPPath]
        if mi == 0 :
            print "running proc:",
            for arg in args :
                print arg,
            print ""

        fout = open ( mapBase + "_%d.log" % mi, "w" )
        foute = open ( mapBase + "_%d_err.log" % mi, "w" )
        p = subprocess.Popen(args, stdout=fout, stderr=foute, cwd=inDir)
        procs.append ( [mi, p, fout, foute] )

    print ""
    print "Waiting...",
    for mi, p, fout, foute in procs :
        p.wait()
        fout.close()
        foute.close()
        print "%d" % mi,
    print ""

    atids = {}
    for r in mol.residues :
        for at in r.atoms :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            atids["%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,altLoc)] = at

    print ""
    print "Getting...",
    for mi, p, fout, foute in procs :
        fin = mapBase + "_%d_out.txt" % mi
        #print " - getting from: ", fin
        fp = open ( fin )
        for l in fp :
            #print " - ", l
            try :
                atId, Q = l.split()
            except :
                print " - err line: ", l
                blah
            at = atids[atId.strip()]
            #at = r.atomsMap[atName][0]
            at.Q = float(Q)
            #at.CC = float(cc)
            at.bfactor = at.Q

        fp.close()

        if mi == 0 :
            print ""
            print ""
            print "__Out %d__" % mi
            foute = open ( mapBase + "_%d.log" % mi, "r" )
            for l in foute :
                print l,
            print ""
            foute.close()


            print "__Err %d__" % mi
            foute = open ( mapBase + "_%d_err.log" % mi, "r" )
            for l in foute :
                print l,
            print ""
            foute.close()

        if 1 :
            #print " - removing..."
            os.remove ( mapBase + "_%d_out.txt" % mi )
            try :
                os.remove ( mapBase + "_%d_stat.txt" % mi )
            except :
                print " - did not find _stat file"
                pass
            os.remove ( mapBase + "_%d.txt" % mi )
            os.remove ( mapBase + "_%d.mrc" % mi )
            os.remove ( mapBase + "_%d.log" % mi )
            os.remove ( mapBase + "_%d_err.log" % mi )
            print "%d" % mi,

    print ""


    end = time.time()
    print ""
    print " - done, time: %f" % ( end-start )
    totSec = end - start
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]

    nname = molPath + "__Q__" + mapName + ".pdb"
    print "Saving pdb with Q-scores:", nname
    chimera.PDBio().writePDBfile ( [mol], nname )

    Qavg = QStats1 ( mol, cid )

    return Qavg




def QStats1 ( mol, chainId ) :

    totQ, totN = 0.0, 0.0
    #QT, QN = { "Protein":0.0, "Nucleic":0.0, "Other":0.0 }, { "Protein":0.0, "Nucleic":0.0, "Other":0.0}
    QT, QN = {}, {}
    QT_, QN_ = {}, {}
    QH, QL = {}, {}

    doRess = []

    for r in mol.residues :
        #if r.id.chainId == chainId or chainId == None :
        doRess.append ( r )

    print ""
    print "Q for %d res..." % ( len(doRess) )
    for r in doRess :

        #if not r.isNA : continue
        #if not r.isProt : continue

        CalcResQ (r, None, None, useOld=True )

        for at in r.atoms :
            if at.element.name == "H" :
                continue

            if hasattr ( at, "Q") :
                totQ += at.Q
                totN += 1.0

                tp = "Other"
                if at.residue.isProt : tp = "Protein"
                elif at.residue.isNA : tp = "Nucleic"
                else : tp = at.residue.type

                if tp in QT :
                    QT[tp] += at.Q;
                    QN[tp] += 1.0;
                    QH[tp] = max(QH[tp], at.Q)
                    QL[tp] = min(QL[tp], at.Q)
                else :
                    QT[tp] = at.Q; QN[tp] = 1.0
                    QH[tp] = at.Q; QL[tp] = at.Q

                tps = r.id.chainId + ":" + tp
                if tps in QT_ :
                    QT_[tps] += at.Q; QN_[tps] += 1.0
                else :
                    QT_[tps] = at.Q; QN_[tps] = 1.0


    #for tp in ["Other", "Protein", "Nucleic"] :
    print ""
    print "Chain\tAvg.Q-score\tEst.Res.(A)"
    tpk = QT_.keys()
    tpk.sort()
    for tp in tpk :
        if QN_[tp] > 0 :
            avgQ = QT_[tp]/QN_[tp]
            avgR = 0
            if "nucleic" in tp.lower() :
                avgR = (avgQ-1.0673)/-0.1574
            else :
                avgR = (avgQ-1.1244)/-0.1794
            print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
        else :
            print " %s\tn/a" % (tp)

    Q__ = { " protein":0, " nucleic":0, " water":0, " ion":0 }

    #for tp in ["Other", "Protein", "Nucleic"] :
    print ""
    print "Type\tAvg.Q-score\tEst.Res.(A)"
    for tp in QT.keys() :
        if QN[tp] > 0 :
            avgQ = QT[tp]/QN[tp]
            avgR = 0
            if "nucleic" in tp.lower() :
                avgR = (avgQ-1.0673)/-0.1574
                Q__[" nucleic"] = avgQ
            elif "protein" in tp.lower() :
                avgR = (avgQ-1.1244)/-0.1794
                Q__[" protein"] = avgQ
            elif "hoh" in tp.lower() :
                avgR = (avgQ-1.1244)/-0.1794
                Q__[" water"] = avgQ
            elif tp.upper() in chargedIons :
                avgR = (avgQ-1.1244)/-0.1794
                Q__[" ion"] = avgQ
            else :
                avgR = (avgQ-1.1244)/-0.1794
                Q__[tp] = avgQ
            print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
        else :
            print " %s\tn/a" % (tp)

    print ""

    for tp in QT.keys() :
        if QN[tp] > 0 :
            print "\t%s" % tp,
    print ""

    print "Avg.Q.",
    for tp in QT.keys() :
        if QN[tp] > 0 :
            avgQ = QT[tp]/QN[tp]
            print "\t%.3f" % avgQ,
    print ""

    print "Max.Q.",
    for tp in QT.keys() :
        if QN[tp] > 0 :
            print "\t%.3f" % QH[tp],
    print ""

    print "Min.Q.",
    for tp in QT.keys() :
        if QN[tp] > 0 :
            print "\t%.3f" % QL[tp],
    print ""

    print ""

    return Q__



def SaveQStats ( mol, chainId, dmap, sigma, RES=3.0 ) :

    cres = {}
    for r in mol.residues :
        try :
            cres[r.id.chainId].append ( [r.id.position, r] )
        except :
            cres[r.id.chainId] = [ [r.id.position, r] ]


    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    nname = molPath + "__Q__" + mapName + "_" + chainId + ".txt"
    #nname = molPath + "__Q__" + mapName + "_" + cid + ".txt"
    fp = open (nname, "w")

    print ""
    print "Saving per-chain & per-residue Q-scores:"
    print " -> res=", RES
    print " -> file:", nname

    fp.write ( "Chain\tQ_chain\tEst.Res.\tExpectedQ@%.2f\n" % RES )

    chains = cres.keys()
    chains.sort()

    for cid in chains :
        if 0 or cid == chainId or chainId == "All" :

            tps = {}
            resAtoms = []
            rs = cres[cid]
            for ri, r in rs :
                resAtoms.extend ( r.atoms )
                tp = "Other"
                if r.isProt : tp = "Protein"
                elif r.isNA : tp = "Nucleic"
                elif r.type.upper() in chargedIons : tp = "Ion"
                elif r.type.upper() == "HOH" : tp = "Water"
                tps[tp] = 1

            ctypes = ""
            for tp in tps.keys() :
                ctypes = (ctypes + tp) if len(ctypes) == 0 else (ctypes + "," + tp)

            cQ = numpy.average ( [at.Q for at in resAtoms if at.element.name != "H"] )

            formula, estRes = None, None
            if "Protein" in ctypes :
                formula = "=-0.1775 * %.2f + 1.1192" % RES
                estRes = (cQ - 1.1192) / -0.1775
            elif "Nucleic" in ctypes :
                formula = "= -0.1377 * %.2f + 0.9973" % RES
                estRes = (cQ - 0.9973) / -0.1377
            elif "Ion" in ctypes :
                formula = "= -0.1103 * %.2f + 1.0795" % RES
                estRes = (cQ - 1.0795) / -0.1103
            elif "Water" in ctypes :
                formula = "= -0.0895 * %.2f + 1.0001" % RES
                estRes = (cQ - 1.0001) / -0.0895

            fp.write ( "%s\t%.2f\t%.2f\t%s\t(%s)\n" % (cid, cQ, estRes, formula, ctypes) )

            #print " - cid: %s - %s - %.2f" % (cid, ctypes, cQ)

    fp.write ( "\n" )
    fp.write ( "Sigma: %g\n" % sigma )
    fp.write ( "\n" )
    fp.write ( "Protein: avgQ = -0.1775 * RES + 1.1192\n" )
    fp.write ( "Nucleic: avgQ = -0.1377 * RES + 0.9973\n" )
    fp.write ( "Ion: avgQ = -0.1103 * RES + 1.0795\n" )
    fp.write ( "Water: avgQ = -0.0895 * RES + 1.0001\n" )
    fp.write ( "\n" )

    avgQrna = -0.1574 * RES + 1.0673 # rna
    avgQprot = -0.1794 * RES + 1.1244 # protein
    avgQIon =  -0.1103 * RES + 1.0795 # ion
    avgQWater =  -0.0895 * RES + 1.0001 # water

    for cid in cres.keys () :

        if cid == chainId or chainId == "All" :

            fp.write ( "Chain %s\t\t\t\t\t\t\t\tAverage over 1 residue\t\t\t\t\tAverage over 2 residues\t\t\t\t\tAverage over 3 residues\t\t\t\t\tAverage over 5 residues\n\n" % cid )

            fp.write ( "Chain\tRes\tRes #\tQ_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-1)\tQ_sideChain(avg-1)\tQ_residue(avg-1)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-2)\tQ_sideChain(avg-2)\tQ_residue(avg-2)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-3)\tQ_sideChain(avg-3)\tQ_residue(avg-3)\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone(avg-5)\tQ_sideChain(avg-5)\tQ_residue(avg-5)\tExpectedQ@%.2f\t\n" % RES )

            #cid = 'A'
            rs = cres[cid]

            #print " - cid: %s - " % (cid)

            rs.sort()
            #for i in range (10) :
            #    print rs[i]

            ress = []
            Qs, AV, CC = [], [], []
            for ri, r in rs :

                #if not r.isProt and not r.isNA :
                #    print " - cid: %s - r %d - not prot or RNA" % (cid, r.id.position)
                #    continue

                ress.append (r)

                r.Q = numpy.average ( [at.Q for at in r.atoms if at.element.name != "H"] )

                r.qBB, r.qSC = 0, 0
                if len(r.bbAtoms) > 0 :
                    r.qBB = numpy.average ( [at.Q for at in r.bbAtoms if at.element.name != "H"] )
                if len(r.scAtoms) > 0 :
                    r.qSC = numpy.average ( [at.Q for at in r.scAtoms if at.element.name != "H"] )
                Qs.append ( [r.qBB, r.qSC, r.Q] )

                if 0 :
                    ad = avgdAts ( r.atoms, dmap )
                    aSC, aBB = 0, 0
                    if len(r.scAtoms) > 0 :
                        aSC = avgdAts ( r.scAtoms, dmap )
                    if len(r.bbAtoms) > 0 :
                        aBB = avgdAts ( r.bbAtoms, dmap )
                    AV.append ( [ad, aBB, aSC] )

                if 0 :
                    cc, ccm = ccAts ( r.atoms, dmap, RES )
                    ccSC, ccmSC = ccAts ( r.scAtoms, dmap, RES )
                    ccBB, ccmBB = ccAts ( r.bbAtoms, dmap, RES )
                    CC.append ( [cc, ccBB, ccSC] )
                    #CC.append ( [ccm, ccmBB, ccmSugar, ccmBase] )



            def N ( A, i, ind, N ) :
                #for i, a in enumerate ( A ) :
                sum, n = 0, 0
                for j in range ( i-N, i+N+1 ) :
                    if j >= 0 and j < len(A) :
                        sum += A[j][ind]
                        n += 1.0
                return sum/n


            last_i = None
            for i, r in enumerate ( ress ) :

                # fills in missing residues in proteins and rna
                if (r.isNA or r.isProt) and last_i != None :
                    ii = last_i+1
                    while ii < r.id.position :
                        avgQ = avgQrna if r.isNA else avgQprot
                        fp.write ( "%s\t%s\t%d\t\t\t\t%f\t\t" % (r.id.chainId, "", ii, avgQ ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\t\t" % (avgQ) )
                        fp.write ( "\t\t\t%f\n" % (avgQ) )
                        ii += 1

                if r.isNA :
                    fp.write ( "%s\t%s\t%d\t%f\t%f\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.qSC, r.Q, avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,1,1), N(Qs,i,2,1), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,1,2), N(Qs,i,2,2), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,1,3), N(Qs,i,2,3), avgQrna ) )
                    fp.write ( "%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,1,5), N(Qs,i,2,5), avgQrna ) )
                elif r.isProt :
                    if len(r.scAtoms) > 0 :
                        fp.write ( "%s\t%s\t%d\t%f\t%f\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.qSC, r.Q, avgQprot ) )
                        fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,1,1), N(Qs,i,2,1), avgQprot ) )
                        fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,1,2), N(Qs,i,2,2), avgQprot ) )
                        fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,1,3), N(Qs,i,2,3), avgQprot ) )
                        fp.write ( "%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,1,5), N(Qs,i,2,5), avgQprot ) )
                    else :
                        fp.write ( "%s\t%s\t%d\t%f\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.Q, avgQprot ) )
                        fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,2,1), avgQprot ) )
                        fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,2,2), avgQprot ) )
                        fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,2,3), avgQprot ) )
                        fp.write ( "%f\t\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,2,5), avgQprot ) )
                elif r.type.upper() in chargedIons :
                    fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.Q, avgQIon ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,1), avgQIon ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,2), avgQIon ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,3), avgQIon ) )
                    fp.write ( "\t\t%f\t%f\n" % (N(Qs,i,2,5), avgQIon ) )
                elif r.type.upper() == "HOH" :
                    fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.Q, avgQWater ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,1), avgQWater ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,2), avgQWater ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,3), avgQWater ) )
                    fp.write ( "\t\t%f\t%f\n" % (N(Qs,i,2,5), avgQWater ) )
                else :
                    fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.Q, avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,1), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,2), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\t\t" % (N(Qs,i,2,3), avgQprot ) )
                    fp.write ( "\t\t%f\t%f\n" % (N(Qs,i,2,5), avgQprot ) )

                last_i = r.id.position


            fp.write ( "\n\n" )


    fp.close()
    print ""






def CalcRadZ ( mol, cid, dmap, allAtTree, useOld=False, log=False ) :


    print "Rad-Z Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")


    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if not useOld :
                ress.append ( r )
            elif not hasattr (r, 'scS' ) :
                ress.append ( r )

    print " - residues to do: %d" % len(ress)


    for ri, r in enumerate ( ress ) :

        r.scZ = RadZ ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2 )
        r.bbZ = RadZ ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2 )

        if log and ri % 10 == 0 :
            status ( "Calculating - res %d/%d" % (ri, len(ress)) )
            print ".",


    scoresBB, scoresSC = [], []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if r.bbZ != None :
                scoresBB.append ( r.bbZ )
            if r.scZ != None :
                scoresSC.append ( r.scZ )

    print " - avg radz - side chain %.1f, backbone %.1f" % (numpy.average(scoresSC), numpy.average(scoresBB) )

    return numpy.average(scoresBB), numpy.average(scoresSC)




def qwork (num, ress, dmap, allAtTree, log):

    print 'qwork %d - %d res, %d - %d' % (num, len(ress), ress[0].id.position, ress[-1].id.position)

    for ri, r in enumerate ( ress ) :
        r.scZ = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )
        r.bbZ = RadAts ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )

        if num == 0 and log :
            status ( "Calculating Q scores - %d/%d" % (ri, len(ress)) )
            print ".",



def CalcSigma ( mol, cid, dmap, allAtTree, useOld=False, log=False ) :


    print "Sigma Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            if not useOld :
                ress.append ( r )
            elif not hasattr (r, 'scS' ) :
                ress.append ( r )

    print " - residues to do: %d" % len(ress)



    if 0 :

        import multiprocessing, threading
        N = 4 # multiprocessing.cpu_count()
        print " - cores: %d" % N
        dn = len(ress) / N

        threads = []
        for i in range(N):
            l = i * dn
            h = (i+1)*dn if i != N-1 else len(ress)
            #print "t %d, %d-%d" % (i, l, h)

            #t = threading.Thread(target=qwork, args=(i,ress[l:h], dmap, allAtTree))
            #threads.append(t)
            #t.start()

            #t = threading.Thread(name='d%d'%i, target=qwork, args=(i,ress[l:h], dmap, allAtTree, log))
            #t.setDaemon(True)
            #t.start()
            #threads.append(t)

            #print __name__
            if 1 or __name__ == '__main__':
                p = ctx.Process(target=qwork, args=(i,ress[l:h], dmap, allAtTree, log))
                p.start()
                threads.append(p)

        for i, t in enumerate(threads) :
            print "j %d" % (i)
            t.join()

    else :

        for ri, r in enumerate ( ress ) :

            r.bbZ = RadAts ( r.bbAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )
            r.scZ = RadAts ( r.scAtoms, dmap, allAtTree=allAtTree, show=0, log=0, numPts=10, toRAD=2, dRAD=0.2 )

            if log and ri % 10 == 0 :
                status ( "Calculating - res %d/%d" % (ri, len(ress)) )
                print ".",



    scoresBB, scoresSC = [], []

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            ress.append ( r )
            if r.bbZ != None : scoresBB.append ( r.bbZ )
            if r.scZ != None : scoresSC.append ( r.scZ )

    #sc = [x for x in scores if x is not None]
    #scSC = [1.0/x for x in scoresSC if x is not None]
    #scBB = [1.0/x for x in scoresBB if x is not None]

    #print " - %d res, SC min %.2f max %.2f, avg %.2f" % (len(ress), min(scSC), max(scSC), numpy.average(scSC))
    print " - avg sigma - side chain %.1f, backbone %.1f" % (numpy.average(scoresSC), numpy.average(scoresBB) )


    if 0 :

        sByType = {}
        rByType = {}
        for r in ress :
            if r.scZ != None :
                if not r.type in sByType :
                    rByType[r.type] = []
                    sByType[r.type] = []
                rByType[r.type].append ( [r.scZ, r] )
                sByType[r.type].append ( [r.scZ] )

        avgs = []
        for rtype, ra in sByType.iteritems () :
            avgs.append ( [numpy.average (ra), rtype] )

        from chimera.resCode import protein3to1
        from chimera.resCode import nucleic3to1
        avgs.sort ( reverse=True, key=lambda x: x[0] )


        mapName = os.path.splitext(dmap.name)[0]
        molName = os.path.splitext(mol.name)[0]
        mdir, mpfile = os.path.split(dmap.data.path)
        foname = mdir + "/" + mapName + "__" + molName + ".txt"


        print " - scores to: " + foname
        fp = open (foname,"w")

        for avgScore, rtype in avgs :

            rscores = rByType[rtype]
            rscores.sort ( reverse=False, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)

            rts = ""
            if R.isProt : rts = protein3to1[rtype]
            else : rts = nucleic3to1[rtype]

            print "%s\t%s\t%d\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, R.id.position, R.id.chainId, highestScore)
            fp.write ( "%s\t%s\t%d\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, R.id.position, R.id.chainId, highestScore) )

        fp.close()


    return numpy.average(scoresBB), numpy.average(scoresSC)


def CalcResQ (r, dmap, sigma, allAtTree=None, numPts=8, toRAD=2.0, dRAD=0.1, minD=0.0, maxD=1.0, useOld=False ) :

    scQ, bbQ, Q, numSC, numBB = 0.0, 0.0, 0.0, 0.0, 0.0
    for at in r.atoms :

        if at.element.name == "H" :
            continue

        if not hasattr ( 'isBB' ) :
            SetBBAts ( at.molecule )

        if not hasattr ( at, 'Q' ) or not useOld :
            #qs = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=numPts, toRAD=toRAD, dRAD=dRAD, minD=minD, maxD=maxD )
            at.Q = 0
            at.CC = 0

        Q += at.Q
        if r.isProt or r.isNA :
            if at.isBB :
                bbQ += at.Q
                numBB += 1.0
            else :
                scQ += at.Q
                numSC += 1.0

    if r.isProt or r.isNA :
        if int(numSC) != len(r.scAtoms) :
            print " - res %d.%s.%s - %.0f/%d sc atoms" % (r.id.position,r.type,r.id.chainId, numSC, len(r.scAtoms))

        if numSC > 0 :
            r.scQ = scQ / numSC
        else :
            r.scQ = None

        if numBB > 0 :
            r.bbQ = bbQ / numBB
        else :
            r.bbQ = None

    r.Q = Q / float ( len(r.atoms) )



def CalcQ_ ( mol, cid, dmap, sigma=0.5, allAtTree=None, useOld=False, log=False ) :

    print "Q Scores - in parallel"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")

    ress = []
    for r in mol.residues :
        if cid == None or r.id.chainId == cid :
            ress.append ( r )

    print " - residues to do: %d" % len(ress)


    import multiprocessing
    threads = multiprocessing.cpu_count() / 2
    print 'calc q using %d threads' % threads

    # Avoid periodic Python context switching.
    import sys
    original_check_interval = sys.getcheckinterval()
    sys.setcheckinterval(1000000000)

    # Define thread class for fitting.
    from threading import Thread
    class Q_Thread(Thread):
        def __init__(self, ress, ti):
            Thread.__init__(self)
            self.ress = ress
            self.ti = ti
        def run(self):
            print "run - %d - %d" % (self.ti, len(ress))
            for ri, r in enumerate ( self.ress ) :
                #CalcResQ (r, dmap, sigma, allAtTree=allAtTree, numPts=2, toRAD=2.0, dRAD=0.2 )
                #print "%d-%d/%d" % (ti,ri/len(self.ress)),
                for at in r.atoms :
                    if at.element.name != "H" :
                        qs = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5 )


    # Starts threads with each calculating an equal number of fits.
    n  = len(ress)
    g = [ress[(n*c)/threads:(n*(c+1))/threads] for c in range(threads)]
    threads = []
    for mi, ml in enumerate(g) :
        #print "%d - %d, %d-%d" % (mi,len(ml),ml[0].id.position,ml[-1].id.position)
        t = Q_Thread(ml,mi)
        threads.append(t)

    for t in threads:
        t.start()
    print ""

    # Wait for all threads to finish
    for t in threads:
        t.join()

    # Restore periodic context switching.
    sys.setcheckinterval(original_check_interval)

    # Collect fit results from all threads.
    #for t in threads:
    #    print "",




def CalcQ ( mol, cid, dmap, sigma, allAtTree=None, useOld=False, log=False ) :

    minD, maxD = MinMaxD ( dmap )

    print ""
    print "Q Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    SetBBAts ( mol )

    atoms = []

    import time
    start = time.time()

    #ress = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            for at in r.atoms :
                if at.element.name == "H" :
                    continue
                atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    #for ai, at in enumerate ( atoms[0:2] ) :
    #    qs = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )

    from chimera import tasks, CancelOperation
    task = tasks.Task('Calculating Q-scores', modal = True)

    try :

        for ai, at in enumerate ( atoms ) :

            at.Q = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
            at.bfactor = at.Q

            end = time.time()
            totSec = end - start

            leftTime = ""
            leftSec = 0.0
            iPerSec = float(ai) / totSec
            if iPerSec > 0 :
                leftSec = float ( len(atoms) - ai ) / iPerSec
                leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                leftSec = leftSec - leftHour * 60.0 * 60.0
                leftMin = numpy.floor ( leftSec / 60.0 )
                leftSec = leftSec - leftMin * 60.0
                leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)


            if (ai+1) % 100 == 0 :
                if log :
                    print "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime)

            task.updateStatus( "Calculating Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )

    except :
        print " - something went wrong..."
        return None

    finally :
        task.finished()


    end = time.time()
    print ""
    print " - done, time: %f" % ( end-start )
    totSec = end - start
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]

    try :
        nname = molPath + "__Q__" + mapName + ".pdb"
        chimera.PDBio().writePDBfile ( [mol], nname )
        print " - saved %s with Q-scores in O column" % nname
    except :
        print " - could not save Q-scores file"
        pass

    Qavg = QStats1 ( mol, cid )

    return Qavg








def QsFromFile ( mol, nname ) :

    rids = {}
    for r in mol.residues :
        rids["%d.%s" % (r.id.position,r.id.chainId)] = r


    # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    try :
        fin = open ( nname, "r" )
    except :
        #print " - file not found"
        return False

    print " - Qs from file: %s" % nname

    for line in fin :
        if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
            aname, aloc, cid, resi, occ, bfac = line[12:16].strip(), line[16:17].strip(), line[21], int(line[22:26]), float ( line[54:60] ), float ( line[60:66] )
            #if occ < 1.0 :
            rid = "%s.%s" % (resi,cid)
            if rid in rids :
                r = rids[rid]

                if aname in r.atomsMap :
                    ats = r.atomsMap[aname]
                    found = False
                    for at in ats :
                        if at.altLoc == aloc :
                            at.Q = bfac
                            at.bfactor = at.Q
                            #at.bfactor = 100.0 * (1.0 - at.Q)
                            found = True
                    if not found :
                        #print " -xx- %s.%s - atom %s - loc %s" % (resi, cid, aname, aloc)
                        continue
                else :
                    #print " -xx- %s.%s - atom %s" % (resi,cid, aname)
                    continue


    fin.close ()

    return True


def QScoreFileName ( mol, dmap ) :

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    nname = molPath + "__Q__" + mapName + ".pdb"

    return nname






def AddSpherePts ( pts, clr, rad, mname = "RAD points" ) :

    from chimera import elements, Coord, Atom, MolResId

    ptsMol = GetMod ( mname )

    res = None
    if ptsMol == None:
        from chimera import Molecule, openModels
        ptsMol = Molecule()
        ptsMol.name = mname
        ptsMol.isRealMolecule = False
        openModels.add ( [ptsMol], noprefs = True )
        res = ptsMol.newResidue('marker', chimera.MolResId('1', 1) )
    else :
        res = ptsMol.residues[0]

    for pt in pts :
        a = ptsMol.newAtom('', elements.H)
        res.addAtom(a)

        a.setCoord ( chimera.Point(*pt) )  # ( chimera.Point(*xyz) )
        a.radius = rad
        a.drawMode = Atom.Sphere
        a.color = chimera.MaterialColor ( *clr )
        a.surfaceCategory = 'markers'



def SpherePts ( ctr, rad, N ) :

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    pts = [None] * N
    for i, theta, phi in zip ( range(N), thetas, phis ):
        v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
        #if numpy.abs ( v.length - 1.0 ) > 1e-3 :
        #    print "x"
        pt = ctr + v * rad
        pts[i] = pt

    return pts



import threading



def Calc_ ( label="", res=0.0 ) :

    print "Calc Q scores:", label

    from VolumeViewer import Volume
    vols = chimera.openModels.list(modelTypes = [Volume])
    if len(vols) == 0 :
        print " - no volumes loaded"
        return
    dmap = vols[0]
    print " - dmap: %s" % dmap.name
    print " - res: %s" % res

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mols = chimera.openModels.list(modelTypes = [Molecule])
    if len(mols) == 0 :
        print " - no molecules loaded"
        return
    mol = mols[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )

    ats = [at for at in mol.atoms if not at.element.name == "H"]
    points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
    #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
    allAtTree = None


    qs, dr, q, qcc, emr = 0,0,0,0,0
    #bbRadZ, scRadZ, scRotaZ = 0,0,0

    sigma = 0.4

    cid = None
    #cid = mol.residues[0].id.chainId

    qs = CalcQp ( mol, cid, dmap, sigma=sigma, allAtTree=allAtTree, useOld=False )

    print ""
    print "Avg. Q scores:"
    print ""
    tps = qs.keys()
    tps.sort()
    for tp in tps :
        print " - %s : %.2f" % (tp, qs[tp])
    print ""


    if 1 :
        at = 30
        fp = None
        if os.path.isdir("/Users/greg/Dropbox/_mapsq") :
            fp = open ( "/Users/greg/Dropbox/_mapsq/scores%d_Q_allc_%s_sig%.0f.txt" % (at, label, sigma*100.0), "a" )
        elif os.path.isdir("/home/greg/Dropbox/_mapsq") :
            fp = open ( "/home/greg/Dropbox/_mapsq/scores%d_Q_allc_%s_sig%.0f.txt" % (at, label, sigma*100.0), "a" )
        elif os.path.isdir("C:/Users/greg/Dropbox/_mapsq") :
            fp = open ( "C:/Users/greg/Dropbox/_mapsq/scores%d_Q_allc_%s_sig%.0f.txt" % (at, label, sigma*100.0), "a" )
        else :
            fp = open ( "scores%d_Q_allc_%s_sig%.0f.txt" % (at, label, sigma*100.0), "a" )

        fp.write ( "%s\t%s\t%s" % (dmap.name, mol.name, res)  )

        for tp in tps :
            fp.write ( "\t%s\t%.2f" % (tp, qs[tp])  )

        fp.write ( "\n" )

        #nProt = len ( [at for at in mol.atoms if at.residue.isProt == True] )
        #nNA = len ( [at for at in mol.atoms if at.residue.isNA == True] )
        #fp.write ( "%s\t%s\t%s\t%d\t%d\n" % (dmap.name, mol.name, res, nProt, nNA)  )

        fp.close ()




def emringer ( dmap, mol ) :

    print "----- %s ____________ EMRINGER ____________ %s -----" % (dmap.name, mol.name)

    cdir = os.getcwd()
    print " - now in: ", cdir

    #print " - splitting " + mol.openedAs[0]
    mpath, mname = os.path.split ( mol.openedAs[0] )
    dpath, dname = os.path.split ( dmap.data.path )

    bs = os.path.splitext ( mol.openedAs[0] )[0]


    print " - copying mol file... removes symmetry/connect stuff"
    fin = open ( mol.openedAs[0], "r" )
    fout = open ( bs + "_.pdb", "w" )
    for line in fin :
        if "ATOM" in line or "HETATM" in line :
            fout.write ( line )
    fin.close ()
    fout.close ()


    phPath = "/Users/greg/_mol/phenix-1.14-3260/build/bin/"
    #phPath = "/Users/greg/_mol/phenix-1.15rc3-3435/build/bin/"

    args = [phPath+'phenix.emringer', dmap.data.path, bs+"_.pdb" ]
    print "running: ",
    for arg in args : print arg,
    print ""

    outf = mpath + '/' + '_out.txt'
    errf = mpath + '/' + '_err.txt'
    fout = open ( outf, "w" )
    ferr = open ( errf, "w" )
    import subprocess
    p = subprocess.Popen(args, stdout=fout, stderr=ferr, cwd=mpath)
    p.wait()
    fout.close()
    ferr.close()

    print " - getting score from " + outf
    score = -100
    fin = open ( outf )
    for l in fin :
        if "EMRinger Score:" in l :
            s = l [ len("EMRinger Score:")+1 : ]
            print "Score: ", s
            score = float( s )
            print " - found score: %.3f" % score

    print " - removing ", bs + "_.pdb"
    import shutil
    try :
        os.remove ( bs + "_.pdb" )
        os.remove ( bs + "__emringer.pkl" )
        os.remove ( bs + "__emringer.csv" )
        shutil.rmtree ( bs + "__emringer_plots" )
        print " - done"
    except :
        print "  -- did not find"

    return score


def refine ( dmap, mol, res ) :

    print "----- %s ____________ REFINE ____________ %s -----" % (dmap.name, mol.name)

    cdir = os.getcwd()
    print " - now in: ", cdir

    #print " - splitting " + mol.openedAs[0]
    mpath, mname = os.path.split ( mol.openedAs[0] )
    dpath, dname = os.path.split ( dmap.data.path )

    bs = os.path.splitext ( mol.openedAs[0] )[0]


    print " - copying mol file... removes symmetry/connect stuff"
    fin = open ( mol.openedAs[0], "r" )
    fout = open ( bs + "_.pdb", "w" )
    for line in fin :
        if "ATOM" in line or "HETATM" in line :
            fout.write ( line )
    fin.close ()
    fout.close ()


    phPath = "/Users/greg/_mol/phenix-1.14-3260/build/bin/"
    phPath = "/Users/greg/_mol/phenix-1.15rc3-3435/build/bin/"

    args = [phPath+'phenix.real_space_refine', dmap.data.path, bs+"_.pdb", "resolution=%.1f"%res ]
    print "running: ",
    for arg in args : print arg,
    print ""

    outf = mpath + '/' + '_out.txt'
    errf = mpath + '/' + '_err.txt'
    fout = open ( outf, "w" )
    ferr = open ( errf, "w" )
    import subprocess
    p = subprocess.Popen(args, stdout=fout, stderr=ferr, cwd=mpath)
    p.wait()
    fout.close()
    ferr.close()

    print " - getting score from " + outf
    score = -100
    fin = open ( outf )
    for l in fin :
        if "EMRinger Score:" in l :
            s = l [ len("EMRinger Score:")+1 : ]
            print "Score: ", s
            score = float( s )
            print " - found score: %.3f" % score

    print " - removing ", bs + "_.pdb"
    import shutil
    try :
        os.remove ( bs + "_.pdb" )
        os.remove ( bs + "__emringer.pkl" )
        os.remove ( bs + "__emringer.csv" )
        shutil.rmtree ( bs + "__emringer_plots" )
        print " - done"
    except :
        print "  -- did not find"

    return score


def refdir ( rdir ) :

    print "Refining in", rdir



def CalcR_ ( label = "" ) :

    print "Calc all scores -", label

    from VolumeViewer import Volume
    dmap = chimera.openModels.list(modelTypes = [Volume])[0]
    print " - dmap: %s" % dmap.name

    #fp = open ( "/Users/greg/_data/_mapsq/scores.txt", "a" )
    #fp.write ( "%s...\n" % dmap.name.split("_")[0]  )
    #fp.close ()

    from chimera import Molecule
    mol = chimera.openModels.list(modelTypes = [Molecule])[0]
    print " - mol: %s" % mol.name
    SetBBAts ( mol )


    mapName = os.path.splitext(dmap.name)[0]
    molName = os.path.splitext(mol.name)[0]
    ddir, dfile = os.path.split(dmap.data.path)

    molFile = mol.openedAs[0]
    mdir, mfile = os.path.split(molFile)

    print "PhFmap -- " + molFile

    RES = 3.0
    print " -- res %.1f -- " % RES

    outFile = molFile + "_r%.0f" % RES + "_fmodel.ccp4"

    if not os.path.isfile ( outFile ) :

        phPath = "/usr/local/phenix-1.14-3260/build/bin/"

        args = [phPath+'phenix.fmodel', "high_resolution=%.1f"%RES, "scattering_table=electron", "generate_fake_p1_symmetry=True", molFile ]
        print "running: ",
        for arg in args : print arg,
        print ""

        fout = open ( mdir + '/' + '_0_fmodel.log', "w" )
        import subprocess
        p = subprocess.Popen(args, stdout=fout, cwd=mdir)
        p.wait()
        fout.close()

        print ""
        args = [phPath+'phenix.mtz2map', "high_resolution=%.1f"%RES, "include_fmodel=true", "scattering_table=electron", molFile, molFile + ".mtz" ]
        print "running: ",
        for arg in args : print arg,
        print ""

        fout = open ( mdir + '/' + '_1_mtz2map.log', "w" )
        p = subprocess.Popen(args, stdout=fout, cwd=mdir)
        p.wait()
        fout.close()

        print " - renaming to:", outFile
        os.rename( molFile + "_fmodel.ccp4", outFile )
        os.remove( molFile + ".mtz" )


    print " - loading map:", outFile
    dm = chimera.openModels.open ( outFile )[0]



    molg = MyMolMapX ( mol, mol.atoms, RES, dmap.data.step[0], chimera.Xform.identity() )
    fpoints, fpoint_weights = fit_points_g ( molg, 0.1 )
    map_values = dmap.interpolated_values ( fpoints, mol.openState.xform )

    mmolap, mmcorr1, mmcorr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    print "Molmap - olap: %f, CC: %f, CCm: %f" % (mmolap, mmcorr1, mmcorr2)

    fpoints, fpoint_weights = fit_points_g ( dm.data, 5.0 )
    map_values = dmap.interpolated_values ( fpoints, dm.openState.xform )
    olap, phcorr1, phcorr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    print "Phenix - olap: %f, CC: %f, CCm: %f" % (olap, phcorr1, phcorr2)

    #fpoints, fpoint_weights = fit_points_g ( dmap.data, -1e6 )
    #map_values = dm.interpolated_values ( fpoints, dmap.openState.xform )
    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print "Phenix - olap: %f, CC: %f, CCm: %f" % (olap, corr1, corr2)


    print "%f\t%f\t%f\t%f" % (mmcorr1, mmcorr2, phcorr1, phcorr2)

    fp = open ( "/Users/greg/Dropbox/_mapsq/scores3_R_%s.txt" % label, "a" )
    fp.write ( "%s\t%f\t%f\t%f\t%f\n" % (dmap.name.split("_")[0], mmcorr1, mmcorr2, phcorr1, phcorr2)  )
    fp.close ()





def MaskMapResize ( atoms, R, dmap, fout=None ) :


    import _multiscale
    import _contour
    import _volume
    from _contour import affine_transform_vertices as transform_vertices
    from VolumeData import grid_indices, zone_masked_grid_data, interpolate_volume_data

    points = _multiscale.get_atom_coordinates ( atoms, transformed = True )
    #print " %d points" % len(points)
    fpoints = points


    if 0 :
        _contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, R )

        #mdata = VolumeData.Array_Grid_Data ( mdata.full_matrix(), dmap.data.origin, dmap.data.step, dmap.data.cell_angles, name = "atom masked" )


        mat = mdata.full_matrix()
        threshold = 1e-3

        points = _volume.high_indices(mat, threshold)
        fpoints = points.astype(numpy.single)
        fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

        #print " %d points" % len(points)


        nz = numpy.nonzero( fpoint_weights )[0]
        #print " %d pts nonzero" % len(nz)
        if len(nz) > 0 and len(nz) < len (fpoint_weights) :
            fpoints = numpy.take( fpoints, nz, axis=0 )
            fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    else :
        _contour.affine_transform_vertices ( fpoints, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        #transform_vertices( fpoints, dmap.data.ijk_to_xyz_transform )
        transform_vertices( fpoints, dmap.data.xyz_to_ijk_transform )

    #print " - %s mask %d atoms, %d nonzero points" % ( dmap.name, len(atoms), len(nz) )

    #transform_vertices( fpoints,  Matrix.xform_matrix( fmap.openState.xform ) )
    #transform_vertices( fpoints,  Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #transform_vertices ( fpoints, dmap.data.xyz_to_ijk_transform )

    bound = 6
    li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li, lj, lk, hi, hj, hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    #print " - new map origin:", nO

    ox = round ( nO[0]/dmap.data.step[0] ) * dmap.data.step[0]
    oy = round ( nO[1]/dmap.data.step[1] ) * dmap.data.step[1]
    oz = round ( nO[2]/dmap.data.step[2] ) * dmap.data.step[2]

    nO = ( ox, oy, oz )

    #print " - new map origin:", nO


    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    npoints = grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    if fout == None :
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
        dmap_base = os.path.splitext(dmap.name)[0]
        dmap_path = os.path.splitext (dmap.data.path)[0]
        nv.name = dmap_base + "_masked"
        nv.openState.xform = dmap.openState.xform
        return nv

    else :

        from VolumeData import save_grid_data
        #d = self.grid_data()
        format = save_grid_data(ndata, fout, None, {}, False)
        #print " - saved data"




def SetBBAts ( mol ) :

    #if hasattr ( mol, "bbats" ) :
    #    return
    #mol.bbats = True

    print " - setting bbAts in %s" % mol.name
    for r in mol.residues :

        #r.isProt = "C" in r.atomsMap and "CA" in r.atomsMap and "N" in r.atomsMap
        #r.isProt = "CA" in r.atomsMap
        #r.isNA = "O3'" in r.atomsMap and "O5'" in r.atomsMap

        from chimera.resCode import nucleic3to1
        from chimera.resCode import protein3to1
        protein3to1['HSD'] = protein3to1['HIS']
        protein3to1['HSE'] = protein3to1['HIS']

        r.isProt = r.type in protein3to1
        r.isNA = r.type in nucleic3to1

        r.score1 = None
        r.score2 = None

        if r.isProt :
            r.rtype = "prot"
        elif r.isNA :
            r.rtype = "na"
        else :
            r.rtype = "?"


        if r.isNA :
            try :
                if nucleic3to1[r.type] == "G" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "C" :
                    r.baseAt = r.atomsMap["N1"][0]
                elif nucleic3to1[r.type] == "A" :
                    r.baseAt = r.atomsMap["N9"][0]
                elif nucleic3to1[r.type] == "U" :
                    r.baseAt = r.atomsMap["N1"][0]
            except :
                #print " - baseAt not found - "
                pass


        r.bbAtoms = []
        r.scAtoms = []

        if r.isProt :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue
                n = a.name
                a.isBB = n=="C" or n=="CA" or n=="O" or n=="N" or n=="OT1" or n=="OT2"
                a.isSC = not a.isBB
                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

                a.isSugar, a.isBase = False, False

        elif r.isNA :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC = False, False
                    continue

                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="O4'" or n=="O2'" or n=="C3'" or n=="C4'"
                a.isBB = a.isBB or a.isSugar

                a.isBase = not a.isBB

                if nucleic3to1[r.type] == "G" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="C6" or n=="O6" or n=="N1" or n=="C2" or n=="N2" or n=="N3"

                elif nucleic3to1[r.type] == "C" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="N4" or n=="C5" or n=="C6"

                elif nucleic3to1[r.type] == "A" :
                    a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="N3" or n=="C2" or n=="N1" or n=="C6" or n=="N6"

                elif nucleic3to1[r.type] == "U" :
                    a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="O4" or n=="C5" or n=="C6"

                else :
                    #print " -x- NA res %d.%s is ?" % (r.id.position, r.type)
                    break

                a.isSC = a.isBase

                #if nucleic3to1[r.type] == "G" :
                #    r.isBase = n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n=="" or n="" or n="" or n=""
                #    r.baseAt = r.atomsMap["N9"][0]

                if a.isBB :
                    r.bbAtoms.append ( a )
                else :
                    r.scAtoms.append ( a )

        else :
            for a in r.atoms :
                a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False
