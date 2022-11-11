

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
import VolumeViewer
import _contour
import _gaussian



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





# attempt to do Q-score with volume-volume CC rather than sphere points
# works ok, but is not faster - main reason why to try
# another difference is that with sphere points, the same number of points
# is used at each distance, so the map values at each radial distance even weigh

def QscoreM ( atoms, dmap, sigma, agrid=None, allAtTree=None, show=0, log=0, toRAD=2.0, step=0.2, minD=None, maxD=None, useMask=False ) :

    xyz = _multiscale.get_atom_coordinates(atoms, transformed = False)

    #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )

    li,lj,lk = numpy.min ( xyz, axis=0 ) - (toRAD, toRAD, toRAD)
    hi,hj,hk = numpy.max ( xyz, axis=0 ) + (toRAD, toRAD, toRAD)
    nO = ( li, lj, lk )
    #print nO
    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, d1,d2,d3 )

    d1, d2, d3 = hi - li, hj - lj, hk - lk
    nstep = (step, step, step)
    #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

    nn1 = int ( numpy.ceil ( float(d1) / step) )
    nn2 = int ( numpy.ceil ( float(d2) / step) )
    nn3 = int ( numpy.ceil ( float(d3) / step) )

    #print " - step %.2f, n: %d %d %d" % (S, nn1, nn2, nn3)

    nmat = numpy.zeros ( (nn3,nn2,nn1), numpy.float32 )

    ii = 1.0 / step
    ni = -ii
    xyz_to_ijk = ((ii, 0.0, 0.0, ni*nO[0]), (0.0, ii, 0.0, ni*nO[1]), (0.0, 0.0, ii, ni*nO[2]))
    ijk_to_xyz = ((step, 0.0, 0.0, nO[0]), (0.0, step, 0.0, nO[1]), (0.0, 0.0, step, nO[2]))

    #print ijk_to_xyz

    #ijk[:] = xyz
    weights = [ 1.0 for a in atoms]
    sdevs = [ [sigma, sigma, sigma] for a in atoms ]
    cutoff_range = 5

    A, B = maxD - minD, minD


    #ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
    #print ndata.xyz_to_ijk_transform
    #print ndata.ijk_to_xyz_transform
    #Matrix.transform_points(xyz, ndata.xyz_to_ijk_transform)

    if useMask == False :
        Matrix.transform_points(xyz, xyz_to_ijk)
        _gaussian.sum_of_gaussians(xyz, weights, sdevs, cutoff_range, nmat)

        #print " -gm max %.3f" % numpy.max ( nmat )
        nmat *= A
        nmat += B
        #print " -gm max %.3f" % numpy.max ( nmat )


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


    if useMask :

        nearAts = []
        if agrid != None :
            for at in atoms :
                nats = agrid.AtsNearPtLocal ( at.coord() )
                for nat, v in nats :
                    if at != nat :
                        nearAts.append ( nat )
                        #print " - %s, %d.%s - %.3f" % (nat.name, nat.residue.id.position, nat.residue.id.chainId, v.length)

        if allAtTree != None :
            for at in atoms :
                opointsNear = allAtTree.searchTree ( at.coord(), toRAD )
                for nat in opointsNear :
                    if nat == at :
                        continue
                    v = at.coord() - nat.coord()
                    if v.length < toRAD :
                        nearAts.append (nat)

        if len(nearAts) == 0 :
            print " - no near ats?"

        #print " - %d near ats" % len(nearAts)


        for k in range(nn3) :
            pz = nO[2] + float(k)*step

            for j in range(nn2) :
                py = nO[1] + float(j)*step

                for i in range(nn1) :
                    px = nO[0] + float(i)*step

                    P = chimera.Point(px, py, pz)

                    minDToAt = 1e9
                    for at in atoms :
                        v = at.coord() - P
                        if v.length < minDToAt :
                            minDToAt = v.length

                    if minDToAt > toRAD :
                        nmat[k,j,i] = B-0.1
                        continue

                    closestToAt = True
                    for nat in nearAts :
                        v = nat.coord() - P
                        if v.length < minDToAt :
                            closestToAt = False
                            #break

                    if not closestToAt :
                        nmat[k,j,i] = minD-0.1
                    else :
                        nmat[k,j,i] = A * numpy.exp ( -0.5 * numpy.power(minDToAt/sigma,2) ) + B



    if 0 and agrid :
        nearAts = []
        for at in atoms :
            nats = agrid.AtsNearPtLocal ( at.coord() )
            for nat, v in nats :
                if at != nat :
                    print " - %s, %d.%s - %.3f" % (nat.name, nat.residue.id.position, nat.residue.id.chainId, v.length)
                    nearAts.append ( at )

        #print "%d near ats" % len(nearAts)
        mat1 = numpy.ones ( (nn1,nn2,nn3), numpy.float32 )
        ndata = VolumeData.Array_Grid_Data ( mat1, nO, nstep, dmap.data.cell_angles )
        points = _multiscale.get_atom_coordinates(nearAts, transformed = False)

        mdata = VolumeData.zone_masked_grid_data ( ndata, points, toRAD, invert_mask=False )
        #nmat = mdata.matrix()
        nv = VolumeViewer.volume.volume_from_grid_data ( mdata )
        nv.openState.xform = dmap.openState.xform
        mdata = mask


    fpoints = VolumeData.grid_indices ( (nn3,nn2,nn1), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( fpoints, ijk_to_xyz )
    fpoint_weights = numpy.ravel(nmat).astype(numpy.single)


    #print " - %d points" % len(fpoints)
    ge = numpy.greater_equal(fpoint_weights, B)
    fpoints = numpy.compress(ge, fpoints, 0)
    fpoint_weights = numpy.compress(ge, fpoint_weights)
    #print " - %d above thr" % len(fpoint_weights)
    #nz = numpy.nonzero( fpoint_weights )[0]
    #print " - %d above thr" % len(nz)

    #map_values, outside = VolumeData.interpolate_volume_data(pts, xyz_to_ijk_tf, darray)
    #olap0, cc0, other = overlap_and_correlation ( wts, map_values )

    map_values = dmap.interpolated_values ( fpoints, atoms[0].molecule.openState.xform )
    #print map_values
    olap, cc, ccm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print olap, cc, ccm


    if show :
        ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )
        nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
        nv.openState.xform = dmap.openState.xform
        nv.name = "bam"




    return ccm



def zone_mask ( grid_data, zone_points, zone_radius, invert_mask = False, zone_point_mask_values = None ):

    from numpy import single as floatc, array, ndarray, zeros, int8, intc

    if not isinstance(zone_points, ndarray):
        zone_points = array(zone_points, floatc)

    if (not zone_point_mask_values is None and not isinstance(zone_point_mask_values, ndarray)):
        zone_point_mask_values = array(zone_point_mask_values, int8)

    shape = tuple(reversed(grid_data.size))
    mask_3d = zeros(shape, int8)
    mask_1d = mask_3d.ravel()

    if zone_point_mask_values is None:
        if invert_mask:
            mask_value = 0
            mask_1d[:] = 1
        else:
            mask_value = 1

    from VolumeData import grid_indices
    from _contour import affine_transform_vertices
    from _closepoints import find_closest_points, BOXES_METHOD

    size_limit = 2 ** 22          # 4 Mvoxels
    if mask_3d.size > size_limit:
        # Calculate plane by plane to save memory with grid point array
        xsize, ysize, zsize = grid_data.size
        grid_points = grid_indices((xsize,ysize,1), floatc)
        affine_transform_vertices(grid_points, grid_data.ijk_to_xyz_transform)
        zstep = [grid_data.ijk_to_xyz_transform[a][2] for a in range(3)]
        for z in range(zsize):
            i1, i2, n1 = find_closest_points(BOXES_METHOD, grid_points, zone_points, zone_radius)
            offset = xsize*ysize*z
            if zone_point_mask_values is None:
                mask_1d[i1 + offset] = mask_value
            else:
                mask_1d[i1 + offset] = zone_point_mask_values[n1]
            grid_points[:,:] += zstep
    else :
        grid_points = grid_indices(grid_data.size, floatc)
        affine_transform_vertices(grid_points, grid_data.ijk_to_xyz_transform)
        i1, i2, n1 = find_closest_points(BOXES_METHOD, grid_points, zone_points, zone_radius)
        if zone_point_mask_values is None:
            mask_1d[i1] = mask_value
        else:
            mask_1d[i1] = zone_point_mask_values[n1]

    return mask_3d




def sseQscores ( mol, dmap, sigma ) :

    print ""
    print "SSE Q-scores"
    print " - sigma: %.2f" % sigma
    print " - in map: %s" % dmap.name
    print " - for mol: %s" % mol.name

    minD, maxD = MinMaxD (dmap)

    cress = {}
    for r in mol.residues :
        if r.id.chainId in cress :
            cress[r.id.chainId].append ( [r.id.position, r] )
        else :
            cress[r.id.chainId] = [ [r.id.position, r] ]

    aqs = []
    for chainId, ress in cress.iteritems() :
        ress.sort()
        ress = [r for ri, r in ress]
        sses = SSEs ( ress )
        print " - chain %s, %d ress, %d sses" % (chainId, len(ress), len(sses)),

        aqsc = []
        for startRi, endRi, ss, ress in sses :
            if ss == "H" :
                #print "a-%d" % startRi,
                aQ = sseQscore ( ress, "H", dmap, sigma, showRes=None, log=0, numPts=6, toRAD=3.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0, mol=mol )
                if aQ == None :
                    continue
                aqs.append ( aQ )
                aqsc.append ( aQ )
            else :
                #print "%d" % startRi,
                pass

        #print ""
        if len(aqsc) > 0 :
            aqsc.sort ( reverse=True )
            print "%.2f, %.2f" % (numpy.average(aqsc), aqsc[0])
        else :
            print " - "

    avg, N = 0.0, 0.0
    for r in mol.residues :
        if hasattr ( r, 'aQ' ) :
            avg += r.aQ
            N += 1.0

    if N > 0.0 :
        avg = avg / N

    avgMax = 0.0
    if len(aqs) > 0 :
        aqs.sort( reverse = True )
        avgMax = numpy.average ( aqs[0:5] )

    print " - avg alpha Q: %.2f, top 5 avgs %.2f" % (avg, avgMax)

    return avg, avgMax



# this method calculates CC between radial points placed around sse

def sseQscore ( ress, sseType, dmap, sigma, showRes=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    if len(ress) == 0 :
        return None

    sseQ = 0.0
    mol = ress[0].atoms[0].molecule

    if sseType == "H" :
        if log :
            print " - helix - %d residues" % len(ress)

        sseQ = 0.0

        if len(ress) < 4 :
            # too small
            return None

        pts = []
        ptsRes = []
        bbAts = []
        avgN = 3
        showResI = -1
        for ri, r in enumerate ( ress[:-avgN+1] ) :
            avgPt = numpy.array ( [0,0,0] )
            P = chimera.Point(0,0,0)
            N = 0.0
            for i in range (avgN) :
                rr = ress[ri+i]
                #avgPt += rr.atomsMap["CA"][0].coord().toVector(); N += 1.0
                #avgPt += rr.atomsMap["C"][0].coord().toVector(); N += 1.0
                #avgPt += rr.atomsMap["N"][0].coord().toVector(); N += 1.0
                P += rr.atomsMap["CA"][0].coord().toVector(); N += 1.0
                #P += rr.atomsMap["C"][0].coord().toVector(); N += 1.0
                #P += rr.atomsMap["N"][0].coord().toVector(); N += 1.0

            #avgPt /= N
            #pts.append ( avgPt )
            p = P.toVector() / N
            pts.append ( chimera.Point (*p) )
            ptsRes.append ( ress[ri+1] )
            if showRes != None and showRes == ress[ri+1] :
                print " - found show res at %d" % ri
                showResI = ri

        if log :
            print "  - %d central points " % len(ptsRes)

        for i in range (3) :
            for pi in range ( 1, len(pts)-1 ) :
                p0, p1, p2 = pts[pi-1], pts[pi], pts[pi+1]
                v = p2 - p0; v.normalize()
                v1 = p1 - p0
                d = v1 * v
                mp = p0 + d * v
                pts[pi] = mp

        if showRes != None :
            mod = GetMod ( "helix ctr points" )
            if mod : chimera.openModels.remove ( mod )
            pmod = AddSphereTubesPts ( pts, (.9,.3,.3,1.0), 0.5, "helix ctr points" )
            pmod.openState.xform = ress[0].atoms[0].molecule.openState.xform


        avgQ, N = 0.0, 0.0
        for i, pt in enumerate ( pts ) :

            res = ptsRes[i]

            p1, p2, p3 = None, None, None
            if i-1 >= 0 and i+1 < len(pts) :
                p1, p2, p3 = pts[i-1:i+2]
            elif i-1 < 0 :
                p2, p3 = pts[i:i+2]; p1 = p2
            else :
                p1, p2 = pts[i-1:i+1]; p3 = p2

            v0 = res.atomsMap["CA"][0].coord() - p2

            show = (i == showResI)
            if show :
                aQ = alphaQ ( p2, p3-p1, v0, dmap, sigma=sigma, show=show, log=0, numPts=6, toRAD=3.0, dRAD=0.5, minD=minD, maxD=maxD, fitg=0, mol=mol )[0]
                print " - show Q: %.2f" % aQ

                if log :
                    aQ = alphaQ ( p2, p3-p1, v0, dmap, sigma=sigma, show=0, log=log, numPts=numPts, toRAD=toRAD, dRAD=dRAD, minD=minD, maxD=maxD, fitg=1, mol=mol )[0]
                    print " - %d, res %d - log Q: %.2f" % (i, res.id.position, aQ)

            aQ = alphaQ ( p2, p3-p1, v0, dmap, sigma=sigma, show=0, log=0, numPts=numPts, toRAD=toRAD, dRAD=dRAD, minD=minD, maxD=maxD, fitg=0, mol=mol )[0]
            if log :
                print " - %d, res %d - Q: %.2f" % (i, res.id.position, aQ)
            res.aQ = aQ
            avgQ += aQ
            N += 1.0

        avgQ /= N
        return avgQ





def alphaQ ( p2, vA, v0, dmap, sigma=3.0, show=0, log=0, numPts=8, toRAD=3.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :

    # calculates Q-score for alpha helix
    vA.normalize()
    v1 = None
    for i in range ( 100 ) :
        vR = v0 + chimera.Vector ( i*.1, 0, 0 )
        v1 = chimera.cross ( vA, vR )
        if v1.length > 0.01 :
            v1.normalize()
            break

    #from _multiscale import get_atom_coordinates
    #pts = get_atom_coordinates(atoms, transformed = False)
    #print " __%s__ " % (atoms[0].name), pts[0]

    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals will hold values taken from the reference gaussian...
    g_vals = (numpy.ones ( [numPts,1] ) * refG).astype(numpy.float64, copy=False)

    # d_vals holds the average map values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( numpy.array ( [p2] ), mol.openState.xform ).astype(numpy.float64, copy=False)
    d_vals = numpy.repeat ( d_vals, numPts )
    avgV = numpy.average ( d_vals )

    # r_avg holds the average values and number of points at each radial distance
    r_avg = [ [0,avgV,numPts] ]

    RAD = dRAD
    while RAD < toRAD + 0.01 :

        radPts = []

        np = int ( numpy.ceil ( numPts * (RAD / dRAD) ) ) if show else numPts

        for i in range ( np ) :
            ang = float(i) * 360.0 / float(np)
            xf = chimera.Xform.rotation ( vA, ang )
            p = p2 + xf.apply (v1) * RAD
            radPts.append ( p )

        if show :
            mod = GetMod ( "helix radial points %.1f" % (RAD) )
            if mod : chimera.openModels.remove ( mod )
            pmod = AddSpherePts ( radPts, (.3,.9,.3,1.0), 0.5, "helix radial points %.1f" % (RAD) )
            pmod.openState.xform = mol.openState.xform

        d_vals_n = dmap.interpolated_values ( radPts, mol.openState.xform )
        d_vals = numpy.append ( d_vals, d_vals_n )
        avg = numpy.average ( d_vals_n )

        A,B = maxD - minD, minD
        gv = A * numpy.exp ( -0.5 * numpy.power(RAD/sigma,2) ) + B
        g_vals = numpy.append ( g_vals, numpy.ones([len(radPts),1]) * gv )

        r_avg.append ( [RAD,avg,len(radPts)] )

        RAD += dRAD


    if log and not fitg :
        # log means output the radial profile
        # fitg means fit gaussian and output profile and fit gaussian
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


    olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals, g_vals )
    # this is the CC between _all_ radial values
    Qs = CCmean
    if log :
        print ""
        print "alpha Q-score - sigma %.2f" % (sigma)
        print " - olap: %.3f cc: %.3f, ccmean: [%.3f], %d values" % (olap, CC, Qs, len(d_vals))
        print ""
        #print "%f\t%f\t%f" % (olap, CC, Qs)


    if fitg :
        if log : print "fitting gaussian : "
        V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        #V, N = [ [x[0],x[1]] for x in r_avg[0:15] ], float(15)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        errp = err / r_avg[0][1] * 100.0
        if log :
            print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )

        if 0 and sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            #print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log :
                print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        mx = 0.0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if y > mx :
                mx = y
            if gv > mx :
                mx = gv
            #if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        #print ""

        # normalized values and fitted gaussian
        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y/mx, gv/mx)
            i += 1

        return Qs, yds, err

    else :
        return [Qs]





# this method calculates CC between radial points placed around the atoms and the map
# - two values are returned - CC and CC about the mean - the latter is the Q-score


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

    # g_vals will hold values taken from the reference gaussian...
    g_vals = (numpy.ones ( [len(pts)*numPts,1] ) * refG).astype(numpy.float64, copy=False)
    #g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False)

    if mol == None :
        mol = atoms[0].molecule

    # d_vals holds the average map values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( pts, mol.openState.xform ).astype(numpy.float64, copy=False)
    #print pts
    #print d_vals
    d_vals = numpy.repeat ( d_vals, numPts )

    # r_avg holds the average values and number of points at each radial distance
    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,len(pts)*numPts] ]
    #d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)


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

    # check if any atoms are too close; ignore those atoms and give them q=0
    if 0 and allAtTree :
        for at in atoms :
            anear = allAtTree.searchTree ( at.coord().data(), 2.0 )
            for nat in anear :
                if nat != at :
                    v = at.coord() - nat.coord()
                    if v.length < 1.0 :
                        print "c"
                        return 0.0


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

            # this loop repeats until we get ~numPts at this radial distance or give up after 50 tries
            # at each try, the number of points distributed on a sphere with radius=RAD is increased by 2
            # from these points, only the ones closer to the atom 'at' than any other atom are kept
            # this way it tries to get evenly distributed points
            for i in range (0, 50) : # --here
                outPts = SpherePts ( at.coord(), RAD, npts+i*2 )
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
                #print " - %d, %d pts" % (i, len(at_pts))
                if at_pts_i >= npts or show :
                    #print " - %.2f - after %d" % (RAD, i)
                    pts.extend ( at_pts[0:at_pts_i] )
                    break

        if show :
            pmod = AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f %s" % (RAD,atoms[0].name) )
            pmod.openState.xform = atoms[0].molecule.openState.xform

        #print " - rad %.1f - %d pts" % ( RAD, len(pts) )

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

            #g_vals_avg = numpy.append ( g_vals_avg, gv )
            #d_vals_avg = numpy.append ( d_vals_avg, avg )

            r_avg.append ( [RAD,avg,len(pts)] )

            #if log :
            #    print "%.1f\t%f\t%f\t%d" % (RAD, avg, gv, len(pts))

        RAD += dRAD
        i+=1

    if log and not fitg :
        # log here means output the radial profile
        # fitg means fit gaussian and output profile and fit gaussian
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

    # this is the CC between averaged radial values - not as robust
    #olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals_avg, g_vals_avg )
    #if log :
    #    print "olap -avg-: %.3f cc: %.3f, Q: %.3f -- %d" % (olap, CC, Qs, len(d_vals_avg))
    #    #print "%f\t%f\t%f" % (olap, CC, Qs)

    olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals, g_vals )
    # this is the CC between _all_ radial values
    Qs = CCmean
    if log :
        print "olap --N--: %.3f cc: %.3f, ccmean (Q-score): %.3f -- %d" % (olap, CC, Qs, len(d_vals))
        #print "%f\t%f\t%f" % (olap, CC, Qs)


    if fitg :
        if log : print "fitting gaussian : "
        V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        #V, N = [ [x[0],x[1]] for x in r_avg[0:15] ], float(15)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        errp = err / r_avg[0][1] * 100.0
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            #print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        mx = 0.0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if y > mx :
                mx = y
            if gv > mx :
                mx = gv
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        print ""

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y/mx, gv/mx)
            i += 1

        return Qs, yds, err

    else :
        return Qs




# qscore using grid instead of atomsTree - cuts time by half

def QscoreG ( atoms, dmap, sigma, agrid=None, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0, mol=None ) :

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

    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        # the 0.9 factor is a small adjustment  so that the atom in question is not found as a 'nearby' atom
        # this should not affect the final score significantly
        outRad = RAD*0.9
        #outRad2 = outRad * outRad
        pts = []
        for at in atoms :
            #npts = numPts # 8 # int ( npts )
            npts = numPts
            #if show : npts = int (numPts * RAD*RAD / (dRAD*dRAD))
            #npts = numPts * (RAD*RAD / (dRAD*dRAD))
            #print RAD, dRAD, numPts, " -> ", npts
            for i in range (0, 50) : # --here
                outPts = SpherePts ( at.coord(), RAD, npts+i*2 )
                at_pts, at_pts_i = [None]*len(outPts), 0
                for pt in outPts :
                    #vPt = [pt[0], pt[1], pt[2]]
                    #apt = numpy.array ( vPt )
                    #P = chimera.Point ( pt[0], pt[1], pt[2] )
                    if agrid != None :
                        #opointsNear = allAtTree.searchTree ( vPt, outRad )
                        nearAts = agrid.AtsNearPtLocal ( pt, outRad )
                        #numNearAts = agrid.NumAtsNearPtLocal ( pt, outRad )
                        #if agrid.NumAtsNearAtLocal(at,D=outRad) < 1 :
                        #if numNearAts == 0 :
                        if len(nearAts) <= 0 :
                            #print ",",
                            at_pts[at_pts_i] = pt
                            at_pts_i += 1
                    else :
                        at_pts[at_pts_i] = pt
                        at_pts_i += 1
                #if log :
                #print " - %d, %d pts" % (i, len(at_pts))
                if at_pts_i >= npts : # or show :
                    pts.extend ( at_pts[0:at_pts_i] )
                    break

        if show :
            pmod = AddSpherePts ( pts, (.6,.6,.6,0.4), 0.1, "RAD points %.1f %s" % (RAD,atoms[0].name) )
            pmod.openState.xform = atoms[0].molecule.openState.xform

        #print " - rad %.1f - %d pts" % ( RAD, len(pts) )

        if len (pts) < 1 :
            # when no points are found that are closer to the atom being considered than other atoms
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


    olap, CC, CCmean = FitMap.overlap_and_correlation ( d_vals, g_vals )
    Qscore = CCmean
    if log :
        print "olap --N--: %.3f cc: %.3f, ccmean (Q-score): %.3f -- %d" % (olap, CC, Qscore, len(d_vals))
        #print "%f\t%f\t%f" % (olap, CC, Qs)


    if fitg :
        if log : print "fitting gaussian : "
        #V, N = [ [x[0],x[1]] for x in r_avg ], float(len(r_avg))
        V, N = [ [x[0],x[1]] for x in r_avg[0:15] ], float(15)

        sdev, A, B = optSGD ( V, 5000, 1.0 )
        sdev, A, B = optSGD ( V, 5000, 0.1, sdev, A, B )
        err = numpy.sqrt(err3(V,sdev,A,B)/N)
        errp = err / r_avg[0][1] * 100.0
        if log : print " sgd - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)
        sdev2, A2, B2 = optGN ( V, 0.0001, sdev, A, B )
        if sdev2 != None :
            sdev, A, B = sdev2, A2, B2
            err = numpy.sqrt(err3(V,sdev,A,B)/N)
            #print "max:", r_avg[0][1]
            errp = err / r_avg[0][1] * 100.0
            if log : print "  gn - sdev: %.4f, A %.4f, B %.4f, err: %f (%.1f%%)" % (sdev, A, B, err, errp)

        yds, i = numpy.zeros ( len(r_avg) ), 0
        mx = 0.0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if y > mx :
                mx = y
            if gv > mx :
                mx = gv
            if log : print "%.1f\t%f\t%f" % (x, y, gv)
            i += 1

        print ""

        yds, i = numpy.zeros ( len(r_avg) ), 0
        for x, y, n in r_avg:
            gv = A * numpy.exp ( -0.5 * numpy.power(x/sdev,2) ) + B
            #yds[i] = y - gv
            yds[i] = y
            if log : print "%.1f\t%f\t%f" % (x, y/mx, gv/mx)
            i += 1

        return Qscore, yds, err

    else :
        return Qscore




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



# calculate Q-score given a point (rather than atom), and using a 'points tree' rather than 'atoms tree'

def QscorePt2 ( atPt, xfI, dmap, sigma, allPtTree = None, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0 ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    #xfI = chimera.Xform()
    atPtC = chimera.Point ( *atPt )
    #print atPtC

    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals should have the reference gaussian...
    g_vals = (numpy.ones ( [numPts,1] ) * refG).astype(numpy.float64, copy=False )
    g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False )

    # r_avg holds the average values and number of points at each radial distance
    d_vals = dmap.interpolated_values ( [atPt], xfI ).astype(numpy.float64, copy=False)
    #print atPt
    #print d_vals
    d_vals = numpy.repeat ( d_vals, numPts )

    avgV = numpy.average ( d_vals )
    r_avg = [ [0,avgV,numPts] ]

    d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)

    # make smaller atom tree, shaves a few ms off running time for each point
    if 1 and allPtTree != None :
        pts_near = []
        anear = allPtTree.searchTree ( atPt, toRAD*2.0 )
        pts_near.extend ( anear )

        #points = _multiscale.get_atom_coordinates ( ats_near, transformed = False )
        if log :
            print " - new search tree: %d pts" % ( len(ats_near) )
        allPtTree = AdaptiveTree ( pts_near, pts_near, 1.0)

    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.01 :
        outRad = RAD*0.9
        outRad2 = outRad * outRad
        #outRad2 = outRad * outRad
        pts = []

        # try to get at least [numPts] points at [RAD] distance
        # from the atom, that are not closer to other atoms
        for i in range (0, 50) :
            # points on a sphere at radius RAD...
            outPts = SpherePts ( atPtC, RAD, numPts+i*2 )
            at_pts, at_pts_i = [None]*len(outPts), 0
            for pt in outPts :
                vPt = [pt[0], pt[1], pt[2]]
                apt = numpy.array ( vPt )
                if allPtTree != None :
                    opointsNear = allPtTree.searchTree ( vPt, outRad )
                    foundNearPt = False
                    for npt in opointsNear :
                        v = apt - npt
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
            if at_pts_i >= numPts : # or show :
                #print " - %.2f - after %d" % (RAD, i)
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





# calculate Q-score given a point (rather than atom), and using a 'points grid' rather than 'atoms tree'
# (points grid is much faster than atoms tree)

def QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=None, log=0, numPts=8, toRAD=2.0, dRAD=0.5, minD=None, maxD=None, fitg=0 ) :

    if minD == None or maxD == None :
        minD, maxD = MinMaxD (dmap)

    A, B = maxD - minD, minD
    refG = A * numpy.exp ( -0.5 * numpy.power(0.0/sigma,2) ) + B
    #print " - refg: ", refG

    # g_vals will have the reference gaussian values
    g_vals = (numpy.ones ( [numPts,1] ) * refG).astype(numpy.float64, copy=False )
    # average values at each radial distance
    g_vals_avg = numpy.array ( [refG] ).astype(numpy.float64, copy=False )

    # d_vals will hold the map average values
    d_vals = dmap.interpolated_values ( [atPt], xfI ).astype(numpy.float64, copy=False)
    #print atPt
    #print d_vals
    d_vals = numpy.repeat ( d_vals, numPts )
    avgV = numpy.average ( d_vals )

    # r_avg holds the average values and number of points at each radial distance
    r_avg = [ [0,avgV,numPts] ]

    d_vals_avg = numpy.array ( [avgV] ).astype(numpy.float64, copy=False)

    #olap, corr1, corr2 = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #dRAD, toRAD, RAD = 0.2, 1.8, 0.1

    atPtC = chimera.Point ( *atPt )
    RAD = dRAD
    i = 1.0
    while RAD < toRAD + 0.001 :

        # this small adjustment makes sure the point being considered is not
        # found in the near-pts list for each point at a given radial distance
        outRad = RAD*0.9
        #outRad2 = outRad * outRad
        pts = []

        # try to get at least [numPts] points at [RAD] distance
        # from the atom, that are not closer to other atoms
        for i in range (0, 50) :
            # points on a sphere at radius RAD...
            outPts = SpherePts ( atPtC, RAD, numPts+i*2 )
            at_pts, at_pts_i = [None]*len(outPts), 0
            for pt in outPts :
                if ptGrid != None :
                    ptsNear = ptGrid.NumPtsNearPt ( pt, outRad )
                    if ptsNear == 0 :
                        at_pts[at_pts_i] = pt.data()
                        at_pts_i += 1
                else :
                    at_pts[at_pts_i] = pt.data()
                    at_pts_i += 1
            #if log :
            #    print " - %d, %d pts" % (i, len(at_pts))
            if at_pts_i >= numPts : # or show :
                #print " - %.2f - after %d" % (RAD, i)
                pts.extend ( at_pts[0:at_pts_i] )
                break

        if len (pts) < 1 :
            if log :
                print " - no points for RAD %.1f - " % RAD,
                #print "SC" if atoms[0].isSC else "BB"

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




def CalcQForAts ( dmap, mol, ats, sigma=0.6 ) :

    minD, maxD = MinMaxD (dmap)

    from _multiscale import get_atom_coordinates
    from CGLutil.AdaptiveTree import AdaptiveTree

    allAts = [at for at in mol.atoms if not at.element.name == "H"]
    points = get_atom_coordinates ( allAts, transformed = False )
    allAtTree = AdaptiveTree ( points.tolist(), allAts, 1.0)

    for at in ats :
        at.Q = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, minD=minD, maxD=maxD )



def Calc ( chimeraPath, mol, numProc, res=3.0, bfactorF=-1, sigma=0.6 ) :

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
    print " - model: %s" % mol.name

    #from chimera import Molecule
    #mols = chimera.openModels.list(modelTypes = [Molecule])
    #if len(mols) == 0 :
    #    print " - no molecules loaded"
    #    return

    SetBBAts ( mol )

    # may be closed by CalcQpn
    dmap_name = dmap.name

    if numProc == 1 :
        CalcQ ( mol, None, dmap, sigma, log=False )
    else :
        #CalcQp ( mol, None, dmap, sigma, numProc=numProc, chimeraPath=chimeraPath )
        CalcQpn ( mol, None, dmap, sigma, numProc=numProc, chimeraPath=chimeraPath, closeMap=True )

    SaveQStats ( mol, "All", dmap_name, sigma, res )

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
    D = chimera.openModels.list(modelTypes = [Volume])
    dmap = D[0]
    #dmapA = D[1]
    print " - dmap: %s" % dmap.data.path
    #print " - dmapA: %s" % dmapA.name

    tempPath, mapNameExt = os.path.split ( dmap.data.path )
    mapName, mapExt = os.path.splitext ( mapNameExt )
    procI = mapName.split("_")[0]
    print " - proc: %s" % procI
    print " - in path: %s" % tempPath

    # read ids and position of all atoms
    aPosMap = {}
    fina = open ( os.path.join(tempPath, "all_atoms.txt") )
    l1 = fina.readline()
    sigma, minD, maxD, numAts = l1.split()
    sigma, minD, maxD, numAts = float(sigma), float(minD), float(maxD), int(numAts)
    #fout.write ( "Sigma: %.1f, minD: %.3f, maxD: %.3f, numAllAts: %d\n" % (sigma, minD, maxD, numAts) )
    print "Sigma: %.1f, minD: %.3f, maxD: %.3f, numAllAts: %d\n" % (sigma, minD, maxD, numAts)

    allPts = [None] * numAts # numpy.array.zeros ( [numAts,3] )
    ati = 0
    for l in fina :
        atIdStr, sx, sy, sz = l.split()
        pt = [ float(sx), float(sy), float(sz) ]
        allPts[ati] = [ float(sx), float(sy), float(sz) ]
        aPosMap[atIdStr] = pt
        ati += 1
    fina.close()

    if ati != numAts :
        print " ---!!--- got only %d of %d atoms" (ati, numAts)
        return

    ##ats = [at for at in mol.atoms if not at.element.name == "H"]
    ##points = _multiscale.get_atom_coordinates ( allPts, transformed = False )
    ##print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )

    #allPtTree = AdaptiveTree ( allPts, allPts, 1.0)
    #print " - points tree with %d points" % len(allPts)

    #points = _multiscale.get_atom_coordinates ( ats, transformed = False )
    from gridm import Grid
    ptGrid = Grid()
    ptGrid.FromPoints ( allPts, 3.0 )
    print " - %d pts grid" % len(allPts)

    fin = open ( os.path.join ( tempPath, "%s_atoms.txt" % procI ) )
    fout = open ( os.path.join ( tempPath, "%s_out.txt" % procI ), "w" )
    fout_status = os.path.join ( tempPath, "%s_stat.txt" % procI )

    # get positions of atoms to do in this process
    doAts = []
    for l in fin :
        #atIdStr, sx, sy, sz = l.split()
        atIdStr = l.strip()
        if not atIdStr in aPosMap :
            fout.write (  "- atid not found: %s\n" % atIdStr )
        #at = atids[atIdStr]
        #pt = [ float(sx), float(sy), float(sz) ]
        doAts.append ( [atIdStr, aPosMap[atIdStr]] )
    fin.close()

    # write status to a temp file
    fs = open ( fout_status, "w" );
    fs.write ( "at atom %d/%d" % (0,len(doAts) ) );
    fs.close()

    import time
    start = time.time()
    xfI = dmap.openState.xform

    i = 1
    for atId, atPt in doAts :
        #print "%d.%s.%s" % (r.id.position,r.id.chainId,at.name),

        ##qs = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
        #qs = QscorePt2 ( atPt, xfI, dmap, sigma, allPtTree=allPtTree, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
        qs = QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=ptGrid, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

        fout.write ( "%s %f\n" % (atId, qs) )

        if i%10 == 0 :

            end = time.time()
            totSec = end - start

            leftTime = ""
            leftSec = 0.0
            iPerSec = float(i) / totSec
            if iPerSec > 0 :
                leftSec = float ( len(doAts) - i ) / iPerSec
                leftHour = numpy.floor ( leftSec / 60.0 / 60.0 )
                leftSec = leftSec - leftHour * 60.0 * 60.0
                leftMin = numpy.floor ( leftSec / 60.0 )
                leftSec = leftSec - leftMin * 60.0
                leftTime = "%.0f:%.0f:%.0f" % (leftHour, leftMin, leftSec)

            # update status
            fs = open ( fout_status, "w" );
            fs.write ( "at atom %d/%d - ETA %s" % (i+1,len(doAts),leftTime) );
            fs.close()

        i += 1

    fout.close()

    fs = open ( fout_status, "w" ); fs.write ( "done" ); fs.close()


def GetChiPath () :
    import sys
    # '/Users/greg/_mol/Chimera.app/Contents/Resources/share/__main__.py'
    chimeraPath = os.path.split ( sys.argv[0] )[0]
    print ""
    print " ------------ ", chimeraPath
    print ""

    chimeraPath, share = os.path.split ( chimeraPath )
    chimeraPath = os.path.join ( chimeraPath, 'bin' )
    chimeraPath = os.path.join ( chimeraPath, 'chimera' )
    if os.path.isfile ( chimeraPath ) :
        print " -- on unix/mac"
    else :
        chimeraPath += ".exe"
        if os.path.isfile ( chimeraPath ) :
            print " -- on windows"
        else :
            print " - chimera path not found..."

    return chimeraPath


def MinMaxCoords ( atoms ) :
    min, max = [None,None,None], [None,None,None]
    for at in atoms :
        C = at.coord()
        if min[0] == None or min[0] > C[0] : min[0] = C[0]
        if min[1] == None or min[1] > C[1] : min[1] = C[1]
        if min[2] == None or min[2] > C[2] : min[2] = C[2]
        if max[0] == None or max[0] < C[0] : max[0] = C[0]
        if max[1] == None or max[1] < C[1] : max[1] = C[1]
        if max[2] == None or max[2] < C[2] : max[2] = C[2]
    return min, max

def MinMaxPoints ( points ) :
    min, max = [None,None,None], [None,None,None]
    for C in points :
        if min[0] == None or min[0] > C[0] : min[0] = C[0]
        if min[1] == None or min[1] > C[1] : min[1] = C[1]
        if min[2] == None or min[2] > C[2] : min[2] = C[2]
        if max[0] == None or max[0] < C[0] : max[0] = C[0]
        if max[1] == None or max[1] < C[1] : max[1] = C[1]
        if max[2] == None or max[2] < C[2] : max[2] = C[2]
    return min, max



# new approach to divide up atoms by slicing the map
# - order atoms by x coordinate
# this way can better limit maximum memory use

def CalcQpn ( mol, cid, dmap, sigma, useOld=False, log=False, numProc=None, chimeraPath=None, closeMap=False, qsfile=None, res=3.0 ) :

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    mapPath = os.path.split ( dmap.data.path )[0]
    mapBase = os.path.splitext ( dmap.data.path )[0]

    dmap_name = dmap.name

    gotQ = False
    Qavg = None

    if useOld :
        SetBBAts ( mol )
        nname = mol.openedAs[0] + "__Q__" + dmap.name + ".pdb"
        if QsFromPdbFile ( mol, nname ) :
            #Qavg = QStats1 ( mol, cid )
            #return Qavg
            print " - got Q from %s" % nname
            gotQ = True
        nname = mol.openedAs[0] + "__Q__" + dmap.name + ".cif"
        if QsFromCifFile ( mol, nname ) :
            #Qavg = QStats1 ( mol, cid )
            #return Qavg
            print " - got Q from %s" % nname
            gotQ = True

    if not gotQ :

        #numProc = 2

        if numProc == None :
            import multiprocessing
            numProc = multiprocessing.cpu_count() / 2

        print "Q Scores - p - %d" % numProc
        print " - map: %s" % dmap.name
        print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
        print " - sigma: %.2f" % sigma
        minD, maxD = MinMaxD ( dmap )
        print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

        import time
        start = time.time()

        tempPath = dmap.data.path + "__Q-scores__mp__calculation__files__"
        print " - making temp path: %s" % tempPath
        try :
            os.mkdir(tempPath)
        except :
            print " - could not make temp path (an old calc may have failed):"
            print "    : check/remove temp path manually and try again"
            print "    : or, check write permissions"


        allAtsFilePath = os.path.join ( tempPath, "all_atoms.txt" )

        from gridm import Grid
        hgrid = Grid()
        hgrid.FromAtomsLocal ( mol.atoms, 1.3 )
        print " - %d ats hgrid" % len(mol.atoms)

        # write all (non-H) atoms to one file
        allAtoms = [at for at in mol.atoms if not at.element.name == "H"]
        fout = open ( allAtsFilePath, "w" )
        print " - all atoms -> %s" % allAtsFilePath
        fout.write ( "%.3f %f %f %d\n" % (sigma, minD, maxD, len(allAtoms)) )
        for at in allAtoms :
            ic = "?" if at.residue.id.insertionCode == " " else at.residue.id.insertionCode
            atId = "%d.%s.%s.%s.%s" % (at.residue.id.position,at.residue.id.chainId,ic,at.name,at.altLoc)
            p = at.coord()
            fout.write ( "%s %f %f %f\n" % (atId,p.x,p.y,p.z) )
        fout.close()

        # atoms for which to calculate Q-scores
        # SetBBAts ( mol )
        ress = []
        atoms = []
        for r in mol.residues :
            if cid == None or cid == "All" or r.id.chainId == cid :
                for at in r.atoms :
                    if not at.element.name == "H" :
                        atoms.append ( at )

        print " - atoms to do: %d" % len(atoms)

        #apoints = _multiscale.get_atom_coordinates ( atoms, transformed = False )
        #min, max = numpy.min ( apoints, axis=0 ), numpy.max ( apoints, axis=0 )

        min, max = MinMaxCoords ( atoms )
        print " -- min: %.2f, %.2f, %.2f" % (min[0], min[1], min[2])
        print " -- max: %.2f, %.2f, %.2f" % (max[0], max[1], max[2])

        dx, dy, dz = max[2] - min[2], max[1] - min[1], max[0] - min[0]

        #xatoms = []
        #for at in atoms :
        #    xatoms.append ( at.coord()[0], at )

        print " - sorting atoms..."
        atoms.sort ( reverse=False, key=lambda at: at.coord()[0] )

        print " - splitting atoms..."
        n = len(atoms)
        g = [atoms[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]

        procs = []
        #dX = (max[0] - min[0]) / float(numProc)
        totAt = 0
        #atX = min[0]
        for pi in range ( numProc ) :
            P = lambda: None
            P.atoms, P.i = [], pi
            procs.append ( P )
            #for at in atoms :
            #    C = at.coord()
            #    if C[0] >= atX and C[0] <= atX + dX :
            #        P.atoms.append ( at )
            #print "Proc %d: %.2f - %.2f -- %d atoms" % (pi, atX, atX+dX, len(P.atoms))
            P.atoms = g[pi]
            minX, maxX = P.atoms[0].coord()[0], P.atoms[-1].coord()[0]
            print "Proc %d: %.2f - %.2f -- %d atoms" % (pi, minX, maxX, len(P.atoms))
            totAt += len(P.atoms)
            #atX += dX

            P.atomsPath = os.path.join ( tempPath, "%d_atoms.txt" % P.i )
            fout = open ( P.atomsPath, "w" )
            for at in P.atoms :
                ic = "?" if at.residue.id.insertionCode == " " else at.residue.id.insertionCode
                fout.write ( "%d.%s.%s.%s.%s\n" % (at.residue.id.position,at.residue.id.chainId,ic,at.name,at.altLoc) )
            fout.close()

            P.dmapPath = os.path.join ( tempPath, "%d_map.mrc" % P.i )

            if 1 :
                nmap = MaskMapResize ( P.atoms, 7.0, dmap, P.dmapPath )
            else :
                import shutil
                shutil.copyfile ( dmap.data.path, P.dmapPath )

        print " - total proc atoms: %d" % totAt

        if closeMap :
            print " - closing %s" % dmap.name
            chimera.openModels.close ( [dmap] )


        import subprocess

        if chimeraPath == None :
            chimeraPath = GetChiPath ()

        print " -- path to Chimera:", chimeraPath

        dir_path = os.path.dirname(os.path.realpath(__file__))
        inDir = os.path.split(dir_path)[0]
        print " -- working dir:", inDir
        #mapQPPath = os.path.join ( inDir, 'Segger' )
        mapQPPath = os.path.join ( dir_path, 'mapqp.py' )
        print " -- path to mapQ script:", mapQPPath

        # start'em up
        print "Starting %d processes" % len(procs)
        for P in procs :
            args = [chimeraPath, '--nogui', '--silent', '--nostatus', P.dmapPath, mapQPPath]
            if P.i == 0 :
                print " - running proc:",
                for arg in args :
                    print arg,
                print ""

            P.fout = open ( os.path.join(tempPath, "%d.log" % P.i), "w" )
            P.foute = open ( os.path.join(tempPath, "%d_err.log" % P.i), "w" )
            P.p = subprocess.Popen(args, stdout=P.fout, stderr=P.foute, cwd=inDir)

        print ""
        print "Waiting...",
        for P in procs :
            P.p.wait()
            P.fout.close()
            P.foute.close()
            print "%d" % P.i,
        print ""

        atids = {}
        for r in mol.residues :
            for at in r.atoms :
                ic = "?" if at.residue.id.insertionCode == " " else at.residue.id.insertionCode
                atids["%d.%s.%s.%s.%s" % (at.residue.id.position,at.residue.id.chainId,ic,at.name,at.altLoc)] = at

        print ""
        print "Getting...",
        for P in procs :
            fin = os.path.join(tempPath, "%d_out.txt" % P.i)
            #print " - getting from: ", fin
            fp = open ( fin )
            for l in fp :
                try :
                    atId, Q = l.split()
                except :
                    print " - err line: ", l
                at = atids[atId.strip()]
                at.Q = float(Q)
                at.bfactor0 = at.bfactor
                at.bfactor = at.Q

                hats = hgrid.AtsNearPtLocal ( at.coord() )
                for hat, v in hats :
                    if hat.element.name == "H" :
                        hat.Q = at.Q
                        hat.bfactor0 = hat.bfactor
                        hat.bfactor = hat.Q

                for hat in at.neighbors :
                    if hat.element.name == "H" :
                        hat.Q = at.Q
                        hat.bfactor0 = hat.bfactor
                        hat.bfactor = hat.Q

            fp.close()

            if P.i == 0 :
                print ""
                print ""
                print "__StdOut for process %d__" % P.i
                foute = open ( os.path.join(tempPath, "%d.log" % P.i), "r" )
                for l in foute :
                    print l,
                print ""
                foute.close()

                print "__StdErr file for process %d__" % P.i
                foute = open ( os.path.join(tempPath, "%d_err.log" % P.i), "r" )
                for l in foute :
                    print l,
                print ""
                foute.close()

        if 1 :
            for P in procs :
                #print "Removing temp files",
                os.remove ( os.path.join(tempPath, "%d_out.txt" % P.i) )
                try :
                    os.remove ( os.path.join(tempPath, "%d_stat.txt" % P.i) )
                except :
                    print " - did not find _stat file"
                    pass
                os.remove ( os.path.join(tempPath, "%d_atoms.txt" % P.i) )
                os.remove ( os.path.join(tempPath, "%d_map.mrc" % P.i) )
                os.remove ( os.path.join(tempPath, "%d.log" % P.i) )
                os.remove ( os.path.join(tempPath, "%d_err.log" % P.i) )
                #print "%d" % P.i,

            print ""
            os.remove ( os.path.join(tempPath, "all_atoms.txt") )
            os.rmdir ( tempPath )


        end = time.time()
        print ""
        print " - done, time: %f" % ( end-start )
        totSec = end - start
        totMin = numpy.floor ( totSec / 60.0 )
        totSec = totSec - totMin * 60.0
        print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )

        SaveQFile ( mol, cid, dmap_name, sigma )
        Qavg = QStats1 ( mol, cid )



    if qsfile :
        print ""
        print "Writing Q-scores to %s" % qsfile
        from chimera.resCode import nucleic3to1
        from chimera.resCode import protein3to1

        typeSum, typeN = {}, {}
        avgq, N = 0.0, 0.0
        numNoQ = 0
        for r in mol.residues :
            tp = r.type
            if r.type in protein3to1 : tp = "Protein"
            if r.type in nucleic3to1 : tp = "Nucleic"
            if not tp in typeSum :
                typeSum[tp], typeN[tp] = 0.0, 0
            for at in r.atoms :
                if at.element.name != "H" :
                    if not hasattr ( at, 'Q' ) :
                        print " - at %s.%d[%s].%s %s - no Q" % (at.name, at.residue.id.position, at.residue.type, at.residue.id.chainId, at.altLoc)
                        numNoQ += 1
                    else :
                        typeSum[tp] += at.Q; typeN[tp] += 1.0
                        if r.type in protein3to1 or r.type in nucleic3to1 :
                            avgq += at.Q; N += 1

        print "%d / %d no q" % (numNoQ, len(mol.atoms))

        try :
            fp = open ( qsfile, "a" )
        except :
            print " - could not open", qsfile
            return

        fp.write ( "%s\t%s\t%.2f" % (dmap_name, mol.name, res) )

        if N > 0.0 :
            fp.write ( "\t%.2f" % (avgq/N) )
        else :
            fp.write ( "\t?" )

        for tp, S in typeSum.iteritems () :
            fp.write ( "\t%s\t%.2f" % (tp, S/typeN[tp]) )

        fp.write ("\n")
        fp.close()

    return Qavg





def CalcQp ( mol, cid, dmap, sigma, useOld=False, log=False, numProc=None, chimeraPath=None ) :

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap.name)[0]
    mapPath = os.path.split ( dmap.data.path )[0]
    mapBase = os.path.splitext ( dmap.data.path )[0]

    if useOld :
        SetBBAts ( mol )
        nname = mol.openedAs[0] + "__Q__" + dmap.name + ".pdb"
        if QsFromPdbFile ( mol, nname ) :
            Qavg = QStats1 ( mol, cid )
            return Qavg
        nname = mol.openedAs[0] + "__Q__" + dmap.name + ".cif"
        if QsFromCifFile ( mol, nname ) :
            Qavg = QStats1 ( mol, cid )
            return Qavg

    #numProc = 2

    if numProc == None :
        import multiprocessing
        numProc = multiprocessing.cpu_count() / 2

    print "Q Scores - p - %d" % numProc
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma

    minD, maxD = MinMaxD ( dmap )
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    import time
    start = time.time()

    tempPath = mapBase + "__Q-scores__temp__"
    print " - making temp path: %s" % tempPath
    try :
        os.mkdir(tempPath)
    except :
        print " - could not make temp path (an old calc may have failed):"
        print "    : check/remove temp path manually and try again"
        print "    : or, check write permissions"


    allAtsFilePath = os.path.join ( tempPath, "all_atoms.txt" )

    # write all (non-H) atoms to one file
    allAtoms = [at for at in mol.atoms if not at.element.name == "H"]
    fout = open ( allAtsFilePath, "w" )
    print " - all atoms -> %s" % allAtsFilePath
    fout.write ( "%.3f %f %f %d\n" % (sigma, minD, maxD, len(allAtoms)) )
    for at in allAtoms :
        r = at.residue
        altLoc = '_' if at.altLoc == '' else at.altLoc
        atId = "%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,altLoc)
        p = at.coord()
        fout.write ( "%s %f %f %f\n" % (atId,p.x,p.y,p.z) )
    fout.close()

    # atoms for which to calculate Q-scores
    SetBBAts ( mol )
    ress = []
    atoms = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            for at in r.atoms :
                if not at.element.name == "H" :
                    atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    import subprocess

    if chimeraPath == None :
        chimeraPath = GetChiPath ()

    print " -- path to Chimera:", chimeraPath

    dir_path = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.split(dir_path)[0]
    print " -- working dir:", inDir
    #mapQPPath = os.path.join ( inDir, 'Segger' )
    mapQPPath = os.path.join ( dir_path, 'mapqp.py' )
    print " -- path to mapQ script:", mapQPPath

    n = len(atoms)
    g = [atoms[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]
    procs = []
    for mi, atoms1 in enumerate(g) :

        ress1 = atoms1[0].residue
        ressN = atoms1[-1].residue
        print " - %d/%d, %d-%d - %d atoms" % (mi+1, numProc, ress1.id.position, ressN.id.position, len(atoms1))

        procAtomsPath = os.path.join ( tempPath, "%d_atoms.txt" % mi )
        fout = open ( procAtomsPath, "w" )
        for at in atoms1 :
            r = at.residue
            altLoc = '_' if at.altLoc == '' else at.altLoc
            p = at.coord()
            #fout.write ( "%d.%s.%s.%s %.3f %.3f %.3f\n" % (r.id.position,r.id.chainId,at.name,altLoc,p.x,p.y,p.z) )
            fout.write ( "%d.%s.%s.%s\n" % (r.id.position,r.id.chainId,at.name,altLoc) )
        fout.close()

        nmap_path = os.path.join ( tempPath, "%d_map.mrc" % mi )

        if 1 :
            nmap = MaskMapResize ( atoms1, 6, dmap, nmap_path )
        else :
            import shutil
            shutil.copyfile ( dmap.data.path, nmap_path )

        #args = [chimeraPath, '--nogui', '--silent', '--nostatus', mol.openedAs[0], nmap_path, mapQPPath]
        #args = [chimeraPath, '--nogui', '--silent', '--nostatus', nmap_path, dmap.data.path, mapQPPath]
        args = [chimeraPath, '--nogui', '--silent', '--nostatus', nmap_path, mapQPPath]
        if mi == 0 :
            print "running proc:",
            for arg in args :
                print arg,
            print ""

        fout = open ( os.path.join(tempPath, "%d.log" % mi), "w" )
        foute = open ( os.path.join(tempPath, "%d_err.log" % mi), "w" )
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
        fin = os.path.join(tempPath, "%d_out.txt" % mi)
        #print " - getting from: ", fin
        fp = open ( fin )
        for l in fp :
            #print " - ", l
            try :
                atId, Q = l.split()
            except :
                print " - err line: ", l
            at = atids[atId.strip()]
            #at = r.atomsMap[atName][0]
            at.Q = float(Q)
            #at.CC = float(cc)
            at.bfactor = at.Q

        fp.close()

        if mi == 0 :
            print ""
            print ""
            print "__StdOut for process %d__" % mi
            foute = open ( os.path.join(tempPath, "%d.log" % mi), "r" )
            for l in foute :
                print l,
            print ""
            foute.close()


            print "__StdErr file for process %d__" % mi
            foute = open ( os.path.join(tempPath, "%d_err.log" % mi), "r" )
            for l in foute :
                print l,
            print ""
            foute.close()

    if 1 :
        for mi, p, fout, foute in procs :
            print "Removing temp files",
            os.remove ( os.path.join(tempPath, "%d_out.txt" % mi) )
            try :
                os.remove ( os.path.join(tempPath, "%d_stat.txt" % mi) )
            except :
                print " - did not find _stat file"
                pass
            os.remove ( os.path.join(tempPath, "%d_atoms.txt" % mi) )
            os.remove ( os.path.join(tempPath, "%d_map.mrc" % mi) )
            os.remove ( os.path.join(tempPath, "%d.log" % mi) )
            os.remove ( os.path.join(tempPath, "%d_err.log" % mi) )
            print "%d" % mi,

        print ""
        os.remove ( os.path.join(tempPath, "all_atoms.txt") )
        os.rmdir ( tempPath )


    end = time.time()
    print ""
    print " - done, time: %f" % ( end-start )
    totSec = end - start
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print " - done, time: %.0f min, %.1f sec" % ( totMin, totSec )

    SaveQFile ( mol, cid, dmap.name, sigma )
    Qavg = QStats1 ( mol, cid )

    return Qavg







def QStats1 ( mol, chainId='All', doCalcResQ=True ) :

    totQ, totN = 0.0, 0.0
    #QT, QN = { "Protein":0.0, "Nucleic":0.0, "Other":0.0 }, { "Protein":0.0, "Nucleic":0.0, "Other":0.0}
    QT, QN = {}, {}
    QT_, QN_ = {}, {}
    QH, QL = {}, {}

    if chainId == None :
        chainId = "All"

    print "Q for %d res, chain %s" % ( len(mol.residues), chainId )
    for r in mol.residues :

        if r.id.chainId == chainId or chainId == "All" :

            if doCalcResQ :
                CalcResQ ( r )

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
                        QT[tp] += at.Q; QN[tp] += 1.0;
                        QH[tp] = max(QH[tp], at.Q); QL[tp] = min(QL[tp], at.Q)
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
    print "Chain\tAvg.Q-score"
    tpk = QT_.keys()
    tpk.sort()
    for tp in tpk :
        if QN_[tp] > 0 :
            avgQ = QT_[tp]/QN_[tp]
            #avgR = 0
            #if "nucleic" in tp.lower() :
            #    avgR = (avgQ-1.0673)/-0.1574
            #else :
            #    avgR = (avgQ-1.1244)/-0.1794
            #print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
            print " %s\t%.3f" % (tp, avgQ )
        else :
            print " %s\tn/a" % (tp)

    Q__ = { " protein":0, " nucleic":0, " water":0, " ion":0 }

    #for tp in ["Other", "Protein", "Nucleic"] :
    print ""
    #print "Type\tAvg.Q-score\tEst.Res.(A)"
    print "Type\tAvg.Q-score"
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
            #print " %s\t%.3f\t%.2f" % (tp, avgQ, avgR )
            print " %s\t%.3f" % (tp, avgQ )
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

    #return Q__
    if totN > 0.0 :
        return totQ/totN
    else :
        return 0.0


def QStatsProt ( mol, dmap, chainId ) :

    SetBBAts ( mol )

    ress = []
    for r in mol.residues :
        if r.id.chainId == chainId and r.isProt :
            ress.append ( r )

    if len(ress) == 0 :
        print "QstatsProt - no protein residues in chain %s" % chainId
        return

    sByType = {}
    rByType = {}

    def addType (tp, r, score) :
        if not tp in sByType :
            rByType[tp] = []
            sByType[tp] = []
        rByType[tp].append ( [score, r] )
        sByType[tp].append ( [score] )


    for r in ress :

        if r.isProt and r.type == "LEU" :
            avg = (r.atomsMap["CD1"][0].Q + r.atomsMap["CD2"][0].Q)/2.0
            addType ( "LEU(CD)", r, avg )
        if r.isProt and r.type == "LEU" and r.id.position==114 :
            avg = (r.atomsMap["CD1"][0].Q + r.atomsMap["CD2"][0].Q)/2.0
            addType ( "LEU_114(CD)", r, avg )

        if r.isProt and r.type == "VAL" :
            avg = (r.atomsMap["CG1"][0].Q + r.atomsMap["CG2"][0].Q)/2.0
            addType ( "VAL(CG)", r, avg )
        if r.isProt and r.type == "VAL" and r.id.position==33 :
            avg = (r.atomsMap["CG1"][0].Q + r.atomsMap["CG2"][0].Q)/2.0
            addType ( "VAL_33(CG)", r, avg )

        if r.isProt and r.type == "ARG" :
            avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
            addType ( "ARG(NH)", r, avg )
        if r.isProt and r.type == "ARG" and r.id.position==76 :
            avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
            addType ( "ARG_76(NH)", r, avg )
        if r.isProt and r.type == "ARG" and r.id.position==9 :
            avg = (r.atomsMap["NH1"][0].Q + r.atomsMap["NH2"][0].Q)/2.0
            addType ( "ARG_9(NH)", r, avg )

        if r.isProt and r.type == "LYS" :
            avg = r.atomsMap["NZ"][0].Q
            addType ( "LYS(NZ)", r, avg )

        if r.isProt and r.type == "ASP" :
            avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
            addType ( "ASP(OD)", r, avg )
        if r.isProt and r.type == "ASP" and r.id.position==42 :
            avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
            addType ( "ASP_42(OD)", r, avg )
        if r.isProt and r.type == "ASP" and r.id.position==131 :
            avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
            addType ( "ASP_131(OD)", r, avg )
        if r.isProt and r.type == "ASP" and r.id.position==171 :
            avg = (r.atomsMap["OD1"][0].Q + r.atomsMap["OD2"][0].Q)/2.0
            addType ( "ASP_171(OD)", r, avg )

        if r.isProt and r.type == "GLU" :
            avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
            addType ( "GLU(OE)", r, avg )
        if r.isProt and r.type == "GLU" and r.id.position==17 :
            avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
            addType ( "GLU_17(OE)", r, avg )
        if r.isProt and r.type == "GLU" and r.id.position==27 :
            avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
            addType ( "GLU_27(OE)", r, avg )
        if r.isProt and r.type == "GLU" and r.id.position==67 :
            avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
            addType ( "GLU_67(OE)", r, avg )
        if r.isProt and r.type == "GLU" and r.id.position==134 :
            avg = (r.atomsMap["OE1"][0].Q + r.atomsMap["OE2"][0].Q)/2.0
            addType ( "GLU_134(OE)", r, avg )


        if r.isProt or r.isNA :
            if r.scQ :
                addType ( r.type, r, r.scQ )
        else :
            addType ( r.type, r, r.Q )

    avgs = []
    for rtype, ra in sByType.iteritems () :
        avgs.append ( [numpy.average (ra), rtype, numpy.std (ra)] )

    from chimera.resCode import protein3to1
    from chimera.resCode import nucleic3to1

    # sort by avg score
    #avgs.sort ( reverse=True, key=lambda x: x[0] )

    # sort by residue type
    avgs.sort ( reverse=False, key=lambda x: x[1] )

    mapName = os.path.splitext(dmap.name)[0]
    molName = os.path.splitext(mol.name)[0]
    mdir, mpfile = os.path.split(dmap.data.path)
    foname = mdir + "/" + mapName + "__" + molName + ".txt"

    print " - scores to: " + foname
    fp = open (foname,"w")

    for avgScore, rtype, sdev in avgs :

        rscores = rByType[rtype]
        if len(rscores) > 0 :
            rscores.sort ( reverse=True, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)

            rts = ""
            if R.isProt : rts = protein3to1[R.type]
            elif R.isNA : rts = nucleic3to1[R.type]

            print "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore)
            fp.write ( "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore) )

    fp.close()



def QStatsRNA ( mol, dmap, chainId ) :

    SetBBAts ( mol )

    ress = []
    for r in mol.residues :
        if r.id.chainId == chainId and r.isNA :
            ress.append ( r )

    if len(ress) == 0 :
        print "Qstats RNA - no RNA residues found in chain %s" % chainId
        return

    print ""
    print "RNA stats for chain %s" % chainId
    print ""


    sByType = {}
    rByType = {}

    def addType (tp, r, score) :
        if not tp in sByType :
            rByType[tp] = []
            sByType[tp] = []
        rByType[tp].append ( [score, r] )
        sByType[tp].append ( [score] )


    scAts = []
    bbAts = []
    allAts = []

    for r in ress :
        if r.isNA :

            avg = numpy.average ( [at.Q for at in r.scAtoms] )
            #addType ( nucleic3to1[r.type] + "_SC", r, avg )
            addType ( r.type + "_SC", r, avg )

            avg = numpy.average ( [at.Q for at in r.bbAtoms] )
            #addType ( nucleic3to1[r.type] + "_BB", r, avg )
            addType ( r.type + "_BB", r, avg )

            scAts.extend ( r.scAtoms )
            bbAts.extend ( r.bbAtoms )
            allAts.extend ( [at for at in r.atoms if at.element.name != "H"] )


    avgQ = numpy.average ( [at.Q for at in allAts] )
    avgQbb = numpy.average ( [at.Q for at in bbAts] )
    avgQsc = numpy.average ( [at.Q for at in scAts] )

    sQ = numpy.std ( [at.Q for at in allAts] )
    sQbb = numpy.std ( [at.Q for at in bbAts] )
    sQsc = numpy.std ( [at.Q for at in scAts] )

    avgs = []
    for rtype, ra in sByType.iteritems () :
        avgs.append ( [numpy.average (ra), rtype, numpy.std (ra)] )


    from chimera.resCode import protein3to1
    from chimera.resCode import nucleic3to1

    # sort by avg score
    #avgs.sort ( reverse=True, key=lambda x: x[0] )

    # sort by residue type
    avgs.sort ( reverse=False, key=lambda x: x[1] )


    mapName = os.path.splitext(dmap.name)[0]
    molName = os.path.splitext(mol.name)[0]
    mdir, mpfile = os.path.split(dmap.data.path)
    foname = mdir + "/" + mapName + "__" + molName + "_rscores.txt"


    print " - scores to: " + foname
    fp = open (foname,"w")

    print ""
    print "Map\tModel\tQ_All\tQ_Backbone\tQ_SideChain\tStdQ_All\tStdQ_Backbone\tStdQ_SideChain"
    print "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f" % (mapName, molName, avgQ, avgQbb, avgQsc, sQ, sQbb, sQsc)
    print ""

    fp.write ( "\n" )
    fp.write ( "Map\tModel\tQ_All\tQ_Backbone\tQ_SideChain\tStdQ_All\tStdQ_Backbone\tStdQ_SideChain\n" )
    fp.write ( "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f" % (mapName, molName, avgQ, avgQbb, avgQsc, sQ, sQbb, sQsc) )

    fp.write ( "\n\n" )
    fp.write ( "Type\Mol\t#\tAvg.Q.\tSDev\tPos\tChain\tMaxQ\n" )
    print "RType\tResidue\t#\tAvg.Q.\tSDev\tPos\tChain\tMaxQ"
    print ""

    for avgScore, rtype, sdev in avgs :

        rscores = rByType[rtype]
        if len(rscores) > 0 :
            rscores.sort ( reverse=True, key=lambda x: x[0] )
            hr = rscores[0]
            R = hr[1]
            highestScore = hr[0]
            numRes = len(rscores)

            rts = ""
            if R.isProt : rts = protein3to1[R.type]
            elif R.isNA : rts = nucleic3to1[R.type]

            print "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore)
            fp.write ( "%s\t%s\t%d\t%f\t%f\t%d\t.%s\t%f\n" % (rtype, rts, numRes, avgScore, sdev, R.id.position, R.id.chainId, highestScore) )

    fp.close()


# expected Q-scores given the resolution of a map
#  - sigma=0.6, for resolutions 1.5 and lower
#  - sigma=0.4, for resolution higher than 1.5

def eQ_protein (RES, sigma=0.6) :
    if abs(sigma-0.6) < 1e-5 :
        return -0.1775 * RES + 1.1192, "-0.1775 * RES + 1.1192"
    elif abs(sigma-0.4) < 1e-5 :
        return -0.1866 * RES + 1.1242, "-0.1866 * RES + 1.1242"
    else :
        return None, "no eqn for protein with sigma=%.2f" % sigma

def eQ_nucleic (RES, sigma=0.6) :
    if abs(sigma-0.6) < 1e-5 :
        return -0.1377 * RES + 0.9973, "-0.1377 * RES + 0.9973"
    elif abs(sigma-0.4) < 1e-5 :
        return -0.1465 * RES + 0.9436, "-0.1465 * RES + 0.9436"
    else :
        return None, "no eqn for nucleic with sigma=%.2f" % sigma

def eQ_ion (RES, sigma=0.6) :
    if abs(sigma-0.6) < 1e-5 :
        return -0.1103 * RES + 1.0795, "-0.1103 * RES + 1.0795"
    elif abs(sigma-0.4) < 1e-5 :
        return -0.1103 * RES + 1.0795, "-0.1103 * RES + 1.0795"
    else :
        return None, "no eqn for ion with sigma=%.2f" % sigma

def eQ_water ( RES, sigma=0.6) :
    if abs(sigma-0.6) < 1e-5 :
        return -0.0895 * RES + 1.0001, "-0.0895 * RES + 1.0001"
    elif abs(sigma-0.4) < 1e-5 :
        return -0.0895 * RES + 1.0001, "-0.0895 * RES + 1.0001"
    else :
        return None, "no eqn for water with sigma=%.2f" % sigma


def ExpectedQScore (RES, sigma) :

    expQ = 1.0
    eqn = "1.0"
    x = RES
    if abs(sigma-0.6) < 1e-5 :

        #expQ = -0.0000027257*pow(RES,6) + 0.0002420020*pow(RES,5) - 0.0058960622*pow(RES,4) + 0.0611398758*pow(RES,3) - 0.2790872603*pow(RES,2) + 0.3548162472*RES + 0.8108389026
        #eqn = "-0.0000027257*POWER(RES,6) + 0.0002420020*POWER(RES,5) - 0.0058960622*POWER(RES,4) + 0.0611398758*POWER(RES,3) - 0.2790872603*POWER(RES,2) + 0.3548162472*RES + 0.8108389026"

        expQ = -0.0019064058*pow(RES,3) + 0.0499875375*pow(RES,2) - 0.4513945578*RES + 1.5361860733
        eqn = "-0.0019064058*POWER(RES,3) + 0.0499875375*POWER(RES,2) - 0.4513945578*RES + 1.5361860733"

    elif abs(sigma-0.4) < 1e-5 :

        #expQ = 0.0000121733*pow(x,6) - 0.0003668610*pow(x,5) + 0.0040966965*pow(x,4) - 0.0225427321*pow(x,3) + 0.0930645229*pow(x,2) - 0.4500934403*x + 1.3652087689
        #eqn = "0.0000121733*POWER(RES,6) - 0.0003668610*POWER(RES,5) + 0.0040966965*POWER(RES,4) - 0.0225427321*POWER(RES,3) + 0.0930645229*POWER(RES,2) - 0.4500934403*RES + 1.3652087689"

        expQ = -0.0020025943*pow(x,3) + 0.0504435042*pow(x,2) - 0.4343663048*x + 1.3945104702
        eqn = "-0.0020025943*POWER(RES,3) + 0.0504435042*POWER(RES,2) - 0.4343663048*RES + 1.3945104702"

        #eqn = -0.0000046557*POWER(RES,4) - 0.0018927526*POWER(RES,3) + 0.0495546734*POWER(RES,2) - 0.4315078614*RES + 1.3914931115


        # o	y = 0.0000121733x6 - 0.0003668610x5 + 0.0040966965x4 - 0.0225427321x3 + 0.0930645229x2 - 0.4500934403x + 1.3652087689

    return expQ, eqn


def SaveQStats ( mol, chainId, dmap_name, sigma, RES=3.0 ) :

    SetBBAts ( mol )

    if chainId == None :
        chainId = "All"

    cres = {}
    for r in mol.residues :
        if r.id.chainId == chainId or chainId == "All" :
            if r.id.chainId in cres :
                cres[r.id.chainId].append ( [r.id.position, r] )
            else :
                cres[r.id.chainId] = [ [r.id.position, r] ]

    molPath = os.path.splitext(mol.openedAs[0])[0]
    mapName = os.path.splitext(dmap_name)[0]

    nname = mol.openedAs[0] + "__Q__" + dmap_name + "_" + chainId + ".txt"
    #nname = molPath + "__Q__" + mapName + "_" + cid + ".txt"

    if hasattr ( mol, 'cif' ) :
        # cif file name updated to have the __Q__ already
        nname = mol.openedAs[0] + "_" + chainId + ".txt"

    print ""
    print "Saving per-chain & per-residue Q-scores:"
    print " -> res=", RES
    print " -> file:", nname
    print " -> chain:", chainId

    fp = open (nname, "w")

    fp.write ( "\n" )
    fp.write ( "Map: %s\n" % dmap_name )
    fp.write ( "Resolution entered (RES): %g\n" % RES )
    fp.write ( "Model: %s\n" % mol.name )
    fp.write ( "Sigma: %g\n" % sigma )
    fp.write ( "\n" )

    #avgQrna, eq_nucleic = eQ_nucleic(RES, sigma)
    #avgQprot, eq_protein = eQ_protein(RES, sigma)
    #avgQIon, eq_ion =  eQ_ion(RES, sigma)
    #avgQWater, eq_water =  eQ_water(RES, sigma)

    #fp.write ( "Protein: expectedQ = %s\n" % eq_protein )
    #fp.write ( "Nucleic: expectedQ = %s\n" % eq_nucleic )
    #fp.write ( "Ion: expectedQ = %s\n" % eq_ion )
    #fp.write ( "Water: expectedQ = %s\n" % eq_water )
    #fp.write ( "\n" )

    expQ, eqn = ExpectedQScore (RES, sigma)
    fp.write ( "expectedQ = %s\n" % eqn )
    print " -- Using expected-Q equation: ", eqn

    #fp.write ( "Chain\tType\t# residues\tAvg. Q\tExpectedQ@%.2f\tEst.Res.\n" % RES )
    fp.write ( "Chain\tType\t# residues\tAvg. Q\tExpectedQ@%.2f\n" % RES )

    chains = cres.keys()
    chains.sort()

    for cid in chains :
        ress = cres[cid]

        type_ats = {}
        type_ress = {}
        resAtoms = []
        for ri, r in ress :
            tp = ""
            if r.isProt : tp = "Protein"
            elif r.isNA : tp = "Nucleic"
            elif r.type.upper() in chargedIons : tp = "Ion"
            elif r.type.upper() == "HOH" : tp = "Water"
            else : tp = r.type

            if tp in type_ats : type_ats[tp].extend (r.atoms)
            else : type_ats[tp] = r.atoms[:]

            if tp in type_ress : type_ress[tp].append ( r )
            else : type_ress[tp] = [r]

        for rtype, atoms in type_ats.iteritems() :

            qs = [at.Q for at in atoms if (at.element.name != "H" and hasattr(at,'Q'))]
            if len(qs) == 0 :
                continue
            avgQ = numpy.average ( qs )
            numR = len ( type_ress[rtype] )

            formula, estRes = None, None
            #if "Protein" in rtype :
            #    formula = "=" + eQ_protein(RES,sigma)[1].replace ("RES",'%.2f') % RES
            #    estRes = (avgQ - 1.1192) / -0.1775
            #elif "Nucleic" in rtype :
            #    formula ="=" + eQ_nucleic(RES,sigma)[1].replace ("RES",'%.2f') % RES
            #    estRes = (avgQ - 0.9973) / -0.1377
            #elif "Ion" in rtype :
            #    formula = "=" + eQ_ion(RES,sigma)[1].replace ("RES",'%.2f') % RES
            #    estRes = (avgQ - 1.0795) / -0.1103
            #elif "Water" in rtype :
            #    formula ="=" + eQ_water(RES,sigma)[1].replace ("RES",'%.2f') % RES
            #    estRes = (avgQ - 1.0001) / -0.0895
            #else :
            #    formula = "?"
            #    estRes = 0.0
            formula = "=" + eqn.replace ("RES","%.2f"%RES)
            #fp.write ( "%s\t%s\t%d\t%.2f\t%s\t%.2f\n" % (cid, rtype, numR, avgQ, formula, estRes) )
            fp.write ( "%s\t%s\t%d\t%.2f\t%s\n" % (cid, rtype, numR, avgQ, formula) )

        #print " - cid: %s - %s - %.2f" % (cid, ctypes, cQ)
    fp.write ( "\n" )

    for cid in cres.keys () :
        rs = cres[cid]
        rs.sort()
        r = rs[0][1]

        if r.isProt :
            fp.write ( "Protein - Chain %s\t\t\t\t\t\t\t\tAverage over 3 residues\t\t\t\t\tAverage over 5 residues\t\t\t\t\tAverage over 7 residues\t\t\t\t\tAverage over 11 residues\n\n" % cid )
            fp.write ( "Chain\tRes\tRes #\tQ_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_backBone\tQ_sideChain\tQ_residue\tExpectedQ@%.2f\t\n" % RES )
        elif r.isNA :
            fp.write ( "Nucleic Acid - Chain %s\t\t\t\t\t\t\t\t\tAverage over 3 nucleotides\t\t\t\t\t\tAverage over 5 nucleotides\t\t\t\t\t\tAverage over 7 nucleotides\t\t\t\t\t\tAverage over 11 nucleotides\n\n" % cid )
            fp.write ( "Chain\tRes\tRes #\tQ_backBone\tQ_sugar\tQ_base\tQ_nucleotide\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_phosphate\tQ_sugar\tQ_base\tQ_nucleotide\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_phosphate\tQ_sugar\tQ_base\tQ_nucleotide\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_phosphate\tQ_sugar\tQ_base\tQ_nucleotide\tExpectedQ@%.2f\t\t" % RES )
            fp.write ( "Q_phosphate\tQ_sugar\tQ_base\tQ_nucleotide\tExpectedQ@%.2f\t\n" % RES )
        else :
            fp.write ( "Molecule - Chain %s\n\n" % cid )
            fp.write ( "Chain\tMolecule\tMol #\t\t\tQ_molecule\tExpectedQ@%.2f\n" % RES )


        ress = []
        Qs, AV, CC = [], [], []
        for ri, r in rs :

            #if not r.isProt and not r.isNA :
            #    print " - cid: %s - r %d - not prot or RNA" % (cid, r.id.position)
            #    continue

            ress.append (r)

            qs = [at.Q for at in r.atoms if (at.element.name != "H" and hasattr(at,'Q'))]
            if len(qs) == 0 :
                continue

            r.Q = numpy.average ( qs )

            r.qBB, r.qSC, r.qSugar, r.qResidue = 0, 0, 0, 0
            if len(r.bbAtoms) > 0 :
                r.qBB = numpy.average ( [at.Q for at in r.bbAtoms if at.element.name != "H"] )
            if len(r.scAtoms) > 0 :
                r.qSC = numpy.average ( [at.Q for at in r.scAtoms if at.element.name != "H"] )
            if len(r.sugarAtoms) > 0 :
                r.qSugar = numpy.average ( [at.Q for at in r.sugarAtoms if at.element.name != "H"] )
            r.qResidue = numpy.average ( [at.Q for at in r.atoms if at.element.name != "H"] )

            Qs.append ( [r.qBB, r.qSC, r.Q, r.qSugar] )

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


        # averages items in a list over N items before and after
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

            if not hasattr ( r, 'Q' ) or not hasattr (r, 'qBB') :
                 continue

            # fills in missing residues in proteins and rna
            if (r.isNA or r.isProt) and last_i != None :
                ii = last_i+1
                while ii < r.id.position :
                    # fill gaps
                    if r.isNA :
                        fp.write ( "%s\t%s\t%d\t" % (r.id.chainId, "", ii ) )
                        fp.write ( "\t\t\t\t%f\t\t" % (expQ) ) # avgQrna
                        fp.write ( "\t\t\t\t%f\t\t" % (expQ) )
                        fp.write ( "\t\t\t\t%f\t\t" % (expQ) )
                        fp.write ( "\t\t\t\t%f\t\t" % (expQ) )
                        fp.write ( "\t\t\t\t%f\n" % (expQ) )
                    else :
                        fp.write ( "%s\t%s\t%d\t\t\t\t%f\t\t" % (r.id.chainId, "", ii, expQ ) )
                        fp.write ( "\t\t\t%f\t\t" % (expQ) )
                        fp.write ( "\t\t\t%f\t\t" % (expQ) ) # avgQprot
                        fp.write ( "\t\t\t%f\t\t" % (expQ) )
                        fp.write ( "\t\t\t%f\n" % (expQ) )
                    ii += 1

            if r.isNA :
                fp.write ( "%s\t%s\t%d\t" % (r.id.chainId, r.type, r.id.position) )
                fp.write ( "%f\t%f\t%f\t%f\t%f\t\t" % (r.qBB, r.qSugar, r.qSC, r.Q, expQ ) )
                fp.write ( "%f\t%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,3,1), N(Qs,i,1,1), N(Qs,i,2,1), expQ ) )
                fp.write ( "%f\t%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,3,2), N(Qs,i,1,2), N(Qs,i,2,2), expQ ) )
                fp.write ( "%f\t%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,3,3), N(Qs,i,1,3), N(Qs,i,2,3), expQ ) )
                fp.write ( "%f\t%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,3,5), N(Qs,i,1,5), N(Qs,i,2,5), expQ ) )
            elif r.isProt :
                if len(r.scAtoms) > 0 :
                    fp.write ( "%s\t%s\t%d\t%f\t%f\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.qSC, r.Q, expQ ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,1,1), N(Qs,i,2,1), expQ ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,1,2), N(Qs,i,2,2), expQ ) )
                    fp.write ( "%f\t%f\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,1,3), N(Qs,i,2,3), expQ ) )
                    fp.write ( "%f\t%f\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,1,5), N(Qs,i,2,5), expQ ) )
                else :
                    fp.write ( "%s\t%s\t%d\t%f\t\t%f\t%f\t\t" % (r.id.chainId, r.type, r.id.position, r.qBB, r.Q, expQ ) )
                    fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,1), N(Qs,i,2,1), expQ ) )
                    fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,2), N(Qs,i,2,2), expQ ) )
                    fp.write ( "%f\t\t%f\t%f\t\t" % (N(Qs,i,0,3), N(Qs,i,2,3), expQ ) )
                    fp.write ( "%f\t\t%f\t%f\n" % (N(Qs,i,0,5), N(Qs,i,2,5), expQ ) )
            elif r.type.upper() in chargedIons :
                fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\n" % (r.id.chainId, r.type, r.id.position, r.Q, expQ) ) # avgQIon
            elif r.type.upper() == "HOH" :
                fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\n" % (r.id.chainId, r.type, r.id.position, r.Q, expQ) ) # avgQWater
            else :
                fp.write ( "%s\t%s\t%d\t\t\t%f\t%f\n" % (r.id.chainId, r.type, r.id.position, r.Q, expQ) ) # avgQprot

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


def CalcResQ ( r ) :

    if not hasattr ( r, 'isProt' ) :
        SetBBAts ( r.atoms[0].molecule )

    r.Q, r.qBB, r.qSC, r.qSugar = None, None, None, None

    qs = [at.Q for at in r.atoms if (at.element.name != "H" and hasattr(at,'Q'))]
    if len(qs) == 0 :
        return

    r.Q = numpy.average ( qs )

    if len(r.bbAtoms) > 0 :
        r.qBB = numpy.average ( [at.Q for at in r.bbAtoms if at.element.name != "H"] )

    if len(r.scAtoms) > 0 :
        r.qSC = numpy.average ( [at.Q for at in r.scAtoms if at.element.name != "H"] )

    if len(r.sugarAtoms) > 0 :
        r.qSugar = numpy.average ( [at.Q for at in r.sugarAtoms if at.element.name != "H"] )

    r.qResidue = r.Q



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
                #CalcResQ ( r )
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




def CalcQ ( mol, cid, dmap, sigma, useOld=False, log=False ) :

    print ""
    print "Q Scores"
    print " - map: %s" % dmap.name
    print " - mol: %s, chain: %s" % (mol.name, cid if cid != None else "_all_")
    print " - sigma: %.2f" % sigma

    minD, maxD = MinMaxD ( dmap )
    print " - mind: %.3f, maxd: %.3f" % (minD, maxD)

    ats = [at for at in mol.atoms if not at.element.name == "H"]

    allAtTree, ptGrid, atGrid = None, None, None
    if 0 :
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - atoms tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
        #allAtTree = None
    if 0 :
        import gridm
        reload(gridm)
        atGrid = gridm.Grid()
        atGrid.FromAtomsLocal ( ats, 3.0 )
        print " - %d ats grid in %.3f sec" % (len(ats), 0.0)
    if 1 :
        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        import gridm
        reload(gridm)
        ptGrid = gridm.Grid()
        ptGrid.FromPoints ( points, 3.0 )
        print " - %d pts grid" % len(points)


    atoms = []

    import time
    start = time.time()

    #ress = []
    for r in mol.residues :
        if cid == None or cid == "All" or r.id.chainId == cid :
            for at in r.atoms :
                if not at.element.name == "H" :
                    atoms.append ( at )

    print " - atoms to do: %d" % len(atoms)

    #for ai, at in enumerate ( atoms[0:2] ) :
    #    qs = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )

    from chimera import tasks, CancelOperation
    task = tasks.Task('Calculating Q-scores', modal = True)

    SetBBAts ( mol )

    modi = 100
    if len(atoms) > 100000 :
        modi = 1000

    if 1 :
    #try :

        for ai, at in enumerate ( atoms ) :

            if allAtTree :
                at.Q = Qscore ( [at], dmap, sigma, allAtTree=allAtTree, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD )
            if atGrid :
                at.Q = QscoreG ( [at], dmap, sigma, agrid=atGrid, show=0, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )
            if ptGrid :
                atPt = at.coord()
                atPt = [atPt.x, atPt.y, atPt.z]
                xfI = dmap.openState.xform
                at.Q = QscorePt3 ( atPt, xfI, dmap, sigma, ptGrid=ptGrid, log=0, numPts=8, toRAD=2.0, dRAD=0.1, minD=minD, maxD=maxD, fitg=0 )

            at.bfactor0 = at.bfactor
            at.bfactor = at.Q

            for hat in at.neighbors :
                if hat.element.name == "H" :
                    hat.Q = at.Q
                    hat.bfactor0 = hat.bfactor
                    hat.bfactor = hat.Q

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


            if ai+1 == 100 :
                if log :
                    print " - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime)

            elif (ai+1) % modi == 0 :
                if log :
                    print " - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime)

            task.updateStatus( " - Q scores - atom %d/%d - eta: %s" % (ai+1, len(atoms), leftTime) )

    try :
        print ""

    except :
        print " - something went wrong while calculating Q-scores"
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


    SaveQFile ( mol, cid, dmap.name, sigma )
    Qavg = QStats1 ( mol, cid )

    return Qavg




def SaveQFile ( mol, cid, dmap_name, sigma ) :

    if not hasattr ( mol, 'openedAs' ) :
        print ""
        print " >>> Could not save file with Q-scores - molecule was not opened from file?"
        print ""
        return

    if hasattr ( mol, 'cif' ) : # and molExt == '.cif' :
        from mmcif import WriteMol
        fout = mol.openedAs[0] + "__Q__" + dmap_name + ".cif"
        WriteMol ( mol, fout )

    else :
        nname_ = mol.openedAs[0] + "__Q__" + dmap_name + "_temp_.pdb"
        try :
            chimera.PDBio().writePDBfile ( [mol], nname_ )
        except :
            print " - could not save Q-scores file"
            return

        nname = mol.openedAs[0] + "__Q__" + dmap_name + ".pdb"

        fpo = open ( nname, "w" )
        fpi = open ( nname_ )

        ver = ""
        try :
            from Segger.mapq import mapqVersion
            ver = mapqVersion
            print " ----1- version: %s" % mapqVersion
        except :
            pass
        try :
            from mapq.mapq import mapqVersion
            ver = mapqVersion
            print " ----2- version: %s" % mapqVersion
        except :
            pass
        try :
            from mapq import mapqVersion
            ver = mapqVersion
            print " ----3- version: %s" % mapqVersion
        except :
            pass

        fpo.write ( "REMARK   0 \n" )
        fpo.write ( "REMARK   0 Q-scores calculated with MapQ\n" )
        fpo.write ( "REMARK   0   - sigma %.1f\n" % sigma )
        fpo.write ( "REMARK   0   - more info: github.com/gregdp/mapq\n"  )
        fpo.write ( "REMARK   0   - Q-scores for each atom are stored in B-factor column\n" )
        fpo.write ( "REMARK   0   - Model: %s\n" % mol.name )
        if cid == None :
            fpo.write ( "REMARK   0   - for all atoms\n" )
        else :
            fpo.write ( "REMARK   0   - for atoms in chain: %s\n" % cid )
            fpo.write ( "REMARK   0   - (other atoms have original B-factor values)\n" )
        fpo.write ( "REMARK   0   - Map: %s\n" % dmap_name )
        fpo.write ( "REMARK   0 \n" )
        for l in fpi :
            fpo.write (l)
        fpo.close()
        fpi.close()

        print " - saved %s with Q-scores" % nname

        from os import remove
        try :
            remove(nname_)
        except :
            print " - could not remove %s" % nname_



def QsFromPdbFile ( mol, qfpath ) :

    rids = {}
    for r in mol.residues :
        rids["%d.%s" % (r.id.position,r.id.chainId)] = r

    # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    try :
        fin = open ( qfpath, "r" )
    except :
        #print " - file not found"
        return False

    print " - Qs from file: %s" % qfpath

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
                            at.bfactor = bfac
                            #at.bfactor = 100.0 * (1.0 - at.Q)
                            #dval = self.cur_dmap.interpolated_values ( [ at.coord()  ], self.cur_mol.openState.xform ).astype(numpy.float64, copy=False)[0]
                            found = True
                    if not found :
                        #print " -xx- %s.%s - atom %s - loc %s" % (resi, cid, aname, aloc)
                        continue
                else :
                    #print " -xx- %s.%s - atom %s" % (resi,cid, aname)
                    continue


    fin.close ()

    return True



def QsFromCifFile ( mol, qfpath ) :


    if not os.path.isfile ( qfpath ) :
        return False

    print " - Qs from file: %s" % qfpath

    from mmcif import ReadMol
    qmol = ReadMol ( qfpath, log=False )

    rids = {}
    for r in qmol.residues :
        rids["%d.%s" % (r.id.position,r.id.chainId)] = r

    numNotFound, numQ, numNoQ = 0, 0, 0
    for at in mol.atoms :
        rid = "%d.%s" % (at.residue.id.position,at.residue.id.chainId)
        if rid in rids :
            qres = rids[rid]
            if at.name in qres.atomsMap :
                found = False
                for qat in qres.atomsMap[at.name] :
                    #print "[%s] [%s]" % (at.altLoc, qat.altLoc)
                    if at.altLoc == qat.altLoc :
                        found = True
                        #print qat.Q
                        if hasattr ( qat, 'Q' ) :
                            at.Q = qat.Q
                            numQ += 1
                        else :
                            numNoQ += 1
                if not found :
                    #print " -xx- %s.%s - atom %s - loc %s" % (resi, cid, aname, aloc)
                    continue
            else :
                #print " -xx- %s.%s - atom %s" % (resi,cid, aname)
                numNotFound += 1


    if numNotFound != 0 :
        print " - %d/%d atoms not found in q-score file" % (numNotFound, len(mol.atoms))

    print " - got Q-scores for %d/%d atoms - %d no Q" % (numQ, len(mol.atoms), numNoQ)

    return True



def QScoreFileName ( mol, dmap ) :

    molPath = mol.openedAs[0]
    mapName = dmap.name

    qfpath = ""
    if hasattr ( mol, 'cif' ) :
        qfpath = molPath + "__Q__" + mapName + ".cif"
    elif ".pdb" in molPath :
        qfpath = molPath + "__Q__" + mapName + ".pdb"
    elif ".ent" in molPath :
        qfpath = molPath + "__Q__" + mapName + ".ent"

    if not os.path.isfile ( qfpath ) :
        #print "-xxx-", qfpath
        molPath = os.path.splitext(mol.openedAs[0])[0]
        mapName = os.path.splitext(dmap.name)[0]
        if hasattr ( mol, 'cif' ) :
            qfpath = molPath + "__Q__" + mapName + ".cif"
        elif ".pdb" in mol.openedAs[0] :
            qfpath = molPath + "__Q__" + mapName + ".pdb"
        elif ".ent" in mol.openedAs[0] :
            qfpath = molPath + "__Q__" + mapName + ".ent"

    return qfpath






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

    return ptsMol




def AddSphereTubesPts ( pts, clr, rad, mname = "RAD points" ) :

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

    lastA = None
    for pt in pts :
        a = ptsMol.newAtom('', elements.H)
        res.addAtom(a)

        a.setCoord ( chimera.Point(*pt) )  # ( chimera.Point(*xyz) )
        a.radius = rad
        a.drawMode = Atom.Sphere
        a.color = chimera.MaterialColor ( *clr )
        a.surfaceCategory = 'markers'

        if lastA != None :
            b = ptsMol.newBond ( lastA, a )
            b.display = b.Smart
            b.drawMode = b.Stick
            b.radius = rad * 0.75

        lastA = a

    return ptsMol


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

    cats = {}

    for r in mol.residues :
        if r.isProt or r.isNA :
            if not r.id.chainId in cats :
                cats[r.id.chainId] = []
            cats[r.id.chainId].extend ( [at for at in r.atoms if not at.element.name == "H"] )

    chains = cats.keys()
    chains.sort()
    cid = chains[0]

    #ats = [at for at in mol.atoms if not at.element.name == "H"]
    ats = cats[cid]

    print " - %d chains, using %s, %d atoms" % ( len(chains), cid, len(ats) )

    qs, dr, q, qcc, emr = 0,0,0,0,0
    #bbRadZ, scRadZ, scRotaZ = 0,0,0
    rs, ccm, cc = 0, 0, 0
    ccm1, cc1, ccm2, cc2 = 0, 0, 0, 0
    qs, tps = None, None


    if 1 :

        points = _multiscale.get_atom_coordinates ( ats, transformed = False )
        print " - search tree: %d/%d ats" % ( len(ats), len(mol.atoms) )
        #allAtTree = AdaptiveTree ( points.tolist(), ats, 1.0)
        allAtTree = None

        sigma = 0.4
        cid = None
        #cid = mol.residues[0].id.chainId
        qs = CalcQp ( mol, cid, dmap, sigma=sigma, useOld=True )

        #print ""
        #print "Avg. Q scores:"
        #print ""
        #tps = qs.keys()
        #tps.sort()
        #for tp in tps :
        #    print " - %s : %.2f" % (tp, qs[tp])
        #print ""

    if 1 :

        rs, cc, ccm = RScore ( ats, dmap, minRes=1.0, maxRes=15.0, dRes = 0.1, xf=None )
        cc1, ccm1 = CCScore ( ats, dmap, 1.0 )
        cc15, ccm15 = CCScore ( ats, dmap, 1.5 )
        cc2, ccm2 = CCScore ( ats, dmap, 2.0 )



    if 1 :

        at = 33
        fp = None

        #fname = "scores%d_Q_allc_%s_sig%.0f.txt" % (at, label, sigma*100.0)
        fname = "scores%d_R_%s.txt" % (at, label)

        if os.path.isdir("/Users/greg/Dropbox/_mapsq") :
            fp = open ( "/Users/greg/Dropbox/_mapsq/%s" % fname, "a" )
        elif os.path.isdir("/home/greg/Dropbox/_mapsq") :
            fp = open ( "/home/greg/Dropbox/_mapsq/%s" % fname, "a" )
        elif os.path.isdir("C:/Users/greg/Dropbox/_mapsq") :
            fp = open ( "C:/Users/greg/Dropbox/_mapsq/%s" % fname, "a" )
        else :
            fp = open ( fname, "a" )

        fp.write ( "%s\t%s\t%s" % (dmap.name, mol.name, res)  )

        if 1 :
            #for tp in tps :
            #    fp.write ( "\t%s\t%.2f" % (tp, qs[tp])  )
            fp.write ( "\t%.3f" % (qs)  )

        if 1 :
            fp.write ( "\t%.3f\t%.3f\t%.3f" % (rs, cc, ccm) )
            fp.write ( "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (cc1, ccm1, cc15, ccm15, cc2, ccm2) )

        fp.write ( "\n" )

        #nProt = len ( [at for at in mol.atoms if at.residue.isProt == True] )
        #nNA = len ( [at for at in mol.atoms if at.residue.isNA == True] )
        #fp.write ( "%s\t%s\t%s\t%d\t%d\n" % (dmap.name, mol.name, res, nProt, nNA)  )

        fp.close ()




def CalcN_ () :

    print ""
    print "________Calc N_______"

    from VolumeViewer import Volume
    vols = chimera.openModels.list(modelTypes = [Volume])
    if len(vols) == 0 :
        print " - no volumes loaded"
        return
    dmap = vols[0]
    print " - dmap: %s" % dmap.name, dmap.data.path

    path1, mname = os.path.split ( dmap.data.path )
    mfname, ext = os.path.splitext ( mname )

    ts = mfname.split("_")
    res = ts[2]
    eid = ts[0] + "_" + ts[1]
    resf = float(res) / 100.0

    path2, mpath = os.path.split ( path1 )
    path3, year = os.path.split ( path2 )

    print "-"
    print " - eid: %s, year: %s, res: %.2f" % (eid, year, resf)
    print "-"


    print " - looking for cif file:"
    mfile = ""
    for f in os.listdir ( path1 ) :
        #print f
        if "__Q__" in f :
            print " - found Q file - ", f
            continue
        fname, fext = os.path.splitext ( f )
        if fext == ".cif" :
            print " - found cif file - ", f
            mfile = path1 + "/" + f

    print ""
    print "Geting mol:", mfile

    import mmcif
    mol = mmcif.LoadMol2 ( mfile )
    SetBBAts ( mol )

    cats = {}

    for r in mol.residues :
        if r.isProt or r.isNA :
            if not r.id.chainId in cats :
                cats[r.id.chainId] = []
            cats[r.id.chainId].extend ( [at for at in r.atoms if not at.element.name == "H"] )

    chains = cats.keys()
    chains.sort()
    cid = chains[0]

    #ats = [at for at in mol.atoms if not at.element.name == "H"]
    #ats = cats[cid]
    #print " - %d chains, using %s, %d atoms" % ( len(chains), cid, len(ats) )

    Q = 0
    R, cc, ccm = 0, 0, 0
    #cc1, cc1m, cc2, cc2m = 0, 0, 0, 0
    cc18, cc18m = 0, 0
    cc5, cc5m = 0, 0

    aQ2, aQ3, aQ4 = 0, 0, 0

    if 0 :

        sigma = 0.4
        cid = None
        #cid = mol.residues[0].id.chainId
        Q = CalcQp ( mol, cid, dmap, sigma=sigma, useOld=True )


    if 0 :

        ats = [at for at in mol.atoms if not at.element.name == "H"]

        print " -- R-score -- "
        R, cc, ccm = RScore ( ats, dmap, minRes=1.0, maxRes=20.0, dRes = 0.1, xf=None )
        #cc1, cc1m = CCScore ( ats, dmap, 1.0 )
        #cc15, cc15m = CCScore ( ats, dmap, 1.5 )
        #cc2, cc2m = CCScore ( ats, dmap, 2.0 )
        print " -- CC 1.777 -- "
        cc18, cc18m = CCScore ( ats, dmap, 1.777 )
        print " -- CC 5.0 -- "
        cc5, cc5m = CCScore ( ats, dmap, 5.0 )

    if 1 :

        aQ3, aQ3max = sseQscores ( mol, dmap, 3.0 )


    if 1 :

        lpath = "."
        if os.path.isdir("/Users/greg/Dropbox/_mapsq2") :
            lpath = "/Users/greg/Dropbox/_mapsq2"
        elif os.path.isdir("/home/greg/Dropbox/_mapsq2") :
            lpath = "/home/greg/Dropbox/_mapsq2"
        elif os.path.isdir("C:/Users/greg/Dropbox/_mapsq2") :
            lpath = "C:/Users/greg/Dropbox/_mapsq2"
        elif os.path.isdir("/home/stanford/_mol/_mapsqN") :
            lpath = "/home/stanford/_mol/_mapsqN"

        lpath += "/" + "%s.txt" % (year)
        print " -- log path:", lpath

        fp = open ( lpath, "a" )
        fp.write ( "%s\t%s\t%.2f" % (eid, mol.name, resf)  )

        if 0 :
            fp.write ( "\t%.3f" % (Q)  )

        if 0 :
            fp.write ( "\t%.3f\t%.3f\t%.3f" % (R, cc, ccm) )
            #fp.write ( "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (cc1, cc1m, cc15, cc15m, cc2, cc2m) )
            fp.write ( "\t%.4f\t%.4f" % (cc18, cc18m) )
            fp.write ( "\t%.4f\t%.4f" % (cc5, cc5m) )

        if 1 :
            fp.write ( "\t%.3f\t%.3f" % (aQ3, aQ3max) )

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


#[ 0.04964269]
#[ 0.08007674]
#[ 0.08772154]
#[ 0.06052513]
#[ 0.05444193]
#[ 0.05091212]
#[ 0.04454869]
#[ 0.03272544]
#[ 0.036254]
#[ 0.02918004]


def MaskMapResize ( atoms, bound, dmap, fout=None ) :


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

    #bound = 10
    #li,lj,lk = numpy.min ( fpoints, axis=0 ) - (bound, bound, bound)
    #hi,hj,hk = numpy.max ( fpoints, axis=0 ) + (bound, bound, bound)

    bx = int ( numpy.round ( bound / dmap.data.step[0] ) )
    by = int ( numpy.round ( bound / dmap.data.step[1] ) )
    bz = int ( numpy.round ( bound / dmap.data.step[2] ) )
    #print " - bound: %.2f, %.2f, %.2f" % (bx, by, bz)

    min, max = MinMaxPoints ( fpoints )
    min = numpy.array(min) - (bx, by, bz)
    max = numpy.array(max) + (bx, by, bz)
    li,lj,lk = min
    hi,hj,hk = max

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

    # todo - don't interpolate

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

    #print " - setting bbAts in %s" % mol.name
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

        r.bbAtoms = []
        r.scAtoms = []
        r.sugarAtoms = []

        if r.isProt :
            for a in r.atoms :
                if a.element.name == "H" :
                    a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False
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
                    a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False
                    continue

                n = a.name

                a.isBB = n=="P" or n=="O1P" or n=="O2P" or n=="OP1" or n=="OP2" or n=="O5'" or n=="C5'" or n=="O3'"
                a.isSugar = n=="C1'" or n=="C2'" or n=="C3'" or n=="C4'" or n=="O4'" or n=="O2'"
                #a.isBB = a.isBB or a.isSugar
                a.isBase = not a.isBB and not a.isSugar
                a.isSC = a.isBase

                #if nucleic3to1[r.type] == "G" : a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="C6" or n=="O6" or n=="N1" or n=="C2" or n=="N2" or n=="N3"
                #elif nucleic3to1[r.type] == "C" : a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="N4" or n=="C5" or n=="C6"
                #elif nucleic3to1[r.type] == "A" : a.isBase = n=="N9" or n=="C8" or n=="N7" or n=="C5" or n=="C4" or n=="N3" or n=="C2" or n=="N1" or n=="C6" or n=="N6"
                #elif nucleic3to1[r.type] == "U" : a.isBase = n=="N1" or n=="C2" or n=="O2" or n=="N3" or n=="C4" or n=="O4" or n=="C5" or n=="C6"
                #else : #print " -x- NA res %d.%s is ?" % (r.id.position, r.type)  break

                if a.isBB :
                    r.bbAtoms.append ( a )
                if a.isSugar :
                    r.sugarAtoms.append ( a )
                if a.isSC :
                    r.scAtoms.append ( a )

        else :
            for a in r.atoms :
                a.isBB, a.isSC, a.isSugar, a.isBase = False, False, False, False




def BBsegs ( ress ) :

    bbs = []

    firstRi, atRi = 0, 1
    for r in ress[1:] :
        if ress[atRi].id.position > ress[atRi-1].id.position + 1 :
            bbs.append ( ress[firstRi:atRi] )
            firstRi = atRi
        atRi += 1

    bbs.append ( ress[firstRi:atRi] )

    return bbs




def SSEs ( allRess ) :

    if len(allRess) < 1 :
        return []

    sses, ss = [], ""

    res, rStart = allRess[0], allRess[0]
    #print "  - at first res / pos: %d " % res.id.position
    if res.isHelix :
        ss = "H"
    elif res.isSheet or res.isStrand :
        ss = "E"
    else :
        ss = "_"

    ress = [ res ]
    lastRes = rStart
    for res in allRess [1:] :

        if res.id.position > lastRes.id.position + 1 :
            #print " - gap at", res.id.position
            sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
            ress = []
            rStart = res
            if res.isHelix :
                ss = "H"
            elif res.isSheet or res.isStrand :
                ss = "E"
            else :
                ss = "_"

        if res.isHelix :
            if ss != "H" :
                #print "%s -> H - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "H"
        elif res.isSheet or res.isStrand :
            if ss != "E" :
                #print "%s -> E - at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "E"
        else :
            if ss == "H" or ss == "E" :
                #print "%s -> _ at %d rid %d | %d->%d, %d res" % (ss, i, res.id.position, rStart.id.position, lastRes.id.position, len(ress))
                sses.append ( [rStart.id.position, lastRes.id.position, ss, ress] )
                ress = []
                rStart = res
                ss = "_"

        ress.append ( res )
        lastRes = res

    #print "Done at rid %d - %s | %d->%d, %d res" % ( res.id.position, ss, rStart.id.position, res.id.position, len(ress))
    sses.append ( [rStart.id.position, res.id.position, ss, ress] )
    return sses





def cubic_interp1d(x0, x, y):
    """
    Interpolate a 1-D function using cubic splines.
      x0 : a float or an 1d-array
      x : (N,) array_like
          A 1-D array of real/complex values.
      y : (N,) array_like
          A 1-D array of real values. The length of y along the
          interpolation axis must be equal to the length of x.

    Implement a trick to generate at first step the cholesky matrice L of
    the tridiagonal matrice A (thus L is a bidiagonal matrice that
    can be solved in two distinct loops).

    additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf
    """
    x = np.asfarray(x)
    y = np.asfarray(y)

    # remove non finite values
    # indexes = np.isfinite(x)
    # x = x[indexes]
    # y = y[indexes]

    # check if sorted
    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    # allocate buffer matrices
    Li = np.empty(size)
    Li_1 = np.empty(size-1)
    z = np.empty(size)

    # fill diagonals Li and Li-1 and solve [L][y] = [B]
    Li[0] = sqrt(2*xdiff[0])
    Li_1[0] = 0.0
    B0 = 0.0 # natural boundary
    z[0] = B0 / Li[0]

    for i in range(1, size-1, 1):
        Li_1[i] = xdiff[i-1] / Li[i-1]
        Li[i] = sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
        Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
        z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    i = size - 1
    Li_1[i-1] = xdiff[-1] / Li[i-1]
    Li[i] = sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
    Bi = 0.0 # natural boundary
    z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    # solve [L.T][x] = [y]
    i = size-1
    z[i] = z[i] / Li[i]
    for i in range(size-2, -1, -1):
        z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]

    # find index
    index = x.searchsorted(x0)
    np.clip(index, 1, size-1, index)

    xi1, xi0 = x[index], x[index-1]
    yi1, yi0 = y[index], y[index-1]
    zi1, zi0 = z[index], z[index-1]
    hi1 = xi1 - xi0

    # calculate cubic
    f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
         zi1/(6*hi1)*(x0-xi0)**3 + \
         (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
         (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
    return f0



def CCScore ( atoms, inMap, R = 2.0, xf=None ) :

    toMol = atoms[0].molecule
    molg = MyMolMapX2 ( atoms, R, min(R/3.0,2.0), xf )
    fpoints, fpoint_weights = fit_points_g ( molg, 1e-2 )
    map_values = inMap.interpolated_values ( fpoints, toMol.openState.xform )
    #print map_values
    olap, CC, CCm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    print " - res %.2f, cc %.3f, ccm %.3f" % (R, CC, CCm)
    return CC, CCm



def RScore ( atoms, inMap, minRes=1.0, maxRes=15.0, dRes = 1.0, xf=None ) :

    R = minRes
    maxR, maxCCm, maxCC = None, None, None
    x, y = [None] * 15, [None] * 15
    while R <= maxRes :
        toMol = atoms[0].molecule
        molg = MyMolMapX2 ( atoms, R, min(R/2.0,2.0), xf )
        fpoints, fpoint_weights = fit_points_g ( molg, 1e-2 )
        map_values = inMap.interpolated_values ( fpoints, toMol.openState.xform )
        #print map_values
        olap, CC, CCm = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
        print "%.2f\t%.3f\t%.3f" % (R, CC, CCm)
        if maxCCm == None or maxCCm < CCm :
            maxCCm = CCm
            #maxR = R
        if maxCC == None or maxCC < CC :
            maxCC = CC
            maxR = R
        R += dRes
        #print ".",

    print ""
    print " - R: %.2f, ccm %.3f, cc %.3f" % (maxR, maxCCm, maxCC)
    return maxR, maxCC, maxCCm


def fit_points_g (fdata, threshold = 1e-5) :

    mat = fdata.full_matrix()

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    from _contour import affine_transform_vertices
    affine_transform_vertices ( fpoints, fdata.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights


def MyMolMapX2 ( atoms, resolution, step=1.0, xf=None ) :

    from math import sqrt, pi

    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution

    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = False)


    # Transform coordinates to local coordinates of the molecule containing
    # the first atom.  This handles multiple unaligned molecules.
    # Or if on_grid is specified transform to grid coordinates.
    #m0 = atoms[0].molecule

    #xf = m0.openState.xform
    if xf :
        import Matrix
        #Matrix.transform_points(xyz, M.xform_matrix(xf.inverse()))
        Matrix.transform_points ( xyz, Matrix.xform_matrix(xf) )

    anum = [a.element.number for a in atoms]

    grid = bounding_grid(xyz, step, pad, [])
    grid.name = ""

    sdev = resolution * sigma_factor
    add_gaussians(grid, xyz, anum, sdev, cutoff_range, [])

    #return grid, molecules
    return grid



# -----------------------------------------------------------------------------
#
def bounding_grid(xyz, step, pad, transforms):

    xyz_min, xyz_max = point_bounds(xyz, transforms)
    origin = [x-pad for x in xyz_min]
    from math import ceil
    shape = [int(ceil((xyz_max[a] - xyz_min[a] + 2*pad) / step)) for a in (2,1,0)]
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from VolumeData import Array_Grid_Data
    grid = Array_Grid_Data(matrix, origin, (step,step,step))
    return grid


# -----------------------------------------------------------------------------
#
def add_gaussians(grid, xyz, weights, sdev, cutoff_range, transforms = []):

    from numpy import zeros, float32, empty
    sdevs = zeros((len(xyz),3), float32)
    for a in (0,1,2):
        sdevs[:,a] = sdev / grid.step[a]

    import Matrix as M
    if len(transforms) == 0:
        transforms = [M.identity_matrix()]
    from _gaussian import sum_of_gaussians
    ijk = empty(xyz.shape, float32)
    matrix = grid.matrix()
    for tf in transforms:
        ijk[:] = xyz
        M.transform_points(ijk, M.multiply_matrices(grid.xyz_to_ijk_transform, tf))
        sum_of_gaussians(ijk, weights, sdevs, cutoff_range, matrix)

    from math import pow, pi
    normalization = pow(2*pi,-1.5)*pow(sdev,-3)
    matrix *= normalization



# -----------------------------------------------------------------------------
#
def point_bounds(xyz, transforms = []):

    from _multiscale import bounding_box
    if transforms :
        from numpy import empty, float32
        xyz0 = empty((len(transforms),3), float32)
        xyz1 = empty((len(transforms),3), float32)
        txyz = empty(xyz.shape, float32)
        import Matrix as M
        for i, tf in enumerate(transforms) :
            txyz[:] = xyz
            M.transform_points(txyz, tf)
            xyz0[i,:], xyz1[i,:] = bounding_box(txyz)
        xyz_min, xyz_max = xyz0.min(axis = 0), xyz1.max(axis = 0)
    else:
        xyz_min, xyz_max = bounding_box(xyz)

    return xyz_min, xyz_max




def GetMod ( name ) :
    for m in chimera.openModels.list() :
        if m.name == name :
            return m
    return None
