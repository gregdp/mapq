import chimera
import numpy

def AxesMod ( COM=[0,0,0], U=None, Extents=[30,30,30], rad=1.0, f=1.0,
              alignTo = None ) :

    import _surface
    mol = _surface.SurfaceModel()
    chimera.openModels.add([mol], sameAs = alignTo)

    axes = AddAxes ( rad, Extents[0]*f, Extents[1]*f, Extents[2]*f, 1.0, mol )
    axes.name = "Axes"

    if U != None :
        R = numpy.array([
            [  U[0,0], U[0,1], U[0,2], 0.0    ],
            [  U[1,0], U[1,1], U[1,2], 0.0    ],
            [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )

        T = numpy.array([
            [  1.0, 0.0, 0.0, COM[0]   ],
            [  0.0, 1.0, 0.0, COM[1]   ],
            [  0.0, 0.0, 1.0, COM[2]   ]  ] )

        import Matrix
        M = Matrix.multiply_matrices ( T, R )

        ps = []
        for p in axes.surfacePieces :
            v, t = numpy.copy(p.geometry[0]), numpy.copy(p.geometry[1])
            ps.append ( [v,t,p.color] )
            axes.removePiece ( p )

        import _contour
        for p in ps :
            _contour.affine_transform_vertices( p[0], M )
            axes.addPiece ( p[0], p[1], p[2] )

        from random import random as rand
        clr = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )
        # for p in axes.surfacePieces : p.color = clr

    #for g in axes.surfacePieces :
    #    g.initial_v = numpy.copy ( g.geometry[0] )

    return axes


def AxesMod1 ( COM=[0,0,0], U=None, length=10.0, exfac=10.0, rad=1.0,
              alignTo = None, axes=None ) :

    import _surface
    toaxes = axes
    if toaxes == None :
        toaxes = _surface.SurfaceModel()
        chimera.openModels.add([toaxes], sameAs = alignTo)
        toaxes.name = "Axes"

    #axes = AddAxes ( rad, Extents[0]*f, Extents[1]*f, Extents[2]*f, 1.0, mol )
    pos = chimera.Vector(0,0,0)
    #mol = AddArrow2 ( pos, chimera.Vector(1,0,0), lX, (cF,0,0,1), rad, mol )
    #mol = AddArrow2 ( pos, chimera.Vector(0,1,0), lY, (0,cF,0,1), rad, mol )
    naxes = AddArrow3 ( pos, chimera.Vector(0,0,1), length, (0,0,1,1), rad, _surface.SurfaceModel() )

    if U != None :
        S = numpy.array([
            [  1.0, 0.0, 0.0, 0   ],
            [  0.0, 1.0, 0.0, 0   ],
            [  0.0, 0.0, 2.0, 0   ]  ] )

        P = numpy.array([
            [  1.0, 0.0, 0.0, 0   ],
            [  0.0, 1.0, 0.0, 0   ],
            [  0.0, 0.0, 1.0, -length   ]  ] )

        R = numpy.array([
            [  U[0,0], U[0,1], U[0,2], 0.0    ],
            [  U[1,0], U[1,1], U[1,2], 0.0    ],
            [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )

        T = numpy.array([
            [  1.0, 0.0, 0.0, COM[0]   ],
            [  0.0, 1.0, 0.0, COM[1]   ],
            [  0.0, 0.0, 1.0, COM[2]   ]  ] )

        import Matrix
        M = Matrix.multiply_matrices ( T, R )
        M = Matrix.multiply_matrices ( M, P )
        M = Matrix.multiply_matrices ( M, S )

        import _contour
        for p in naxes.surfacePieces :
            v, t = numpy.copy(p.geometry[0]), numpy.copy(p.geometry[1])
            _contour.affine_transform_vertices( v, M )
            toaxes.addPiece ( v,t,p.color )


    return toaxes



def AddArrow2 ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :

    if mol == None :
        import _surface
        mol = _surface.SurfaceModel()
        chimera.openModels.add ( [mol] ) # , sameAs = alignTo)
        mol.name = "Box"

    xf = AlignXf ( pos, v )
    mol = CylinderMesh2 (rad, rad, d-(rad*2), 40, clr, xf, mol )

    xf = AlignXf ( pos+(v*(d-(rad*3))), v )
    mol = CylinderMesh2 (rad*2, 0.01, rad*2, 40, clr, xf, mol )

    return mol



def AddArrow3 ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :

    xf = AlignXf ( pos, v )
    mol = CylinderMesh2 (rad, rad, d-(rad*1.5), 40, clr, xf, mol )

    xf = AlignXf ( pos+(v*(d-(rad*1.5))), v )
    mol = CylinderMesh2 (rad*3.0, 0.01, rad*1.5, 40, clr, xf, mol )

    return mol



def AddArrow4 ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None, hrad=3.0, hlen=3.0 ) :

    xf = AlignXf ( pos, v )
    mol = CylinderMesh2 (rad, rad, d-hlen, 10, clr, xf, mol )

    xf = AlignXf ( pos+(v*(d-(hlen))), v )
    mol = CylinderMesh2 (hrad, 0.01, hlen, 10, clr, xf, mol )

    return mol



def AddCylinderSolid ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :
    
    xf = AlignXf ( pos, v )
    mol = CylinderMesh2 (rad, rad, d, 10, clr, xf, mol )
    
    #xf = AlignXf ( pos+v, pos+(v*1.00001) )
    #mol = CylinderMesh2 (rad, 0.00001, 0, 10, clr, xf, mol )
    
    #xf = AlignXf ( pos, pos+(v*0.00001) )
    #mol = CylinderMesh2 (rad, 0.00001, rad, 10, clr, xf, mol )

    return mol



def CylinderSurf ( pos, v, d, clr=(0,1,1,1), rad=0.2, mol=None ) :
    
    xf = AlignXf ( pos, v )
    v, vi = CylMesh ( rad, rad, d, 8, clr, xf )
    sp = mol.addPiece ( v, vi, clr )
    return sp



def AddAxes ( rad, lX, lY, lZ, cF, mol ) :

    pos = chimera.Vector(0,0,0)
    mol = AddArrow2 ( pos, chimera.Vector(1,0,0), lX, (cF,0,0,1), rad, mol )
    mol = AddArrow2 ( pos, chimera.Vector(0,1,0), lY, (0,cF,0,1), rad, mol )
    mol = AddArrow2 ( pos, chimera.Vector(0,0,1), lZ, (0,0,cF,1), rad, mol )
    mol.name = "XYZ (RGB) Axes"

    return mol



def AlignXf ( pos, v ) :
    Z = v
    Z.normalize()
    from random import random as rand
    dZ = chimera.Vector( rand(), rand(), rand() )
    dZ.normalize()
    X = chimera.cross ( Z, dZ )
    X.normalize ()
    Y = chimera.cross ( Z, X )
    Y.normalize ()

    xf = chimera.Xform.xform (
        X.x, Y.x, Z.x, pos[0],
        X.y, Y.y, Z.y, pos[1],
        X.z, Y.z, Z.z, pos[2] )

    #xf3 = chimera.Xform.xform (
    #    d, 0, 0, 0,
    #    0, d, 0, 0,
    #    0, 0, d, 0 )
    #print xf3

    return xf


def Quad2Tri ( vi ) :
    t1 = (vi[0], vi[1], vi[2])
    t2 = (vi[0], vi[2], vi[3])
    return t1, t2



def AddVert ( verts, v ) :
    
    if verts == None :
        #return numpy.array( [ [v[0], v[1], v[2]], ], numpy.float32 )
        return numpy.array( [ v ], numpy.float32 )
    else :
        #pt = numpy.array( [ [v[0], v[1], v[2]], ], numpy.float32 )
        #return numpy.concatenate ( [verts, v] )
        return numpy.concatenate ( [verts, numpy.array( [ v ], numpy.float32 ) ] )
        




def TriangleMesh ( p1, p2, p3, color, xf, mol ) :

    if color == None :
        from random import random as rand
        color = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )


    v = None
    vi = []

    v = AddVert ( v, p1 )
    v = AddVert ( v, p2 )
    v = AddVert ( v, p3 )

    tris = [(0,1,2)]
    vi.extend(tris)

    sph = mol.addPiece ( v, vi, color )
    return mol



def TriangleMeshDiv ( p1, p2, p3, div, color, xf, mol ) :

    if color == None :
        from random import random as rand
        color = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )


    v = None
    vi = []

    v1 = p2 - p1
    v2 = p3 - p1

    v1l = numpy.sqrt ( numpy.sum ( v1*v1 ) )
    v2l = numpy.sqrt ( numpy.sum ( v2*v2 ) )
    N = numpy.ceil ( v1l / div )
    d1 = v1l / N
    d2 = v2l / N
    
    print " - v1 len: %.3f, d: %.3f, numdiv: %d, reald1: %.3f, reald2: %.3f" % (v1l, div, int(N), d1, d2 )
    
    if ( N <= 1 ) :
        v = AddVert ( v, p1 )
        v = AddVert ( v, p2 )
        v = AddVert ( v, p3 )
        vi.extend( [(0,1,2)] )
        sph = mol.addPiece ( v, vi, color )
        return mol

    v1n = v1 / v1l
    v2n = v2 / v2l
    
    v = AddVert ( v, p1 )
    prev_line_i = 0

    for i in range (1, int(N)+1 ) :
    
        #print "row ", i 

        side_p1 = p1 + i * d1 * v1n
        side_p2 = p1 + i * d2 * v2n
        
        v = AddVert ( v, side_p1 )
        prev_i = prev_line_i + i

        sv = side_p2 - side_p1
        svl = numpy.sqrt ( numpy.sum ( sv*sv ) )
        svn = sv / svl

        Nmid = i-1
        if Nmid == 0 :
            v = AddVert ( v, side_p2 )
            vi.extend( [(0,1,2)] )
            #print " v "
            prev_line_i += 1

        else :
            ds = svl / (Nmid+1.0)

            for j in range (1, Nmid+1) :

                v = AddVert ( v, j*ds*svn + side_p1 )

                tri = [(prev_line_i,prev_i,prev_i+1)]
                vi.extend( tri  )
                #print " -", j, tri

                tri = [(prev_line_i,prev_i+1,prev_line_i+1)]
                vi.extend( tri  )
                #print " -", j, tri

                prev_i += 1
                prev_line_i += 1

            v = AddVert ( v, side_p2 )
            tri = [(prev_line_i,prev_i,prev_i+1)]
            #print " +", tri
            vi.extend( tri )

            prev_line_i += 1
            
                

    sph = mol.addPiece ( v, vi, color )
    sph.displayStyle = sph.Mesh
    return sph




def CylMesh (r1, r2, Length, div, color, xf) :

    v = None
    vi = []

    # print "CylinderMesh:", div

    at = 0
    for psi_i in range(div) :

        psi = float(psi_i) * 360.0/float(div)

        #print "%.0f(%d)(%d)" % (psi, at, at-l),
        x1 = r1 * numpy.sin(psi * numpy.pi/180)
        y1 = r1 * numpy.cos(psi * numpy.pi/180)

        x2 = r2 * numpy.sin(psi * numpy.pi/180)
        y2 = r2 * numpy.cos(psi * numpy.pi/180)

        p = chimera.Point ( x1,y1,0 );
        if xf : p = xf.apply ( p )

        if psi_i == 0 :
            v = numpy.array( [ [p[0], p[1], p[2]], ], numpy.float32 )
        else :
            pt1 = numpy.array( [ [p[0], p[1], p[2]], ], numpy.float32 )
            v = numpy.concatenate ( [v, pt1] )

        p = chimera.Point ( x2,y2,Length );
        if xf : p = xf.apply ( p )
        pt2 = numpy.array( [ [p[0], p[1], p[2]], ], numpy.float32 )
        v = numpy.concatenate ( [v, pt2] )

        at = at + 2                

        if psi_i == 0 :
            pass
        else :
            tris = Quad2Tri ( [at-3, at-1, at-2, at-4] )
            vi.extend(tris)

        if psi_i == div-1 :
            tris = Quad2Tri ( [at-1, 1, 0, at-2] )
            vi.extend(tris)


    p = chimera.Point ( 0,0,0 );
    if xf : p = xf.apply ( p )
    pt1 = numpy.array( [ [p[0], p[1], p[2]], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt1] )

    p = chimera.Point ( 0,0,Length );
    if xf : p = xf.apply ( p )
    pt1 = numpy.array( [ [p[0], p[1], p[2]], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt1] )
    
    return v, vi



def CylinderMesh2 (r1, r2, Length, div, color, xf, mol) :
    
    v, vi = CylMesh ( r1, r2, Length, div, color, xf )
    sph = mol.addPiece ( v, vi, color )
    return mol






def PlaneMesh ( w, h, d, color, xf, mol ) :

    if mol == None :
        import _surface
        mol = _surface.SurfaceModel()
        chimera.openModels.add ( [mol] ) # , sameAs = alignTo)
        mol.name = "Box"

    atX = -w/2
    atY = -h/2
    
    if xf == None :
        xf = chimera.Xform()
        
    numx = int ( max ( numpy.ceil ( w / d ), 2 ) )
    numy = int ( max ( numpy.ceil ( h / d ), 2 ) )
    
    dx = w / float(numx-1)
    dy = h / float(numx-1)
    
    print " - plane - w %.2f, h %.2f, %d/%d" % (w, h, numx, numy)

    v = numpy.zeros ( [numx*numy,3] )
    vi = []
    for j in range ( numy ) :
        for i in range ( numx ) :
            v[j*numx+i] = xf.apply ( chimera.Point ( atX + dx*i, atY + dy*j, 0 ) )
            #vs.append ( p.data() )
            
            if i > 0 and j > 0 :
                p1 = j*numx+i
                p2 = j*numx+i-1
                p3 = (j-1)*numx+i-1
                p4 = (j-1)*numx+i
                vi.extend( Quad2Tri ( [p1,p2,p3,p4] ) )

    sph = mol.addPiece ( v, vi, color )
    return sph




def BoxMesh (w, h, l, color, xf, mol) :

    if mol == None :
        import _surface
        mol = _surface.SurfaceModel()
        chimera.openModels.add ( [mol] ) # , sameAs = alignTo)
        mol.name = "Box"


    v = numpy.zeros ( [8,3] )

    # print "BoxMesh:", div

    at = 0

    w = w / 2.0
    h = h / 2.0
    l = l / 2.0
    
    if xf == None :
        xf = chimera.Xform()

    
    v[0] = xf.apply ( chimera.Point ( -w,-h,-l ) )
    v[1] = xf.apply ( chimera.Point ( w,-h,-l ) )
    v[2] = xf.apply ( chimera.Point ( w,h,-l ) )
    v[3] = xf.apply ( chimera.Point ( -w,h,-l ) )
    
    v[4] = xf.apply ( chimera.Point ( -w,-h,l ) )
    v[5] = xf.apply ( chimera.Point ( w,-h,l ) )
    v[6] = xf.apply ( chimera.Point ( w,h,l ) )
    v[7] = xf.apply ( chimera.Point ( -w,h,l ) )
    
    vi = []
    vi.extend( Quad2Tri ( [0,3,2,1] ) )
    vi.extend( Quad2Tri ( [1,5,6,2] ) )
    vi.extend( Quad2Tri ( [2,6,7,3] ) )
    vi.extend( Quad2Tri ( [0,4,5,1] ) )
    vi.extend( Quad2Tri ( [0,4,7,3] ) )
    vi.extend( Quad2Tri ( [5,4,7,6] ) )

    sph = mol.addPiece ( v, vi, color )
    return sph




def BoxArrowMesh (w, h, l, al, color, xf, mol) :

    if mol == None :
        import _surface
        mol = _surface.SurfaceModel()
        chimera.openModels.add ( [mol] ) # , sameAs = alignTo)
        mol.name = "Box"


    v = numpy.zeros ( [14,3] )

    # print "BoxMesh:", div

    at = 0

    w = w / 2.0
    h = h / 2.0
    l = l / 2.0
    
    if xf == None :
        xf = chimera.Xform()

    
    v[0] = xf.apply ( chimera.Point ( -w,-h,-l ) )
    v[1] = xf.apply ( chimera.Point ( w,-h,-l ) )
    v[2] = xf.apply ( chimera.Point ( w,h,-l ) )
    v[3] = xf.apply ( chimera.Point ( -w,h,-l ) )
    
    v[4] = xf.apply ( chimera.Point ( -w,-h,l-al ) )
    v[5] = xf.apply ( chimera.Point ( w,-h,l-al ) )
    v[6] = xf.apply ( chimera.Point ( w,h,l-al ) )
    v[7] = xf.apply ( chimera.Point ( -w,h,l-al ) )
    
    vi = []
    vi.extend( Quad2Tri ( [0,3,2,1] ) )
    vi.extend( Quad2Tri ( [1,5,6,2] ) )
    vi.extend( Quad2Tri ( [2,6,7,3] ) )
    vi.extend( Quad2Tri ( [0,4,5,1] ) )
    vi.extend( Quad2Tri ( [0,4,7,3] ) )
    vi.extend( Quad2Tri ( [5,4,7,6] ) )

    v[8] = xf.apply ( chimera.Point ( -w,-h-al/2.0,l-al ) )
    v[9] = xf.apply ( chimera.Point ( w,-h-al/2.0,l-al ) )
    v[10] = xf.apply ( chimera.Point ( w,h+al/2.0,l-al ) )
    v[11] = xf.apply ( chimera.Point ( -w,h+al/2.0,l-al ) )

    v[12] = xf.apply ( chimera.Point ( -w,0,l ) )
    v[13] = xf.apply ( chimera.Point ( w,0,l ) )

    vi.extend( Quad2Tri ( [8,11,10,9] ) )
    vi.extend( Quad2Tri ( [8,9,13,12] ) )
    vi.extend( Quad2Tri ( [11,12,13,10] ) )
    vi.extend( [(8,12,11)] )
    vi.extend( [(9,10,13)] )
    
    sph = mol.addPiece ( v, vi, color )
    return sph





def SphereMesh (r, div, color, pos, mol) :

    posx,posy,posz = pos[0], pos[1], pos[2]
    v = numpy.array( [ [0+posx,0+posy,r+posz], ], numpy.float32 )
    vi = ()

    at = 1
    l = int ( numpy.ceil (float(div)*3.0/2.0) )
    if div < 10 : l = div*2
    #print "SphereMesh:", div, 'x', l
    lat = 0

    for phi_i in range(div) :

        phi = 90.0 - ( float(phi_i+1) * 180.0/float(div+1) )
        #print "%.2f: " % phi,
        z = r * numpy.sin(phi * numpy.pi/180)
        s = r * numpy.cos(phi * numpy.pi/180)

        for psi_i in range (l) :
            psi = float(psi_i) * 360.0/float(l)

            #print "%.0f(%d)(%d)" % (psi, at, at-l),
            x = s * numpy.sin(psi * numpy.pi/180)
            y = s * numpy.cos(psi * numpy.pi/180)

            pt = numpy.array( [ [x+posx,y+posy,z+posz], ], numpy.float32 )
            v = numpy.concatenate ( [v, pt] )

            if phi_i == 0 :
                if psi_i > 0 :
                    vi = vi + ( (at-1, at, 0), )
                if psi_i == l-1 :
                    vi = vi + ( (at, 1, 0), )
            else :
                if psi_i > 0 :
                    tris = Quad2Tri ( [at-1, at, at-l, at-l-1] )
                    vi = vi + tris
                if psi_i == l-1 :
                    tris = Quad2Tri ( [at, at-l+1, at-l*2+1, at-l] )
                    vi = vi + tris

            if phi_i == div-1 :
                if psi_i > 0 :
                    vi = vi + ( (at, at-1, lat+l), )
                if psi_i == l-1 :
                    vi = vi + ( (at-l+1, at, lat+l), )

            at = at + 1                
                

        lat = len ( v )

    pt = numpy.array( [ [0+posx,0+posy,-r+posz], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt] )


    sph = mol.addPiece( v, vi, color )
    return sph




def prAxes ( points ) :

    com = numpy.sum(points, axis=0) / len(points)
    C = chimera.Vector ( com[0], com[1], com[2] )

    comv = numpy.ones_like ( points ) * com
    points = points - comv

    i = numpy.matrix ( [[1,0,0], [0,1,0], [0,0,1]] )
    ii = i * numpy.sum ( numpy.multiply ( points, points ) )
    p_t = numpy.transpose(points)
    td = numpy.tensordot ( points, p_t, axes=[0,1] )

    I0 = ii - td

    try :
        U, S, V = numpy.linalg.svd( I0 )
    except :
        print "- error computing SVD - prob. singular matrix"
        return []

    #U[0,0] = U[0,0] * -1.0
    #U[1,0] = U[1,0] * -1.0
    #U[2,0] = U[2,0] * -1.0

    #U[0,2] = U[0,2] * -1.0
    #U[1,2] = U[1,2] * -1.0
    #U[2,2] = U[2,2] * -1.0

    return [C, U, S, V]



def map_points (fmap, useThreshold = True):

    from _contour import affine_transform_vertices as transform_vertices

    mat = fmap.data.full_matrix()
    threshold = fmap.surface_levels[0]
    
    if useThreshold == False :
        #threshold = -1e9
        threshold = 1e-5
        #print " - not using threshold"

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    transform_vertices( fpoints, fmap.data.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights



def AxesModOffset ( COM=[0,0,0], U=None, Extents=[30,30,30], rad=1.0, f=1.0,
			 alignTo = None ) :
	
    import _surface
    mol = _surface.SurfaceModel()
    chimera.openModels.add([mol], sameAs = alignTo)
	
    pos = chimera.Vector(0,0,0)
    axes = AddArrow2 ( pos, chimera.Vector(0,1,0), lY, (cF,.3,.3,1), rad, mol )
	
    axes.name = "Riboarrow"
	
    if U != None :
        R = numpy.array([
						 [  U[0,0], U[0,1], U[0,2], 0.0    ],
						 [  U[1,0], U[1,1], U[1,2], 0.0    ],
						 [  U[2,0], U[2,1], U[2,2], 0.0    ]  ] )
		
        T = numpy.array([
						 [  1.0, 0.0, 0.0, COM[0]   ],
						 [  0.0, 1.0, 0.0, COM[1]   ],
						 [  0.0, 0.0, 1.0, COM[2]   ]  ] )
		
        Ti = numpy.array([
						  [  1.0, 0.0, 0.0, Extents[0]*0.7  ],
						  [  0.0, 1.0, 0.0, -Extents[1]/2.0   ],
						  [  0.0, 0.0, 1.0, -Extents[0]*0.7   ]  ] )
		
        import Matrix
        M = Matrix.multiply_matrices ( R, Ti )
        M = Matrix.multiply_matrices ( T, M )
		
        ps = []
        for p in axes.surfacePieces :
            v, t = numpy.copy(p.geometry[0]), numpy.copy(p.geometry[1])
            ps.append ( [v,t,p.color] )
            axes.removePiece ( p )
		
        import _contour
        for p in ps :
            _contour.affine_transform_vertices( p[0], M )
            axes.addPiece ( p[0], p[1], p[2] )
		
        from random import random as rand
        clr = ( rand()*.7, rand()*.7, rand()*.7, 1.0 )
	# for p in axes.surfacePieces : p.color = clr
	
    #for g in axes.surfacePieces :
    #    g.initial_v = numpy.copy ( g.geometry[0] )
	
    return axes

