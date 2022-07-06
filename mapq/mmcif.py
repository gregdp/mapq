




# execfile ( "Segger/mmcif.py" ); r = ReadCif ( "/Users/greg/Box Sync/_data/problems/emd_30342/7cec_s.cif" )

def ReadCif ( fpath, log=False ) :

    print "Reading %s" % fpath

    from os import path

    if not path.isfile ( fpath ) :
        print " - could not find"
        return


    import time
    start = time.time()

    cif_name = ""
    cif = []
    loops = {}

    fp = open ( fpath )

    li = 0
    getNext = True
    while 1 :

        if getNext :
            atLine, getNext = fp.readline(), True
            li += 1
        else :
            getNext = True

        if not atLine :
            print " - done"
            break

        ls = atLine.strip ()

        #print "%d %s" % (li, ls)

        if len(ls) == 0 :
            continue

        elif ls[0:4] == "data" :
            cif_name = ls
            print " - cif name:", cif_name
            cif.append ( atLine )

        elif ls[0] == "#" :
            cif.append ( atLine )

        elif ls[0] == "_" :
            #cif.append ( ls )
            li, name, data, lines = GetData1 ( fp, atLine, ls, li )
            cif.append ( lines )

        elif ls[0:5] == "loop_" :

            li, atLine, ls, getNext, name, labels, data = GetLoop ( fp, atLine, ls, li )
            cif.append ( [name, labels, data] )
            loops[name] = { 'labels':labels, 'data':data }
            #print " - got loop: %s" % name

        else :
            #cif.append ( atLine )
            print " - ? %d -.- %s" % (li, ls)

    fp.close()

    print ( " - done CIF - %d lines -- %.1fs" % (li, time.time()-start) )

    if log :
        print ""
        print "Loops:"
        for name, ld in loops.iteritems() :
            labels, data = ld['labels'], ld['data']
            print " - %s [%d x %d]" % ( name, len(labels), len(data) )

    #outf = path.splitext ( fpath )[0] + "_w.cif"
    #WriteCif ( cif, outf )

    return cif, loops


def GetLoop ( fp, atLine, ls, li ) :

    loopLabels, loopData, loopName = [], [], None

    getData = False
    getNext = True

    # get loop labels first, then data
    while 1 :

        atLine, getNext = fp.readline(), True
        if not atLine :
            print " - %d - eof while getting loop" % li
            break
        li += 1
        ls = atLine.strip()
        #print " _ %d %s %d" % (li, ls, len(ls))

        if len(ls) == 0 :
            #print " - done loop %s" % li
            #print "[" + ls + "]"
            getNext = False
            break

        elif ls[0] == "_" :
            if getData :
                print " - done loop? / %d" % li, ls
                print loopLabels
                #print loopData
                print loopName
                print ""
                getNext = False
                break
            else :
                # loop label
                name, label = ls.split(".")
                loopLabels.append ( label )
                loopName = name

        elif ls[0] == "#" :
            #print " - done loop %s" % li
            #print "[" + ls + "]"
            getNext = False
            break

        else :
            getData = True

            if len(loopLabels) == 0 :
                print " - ? %d - no labels for loop" % li
                return

            li, data = GetData ( fp, atLine, ls, li, loopLabels )

            if len(data) != len(loopLabels) :
                print " - ? %d - labels/data different sizes %d/%d" % (li, len(data), len(loopLabels))
                print data
                return

            else :
                mdata = {}
                for i, label in enumerate ( loopLabels ) :
                    mdata[label] = data[i]
                #loopData.append ( [data, mdata] )
                loopData.append ( {'asArray':data, 'asMap':mdata} )

                if 0 :
                    print "[%d]" % li
                    for label, data in mdata.iteritems () :
                        if type(data) is list:
                            print " - %s:%s" % (label, data[0])
                        else :
                            print " - %s:%s" % (label, data)
                    #print pdict


    #print " - returning from getloop - ", ls
    return li, atLine, ls, getNext, loopName, loopLabels, loopData


def GetData1 ( fp, atLine, ls, li ) :
    # get data for a single label

    name, data, lines = "", [], ""

    tsi = splitm ( li, ls )

    if len(tsi) == 2 :
        lines += atLine
        name, data = tsi

    elif len(tsi) == 1 :
        name = tsi[0]
        lines += atLine
        liStart = li

        atLine, getNext = fp.readline(), True
        if not atLine :
            print " - ? %d - eof while getting single value" % li
            return li, name, data, lines

        lines += atLine
        li += 1
        ls = atLine.strip ()
        if len(ls) == 0 :
            data = ""
        elif ls[0] == ';' :
            data = atLine

            # look for end ;
            while 1 :
                atLine, getNext = fp.readline(), True
                if not atLine :
                    print " - ? %d - eof while getting single value starting at %d" % (li, liStart)
                    return li, name, data, lines

                li += 1
                lines += atLine
                ls = atLine.strip ()
                if len(ls) == 0 :
                    data += atLine
                elif ls[0] == ';' :
                    # done
                    break
                else :
                    data += atLine
        else :
            data = ls

    else :
        print " - ? getdata1 line %d - " % li, ls
        print " - %d tokens" % len(tsi)
        print tsi
        print ls.split()
        print ""

    return li, name, data, lines


def GetData ( fp, atLine, ls, li, labels ) :
    # get (multiple) data for loop

    data = []
    liStart = li

    while 1 :

        if len(ls) == 0 :
            print " - blank line %d" % li, atLine

        elif ls[0] == ';' :
            block = atLine
            t = ls[1:]
            while 1 :
                atLine, getNext = fp.readline(), True
                li += 1
                if not atLine :
                    print " - ? %d - reached eof while scanning block" % liStart
                    break
                block += atLine
                ls = atLine.strip()
                #print " ; %d %s" % (li, ls)
                if len(ls) == 0 :
                     continue
                if ls[0] == ';' :
                    break
                else :
                    t += " " + ls
            #data.append ( [t, block] )
            data.append ( {'string':t, 'lines':block} )

        else :
            tsi = splitm ( li, ls )
            data.extend ( tsi )

        if len(data) >= len(labels) :
            # done
            break

        else :
            # keep going...
            atLine, getNext = fp.readline(), True
            li += 1
            if not atLine :
                print " - ? %d - reached eof while getting data" % liStart
                break
            ls = atLine.strip()


    return li, data


def splitm ( li, l ) :
    ts = l.split()
    tsr = []
    addTo = None
    endChar = None
    for t in ts :
        #print " - %s, %d, %s" % (t, len(t), addTo)
        if addTo != None :
            addTo += " " + t
            if t[-1] == endChar :
                tsr.append ( addTo[0:-1] )
                addTo = None
        elif t[0] == "'" :
            if len(t)>1 and t[-1] == "'" :
                tsr.append ( t[1:-1] )
            else :
                addTo = t[1:]
                endChar = "'"
        elif t[0] == '"' :
            if len(t)>1 and t[-1] == '"' :
                tsr.append ( t[1:-1] )
            else :
                addTo = t[1:]
                endChar = '"'
        else :
            tsr.append ( t )
    if addTo :
        print " - ? %d - unmatched '" % li
        tsr += addTo

    #t2 = []
    for i,t in enumerate(tsr) :
        if len(t) == 0 :
            #print " - length 0 on line %d" % li
            tsr[i] = t[1:-1]
            continue
        if t[0] == '"' and t[-1] == '"' :
            tsr[i] = t[1:-1]

    return tsr



def WriteCif ( cif, fout ) :

    print ""
    print " - writing to %s" % fout

    fp = open ( fout, "w" )

    for ls in cif :
        if type(ls) == list :
            name, labels, data = ls
            fp.write ( "loop_\n" )
            for label in labels :
                fp.write ( "%s.%s\n" % (name, label) )

            # array of column widths - each entry has the max width for the column
            cws = [0] * len(labels)
            for d in data :
                adata, mdata = d['asArray'], d['asMap']

                for i, ds in enumerate(adata) :
                    dd = "'%s'" % ds if ' ' in ds else ds
                    cws[i] = max ( cws[i], len(dd) )

            for d in data :
                adata, mdata = d['asArray'], d['asMap']
                first, lastWasLines = True, False
                for di, dd in enumerate(adata) :
                    if type(dd) == dict :
                        # write original block starting with ; on start and end lines
                        if not first :
                            fp.write ( "\n" )
                        fp.write ( dd['lines'] )
                        first = True
                        lastWasLines = True
                    else :
                        # write in columns
                        dd = "'%s'" % dd if ' ' in dd else dd
                        #dd = dd if first else ("\t" + dd)
                        padn = cws[di] - len(dd) + 1
                        fp.write ( "%s%s" % (dd, " "*padn) )
                        first = False
                        lastWasLines = False
                if not lastWasLines :
                    fp.write ("\n")
        else :
            fp.write ( ls )

    fp.close()




try :
    import chimera
    import numpy

    from chimera.resCode import nucleic3to1
    from chimera.resCode import protein3to1, protein1to3
    protein3to1['HSD'] = protein3to1['HIS']
    protein3to1['HSE'] = protein3to1['HIS']

    nucleic1to3 = { 'T':'THY', 'C':'CYT', 'G':'GUA', 'A':'ADE', 'U':'URA'}
    nucleic3to1['GDP'] = nucleic3to1['GUA']

    atomColors = {'C' : chimera.MaterialColor (0.565,0.565,0.565),
                'Cbb' : chimera.MaterialColor (0.2,0.6,0.2),
                'S' : chimera.MaterialColor (1.000,1.000,0.188),
                'O' : chimera.MaterialColor (1.000,0.051,0.051),
                'N' : chimera.MaterialColor (0.188,0.314,0.973),
                'P' : chimera.MaterialColor (1.0, 0.502, 0.0),
                'H' : chimera.MaterialColor (0.9,.9,.9),
                ' ' : chimera.MaterialColor (0.2,1,.2),
                "MG" : chimera.MaterialColor (0,1,0),
                "NA" : chimera.MaterialColor (.6,.3,.6),
                "CL" : chimera.MaterialColor (.2,.6,.2),
                "CA" : chimera.MaterialColor (.4,.4,.6),
                "ZN" : chimera.MaterialColor (.2,.8,.2),
                "MN" : chimera.MaterialColor (.4,.4,.6),
                "FE" : chimera.MaterialColor (.4,.4,.6),
                "CO" : chimera.MaterialColor (.4,.4,.6),
                "NI" : chimera.MaterialColor (.4,.4,.6)
    }


except :
    pass


# this makes a molecule model but does not add bonds
def LoadMol ( fpath, log=False ) :

    mol = ReadMol ( fpath, log )

    from time import time
    startt = time()

    chimera.openModels.add ( [mol] )
    #return mol

    print " - added %s in %.2f sec" % (mol.name, time() - startt)

    for at in mol.atoms :
        at.display = True
        at.drawMode = at.Sphere
        at.color = mol.chainColors[at.residue.id.chainId]

    for res in mol.residues :
        res.ribbonDisplay = False # drawRib
        res.ribbonDrawMode = 2
        res.ribbonColor = mol.chainColors[at.residue.id.chainId]

    #if hasattr ( mol, 'chainDescr' ) :
    #    for cid, descr in mol.chainDescr.iteritems() :
    #        print " - %s - %s" % (cid, ", ".join(descr))

    return mol


def ParamPathPdb ( rtype ) :

    #ppath = "/Users/greg/Dropbox/_mol/Segger/_param/"

    from os import path
    dir_path = path.dirname ( path.realpath(__file__) )
    inDir = path.split(dir_path)[0]
    #print " -- working dir:", inDir
    #mapQPPath = os.path.join ( inDir, 'Segger' )
    ppath = path.join ( dir_path, '_param' )
    #print " -- path to param:", ppath
    #fname = ppath + "%s.pdb" % rtype
    fname = path.join ( ppath, "%s.pdb" % rtype )
    return fname


def GetResMol ( rtype ) :

    from os import path
    dir_path = path.dirname ( path.realpath(__file__) )
    inDir = path.split(dir_path)[0]
    ppath = path.join ( dir_path, '_param' )
    fname = path.join ( ppath, "%s.pdb" % rtype )

    if not path.isfile(fname) :
        print " - did not find %s" % fname

        phPath = "/Users/greg/_mol/phenix-installer-1.19.2-4158-mac-intel-osx-x86_64/build/bin/phenix.elbow"
        if not path.isfile ( phPath ) :
            print " - %s - phenix.elbow not found" % rtype
            return None

        #args = ["/Users/greg/_mol/phenix-1.19.2-4158/build/bin/phenix.elbow", "--chemical_component", rtype]
        args = [phPath, "--chemical_component", rtype]

        print "Running elbow:"
        print args

        fname_log = path.join ( ppath, "%s.log" % rtype )
        fname_err = path.join ( ppath, "%s_err.log" % rtype )

        fout = open ( fname_log, "w" )
        foute = open ( fname_err, "w" )
        import subprocess
        p = subprocess.Popen(args, stdout=fout, stderr=foute, cwd=ppath)

        print " - waiting..."
        p.wait()
        fout.close()
        foute.close()

        if not path.isfile(fname) :
            print " - elbow file not found %s" % fname
            return None

    import chimera
    nmol = chimera.PDBio().readPDBfile ( fname )[0]
    #print " - read %s - %d atoms - %d res" % ( fname, len(nmol.atoms), len(nmol.residues) )
    #addRes = nmol.residues[0]
    return nmol



def LoadMol2 ( fpath, log=False ) :

    if 0 :
        from chimera import tasks, CancelOperation
        task = tasks.Task('...', modal = True)

        try :
            mol = LoadMol2_ ( fpath, log, task )
            #mol = LoadMolCh_ ( fpath, log, task )

        except Exception, err:
            #umsg ( "Something went wrong..." )
            print Exception, err
            traceback.print_exc()
            return

        finally :
            task.finished()

        return mol

    else :
        mol = LoadMol2_ ( fpath, log, None )
        #mol = LoadMolCh_ ( fpath, log, None )
        return mol



# make molecule model and add bonds using phenix.elbow

def LoadMol2_ ( fpath, log=False, task=None ) :

    mol = ReadMol ( fpath, log=False )

    crmap = {}
    rmolmap = {}

    print " - %d residues" % len(mol.residues)

    import time
    start = time.time()

    for r in mol.residues :
        if not r.id.chainId in crmap :
            crmap[r.id.chainId] = { r.id.position : r }
        else :
            crmap[r.id.chainId][r.id.position] = r

    bonds = []
    for ri, r in enumerate ( mol.residues ) :
        if task :
            if ri % 100 == 0 :
                task.updateStatus( "%s - residue %d/%d" % ( mol.name, ri, len(mol.residues) ) )

        rmol = None

        rtype = r.type.upper()
        #if rtype in nucleic1to3 :
        #    rtype = nucleic1to3[rtype]
        #    #print r.type, "->", rtype

        #print "%d %d.%s %s" % (ri, r.id.position, r.id.chainId, r.type)

        if rtype.lower() in rmolmap :
            rmol = rmolmap[ rtype.lower() ]
        else :
            rmol = GetResMol ( rtype.lower() )
            #if rmol != None :
            rmolmap[rtype.lower()] = rmol

        if 1 and rmol != None :
            for b in rmol.bonds :
                a1n, a2n = b.atoms[0].name, b.atoms[1].name
                if a1n in r.atomsMap and a2n in r.atomsMap :
                    for a1 in r.atomsMap[a1n] :
                        for a2 in r.atomsMap[a2n] :
                            #print "%s - %s" % ( At(a1), At(a2) )
                            #nb = mol.newBond ( a1, a2 )
                            if a1.altLoc == a2.altLoc :
                                bonds.append ( [a1,a2] )

        if 1 :
            if r.type.upper() in protein3to1 :
                if r.id.position-1 in crmap[r.id.chainId] :
                    pres = crmap[r.id.chainId][r.id.position-1]
                    if pres.type.upper() in protein3to1 :
                        #GetSS ( pres, r )
                        if "C" in pres.atomsMap and "N" in r.atomsMap :
                            for a1 in pres.atomsMap["C"] :
                                for a2 in r.atomsMap["N"] :
                                    #print a1.name, pres.id.position, a2.name, r.id.position
                                    #nb = mol.newBond ( a1, a2 )
                                    if a1.altLoc == a2.altLoc :
                                        bonds.append ( [a1,a2] )

            if r.type.upper() in nucleic1to3 or r.type.upper() in nucleic3to1 :
                if r.id.position-1 in crmap[r.id.chainId] :
                    pres = crmap[r.id.chainId][r.id.position-1]
                    if pres.type.upper() in nucleic1to3 or pres.type.upper() in nucleic3to1 :
                        if "O3'" in pres.atomsMap and "P" in r.atomsMap :
                            for a1 in pres.atomsMap["O3'"] :
                                for a2 in r.atomsMap["P"] :
                                    #print a1.name, pres.id.position, a2.name, r.id.position
                                    #nb = mol.newBond ( a1, a2 )
                                    if a1.altLoc == a2.altLoc :
                                        bonds.append ( [a1,a2] )

    print ( " - %d bonds, %.1fs" % (len(bonds), time.time()-start)  )

    # adding bonds at end was faster?
    start = time.time()
    for a1, a2 in bonds :
        nb = mol.newBond ( a1, a2 )

    print ( " - added bonds in %.1fs" % (time.time()-start)  )

    start = time.time()
    chimera.openModels.add ( [mol] )
    print " - added mol, %.1fs" % (time.time()-start)

    return mol


def ColorMol ( mol ) :

    resDisp = {}

    if hasattr ( mol, "cifLoops" ) :
        if  "_mapq_comp_display_attributes" in mol.cifLoops :
            print " - found res disp"
            rcols = mol.cifLoops['_mapq_comp_display_attributes']['data']
            for d in rcols :
                mp = d['asMap']
                cid = mp["label_asym_id"]
                rpos = mp["label_seq_id"]
                rdisp = int ( mp["ribbonDisplay"] )
                r = float ( mp["ribbonColor_R"] )
                g = float ( mp["ribbonColor_G"] )
                b = float ( mp["ribbonColor_B"] )
                a = float ( mp["ribbonColor_A"] )
                resDisp[cid+rpos] = [ True if rdisp==1 else False, (r,g,b,a) ]




    if not hasattr ( mol, 'chainColors' ) :
        from random import random
        mol.chainColors = {}
        for r in mol.residues :
            if not r.id.chainId in mol.chainColors :
                clr = chimera.MaterialColor ( random(), random(), random(), 1.0 )
                mol.chainColors[r.id.chainId] = clr

    for r in mol.residues :
        rt = r.type.upper()
        rid = "%s%d" % (r.id.chainId,r.id.position)
        if rt in protein3to1 or rt in nucleic3to1 or rt in nucleic1to3 :
            r.ribbonDrawMode = 2
            if rid in resDisp :
                rd = resDisp[rid]
                r.ribbonDisplay = rd[0]
                r.ribbonColor = chimera.MaterialColor( *rd[1] )
            else :
                r.ribbonDisplay = True
                r.ribbonColor = mol.chainColors[r.id.chainId]
            for at in r.atoms :
                at.display = False
                at.drawMode = at.EndCap
                if at.element.name.upper() in atomColors :
                    at.color = atomColors[at.element.name.upper()]
                else :
                    at.color = mol.chainColors[r.id.chainId]
        else :
            for at in r.atoms :
                at.display = True
                at.drawMode = at.EndCap
                if at.element.name.upper() in atomColors :
                    at.color = atomColors[at.element.name.upper()]
                else :
                    at.color = mol.chainColors[r.id.chainId]

    for b in mol.bonds :
        b.drawMode = b.Stick
        b.display = b.Smart


def At ( at ) :
    return "%d.%s(%s)_%s" % (at.residue.id.position, at.residue.id.chainId, at.residue.type, at.name)



# this makes one molecule for each chain
def LoadMolCh_ ( fpath, log=False, task=None ) :

    mol = ReadMol ( fpath, log=False )

    from chimera.resCode import nucleic3to1
    from chimera.resCode import protein3to1, protein1to3
    protein3to1['HSD'] = protein3to1['HIS']
    protein3to1['HSE'] = protein3to1['HIS']

    nucleic1to3 = { 'T':'THY', 'C':'CYT', 'G':'GUA', 'A':'ADE', 'U':'URA'}
    nucleic3to1['GDP'] = nucleic3to1['GUA']

    crmap = {}
    rmolmap = {}

    print " - adding bonds"

    import time
    start = time.time()

    for r in mol.residues :
        if not r.id.chainId in crmap :
            crmap[r.id.chainId] = { r.id.position : r }
        else :
            crmap[r.id.chainId][r.id.position] = r

    chains = mol.chainColors.keys()

    print "%d residues - %d chains" % ( len(mol.residues), len(chains) )
    from os.path import splitext

    chMols = {}
    for ch in chains :

        start = time.time()

        chMol = chimera.Molecule()
        chMol.name = splitext ( mol.name )[0] + "_" + ch
        rmap0 = crmap[ch]
        chMols[ch] = chMol
        print " - %s - %d residues" % (ch, len(rmap0)),

        rmap = {}
        for ri, res in rmap0.iteritems() :
            nr = chMol.newResidue ( res.type, chimera.MolResId(ch, res.id.position) )
            rmap[nr.id.position] = nr
            for at in res.atoms :
                nat = chMol.newAtom ( at.name, at.element )
                nr.addAtom ( nat )
                nat.setCoord ( at.coord() )

        print ", %d atoms" % len(chMol.atoms),

        for ri, r in rmap.iteritems() :
            #if ri % 100 == 0 :
            #    print "%d/%d" % ( ri, len(mol.residues) )
            #    if task :
            #        task.updateStatus( "%d/%d" % ( ri, len(mol.residues) ) )

            rmol = None
            rtype = r.type.upper()
            if rtype.lower() in rmolmap :
                rmol = rmolmap[ rtype.lower() ]
            else :
                rmol = GetResMol ( rtype.lower() )
                if rmol != None :
                    rmolmap[rtype.lower()] = rmol

            if rmol != None :
                for b in rmol.bonds :
                    a1n, a2n = b.atoms[0].name, b.atoms[1].name
                    if a1n in r.atomsMap and a2n in r.atomsMap :
                        for a1 in r.atomsMap[a1n] :
                            for a2 in r.atomsMap[a2n] :
                                #print "%s - %s" % ( At(a1), At(a2) )
                                nb = chMol.newBond ( a1, a2 )
                                pass
            else :
                print " - rmol %s not found" % rtype

            if 1 :
                if r.type.upper() in protein3to1 :
                    if r.id.position-1 in rmap :
                        pres = rmap[r.id.position-1]
                        if pres.type.upper() in protein3to1 :
                            #GetSS ( pres, r )
                            if "C" in pres.atomsMap and "N" in r.atomsMap :
                                for a1 in pres.atomsMap["C"] :
                                    for a2 in r.atomsMap["N"] :
                                        #print a1.name, pres.id.position, a2.name, r.id.position
                                        nb = chMol.newBond ( a1, a2 )
                                        pass

                if r.type.upper() in nucleic1to3 or r.type.upper() in nucleic3to1 :
                    if r.id.position-1 in rmap :
                        pres = rmap[r.id.position-1]
                        if pres.type.upper() in nucleic1to3 or pres.type.upper() in nucleic3to1 :
                            if "O3'" in pres.atomsMap and "P" in r.atomsMap :
                                for a1 in pres.atomsMap["O3'"] :
                                    for a2 in r.atomsMap["P"] :
                                        #print a1.name, pres.id.position, a2.name, r.id.position
                                        nb = chMol.newBond ( a1, a2 )
                                        pass

        print ( ", %d bonds, %.1fs" % (len(chMol.bonds), time.time()-start)  )

        #start = time.time()
        #chimera.openModels.add ( [chMol] )
        #print " - added mol %ss, %.1fs" % (chMol.name, time.time()-start)

        start = time.time()
        for r in chMol.residues :
            rt = r.type.upper()
            if rt in protein3to1 or rt in nucleic3to1 or rt in nucleic1to3 :
                r.ribbonDisplay = True
                r.ribbonDrawMode = 2
                r.ribbonColor = mol.chainColors[r.id.chainId]
                for at in r.atoms :
                    at.display = False
                    at.drawMode = at.EndCap
                    if at.element.name.upper() in atomColors :
                        at.color = atomColors[at.element.name.upper()]
                    else :
                        at.color = mol.chainColors[r.id.chainId]
            else :
                for at in r.atoms :
                    at.display = True
                    at.drawMode = at.EndCap
                    if at.element.name.upper() in atomColors :
                        at.color = atomColors[at.element.name.upper()]
                    else :
                        at.color = mol.chainColors[r.id.chainId]

        for b in mol.bonds :
            b.drawMode = b.Stick
            b.display = b.Smart

        #print " - changed mol disp %s - %.1fs" % (chMol.name, time.time()-start)

    return mol



def GetSS ( pres, r ) :

    if not "N" in r.atomsMap or not "CA" in r.atomsMap or not "C" in r.atomsMap :
        return

    n2 = r.atomsMap["N"][0]
    ca2 = r.atomsMap["CA"][0]
    c2 = r.atomsMap["C"][0]

    if not "N" in pres.atomsMap or not "CA" in pres.atomsMap or not "C" in pres.atomsMap :
        return

    n1 = pres.atomsMap["N"][0]
    ca1 = pres.atomsMap["CA"][0]
    c1 = pres.atomsMap["C"][0]
    #o1 = pres.atomsMap["O"][0]

    phi = diha ( n1, ca1, c1, n2 )
    om = diha ( ca1, c1, n2, ca2 )
    psi = diha ( c1, n2, ca2, c2 )
    #oo = diha ( n1, ca1, c1, o1 )

    r.isHelix = phi > -90.0 and phi < -30.0 and psi > -80.0 and psi < -20.0
    r.isHet = False
    r.isSheet = False
    r.isStrand = phi > -150.0 and phi < -90.0 and psi > 90 and psi < -150.0



def diha ( a1, a2, a3, a4 ) :
    #n1 = vnorm ( a1.coord(), a2.coord(), a3.coord() )
    #n2 = vnorm ( a2.coord(), a3.coord(), a4.coord() )
    #return numpy.arccos ( n2 * n1 * -1.0 ) * 180.0 / numpy.pi

    # http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = a2.coord() - a1.coord()
    b2 = a3.coord() - a2.coord()
    b3 = a4.coord() - a3.coord()

    n1 = chimera.cross ( b1, b2 ); n1.normalize()
    n2 = chimera.cross ( b2, b3 ); n2.normalize()
    m1 = chimera.cross ( n1, b2 ); m1.normalize()

    x = n1 * n2
    y = m1 * n2

    return -1.0 * numpy.arctan2 ( y, x) * 180.0 / numpy.pi




def ReadMol ( fpath, log=False ) :

    from random import random

    res = ReadCif ( fpath, log )
    if res == None :
        print " - could not read cif"
        return

    cif, loops = res

    # descriptions by chain id:
    descrByEntityId = GetEntityDescr ( cif, loops )

    try :
        atoms = loops['_atom_site']['data']
        print " - %d atom records" % len(atoms)
    except :
        print " - no atoms in cif?"
        return None

    labels = loops['_atom_site']['labels']
    if 0 :
        print "Labels:"
        for l in labels :
            print " : ", l


    import time
    start = time.time()

    rmap = {}

    nmol = chimera.Molecule()
    from os import path
    #nmol.name = path.splitext ( path.split (fpath)[1] )[0]
    nmol.name = path.split (fpath) [1]
    nmol.openedAs = [ fpath, [] ]
    nmol.cif = cif
    nmol.cifLoops = loops

    nmol.chainColors = {}
    nmol.chainDescr = {}

    numQ = 0
    first = True
    for at in atoms :
        mp = at['asMap']

        if log and first :
            #for label, val in mp.iteritems () :
            for li, label in enumerate ( labels ) :
                print "   %d : %s : %s" % (li+1, label, mp[label])

        first = False

        atType = mp['type_symbol']
        atName = mp['label_atom_id']
        rtype = mp['label_comp_id']
        chainId = mp['label_asym_id']
        chainEId = mp['label_entity_id']
        px = mp['Cartn_x']
        py = mp['Cartn_y']
        pz = mp['Cartn_z']
        occ = mp['occupancy']
        bfactor = mp['B_iso_or_equiv']
        altLoc = mp['label_alt_id']
        if altLoc == "." : altLoc = ''

        if chainEId in descrByEntityId :
            nmol.chainDescr [chainId] = descrByEntityId [chainEId]

        resId = ResId ( mp )
        if resId == None :
            continue

        ris = "%s%d" % (chainId, resId)
        res = None
        if not ris in rmap :
            res = nmol.newResidue ( rtype, chimera.MolResId(chainId, resId) )
            rmap[ris] = res
            res.chainEId = chainEId
        else :
            res = rmap[ris]

        clr = None
        if not chainId in nmol.chainColors :
            clr = chimera.MaterialColor ( random(), random(), random(), 1.0 )
            nmol.chainColors[chainId] = clr
            if 0 and log :
                print " - chain %s" % chainId
        else :
            clr = nmol.chainColors [chainId]

        nat = nmol.newAtom ( atName, chimera.Element(atType) )

        drawRib = rtype in protein3to1 or rtype in nucleic3to1

        #aMap[at] = nat
        res.addAtom( nat )
        nat.setCoord ( chimera.Point( float(px), float(py), float(pz) ) )
        nat.altLoc = altLoc
        nat.occupancy = float(occ)
        nat.bfactor = float(bfactor)

        if 'Q-score' in mp :
            try :
                Q = float ( mp['Q-score'] )
                nat.Q = Q
                numQ += 1
            except :
                #print " - q score is",  mp['Q-score']
                pass

    end = time.time()
    print " - created mol with %d atoms, %.1fs, %d q-scores" % ( len(nmol.atoms), end-start, numQ )

    return nmol


def DataStr ( data ) :
    if type(data) == dict :
        return data['string']
    return data

def GetEntityDescr ( cif, loops ) :

    descrByEntityId = {}
    if  '_entity' in loops :
        elabels, entities = loops['_entity']['labels'], loops['_entity']['data']
        print " - found _entity - %d labels, %d records" % ( len(elabels), len(entities) )
        for ent in entities :
            entMap = ent['asMap']

            descr = []
            #if 'type' in entMap :
            #    #descr.append ( "Type: " + entMap['type'] )
            #    descr.append ( DataStr(entMap['type']) )

            if 'pdbx_description' in entMap :
                #descr.append ( "Descr: " + entMap['pdbx_description'] )
                descr.append ( DataStr(entMap['pdbx_description']) )

            if 'id' in entMap :
                descrByEntityId[ entMap['id'] ] = descr

    if  '_entity_name_com' in loops :
        entities = loops['_entity_name_com']['data']
        print " - found _entity_name_com - %d records" % ( len(entities) )
        for ent in entities :
            entMap = ent['asMap']

            descr = []
            if 'entity_id' in entMap :
                eid = entMap['entity_id']
                if eid in descrByEntityId :
                    descr = descrByEntityId[eid]
                else :
                    descrByEntityId[eid] = descr

            if 'name' in entMap :
                #descr.append ( "Name: " + entMap['name'] )
                descr.append ( DataStr(entMap['name']) )

    return descrByEntityId


def ResId ( mp ) :

    resId = mp['label_seq_id']
    try :
        resId = int(resId)
    except :
        resId = None

    if resId == None :
        try :
            resId = int( mp['auth_seq_id'] )
        except :
            print " - atom resId not numeric: %s/%s" % ( mp['label_seq_id'], mp['auth_seq_id'] )
            resId = None

    return resId


def ConnectAtoms () :

    from chimera.resCode import nucleic3to1
    from chimera.resCode import protein3to1
    protein3to1['HSD'] = protein3to1['HIS']



def UpdateAtoms ( cif, mol, dmap=None ) :

    print "Updating atoms in cif - %s" % mol.name

    addQ = False
    amap = {}
    for r in mol.residues :
        for at in r.atoms :
            atId = "%d.%s.%s.%s" % (r.id.position,r.id.chainId,at.name,at.altLoc)
            amap[atId] = at
            if not addQ and hasattr ( at, 'Q' ) :
                addQ = True

    #print " - %d items in cif" % len(cif)

    for ls in cif :
        if type(ls) == list :
            name, labels, data = ls

            if name == "_atom_site" :
                #print " - found atoms - %d" % len(data)

                ilabels = {}
                for i, l in enumerate ( labels ) :
                    ilabels [ l ] = i
                    #print " -- %s - %d" % (l, i)

                addQatI = None
                if addQ :
                    if not 'Q-score' in ilabels :
                        if 'B_iso_or_equiv' in ilabels :
                            addQatI = ilabels['B_iso_or_equiv'] + 1
                        else :
                            addQatI = len (labels)
                        labels.insert ( addQatI, "Q-score" )
                        print " - added Q-score column %d" % addQatI
                        ilabels = {}
                        for i, l in enumerate ( labels ) :
                            ilabels [ l ] = i
                            #print " -- %s - %d" % (l, i)

                deli = []
                for di, d in enumerate ( data ) :
                    adata, mp = d['asArray'], d['asMap']

                    resId = ResId ( mp )
                    if resId == None : continue
                    atName = mp['label_atom_id']
                    chainId = mp['label_asym_id']
                    altLoc = mp['label_alt_id']
                    if altLoc == "." : altLoc = ''

                    datId = "%d.%s.%s.%s" % (resId,chainId,atName,altLoc)
                    if datId in amap :
                        at = amap[datId]
                        #adata[ ilabels['B_iso_or_equiv'] ] = "%.3f" % at.bfactor

                        if addQ :
                            qs = ("%.3f" % at.Q) if hasattr ( at, 'Q' ) else "?"
                            if addQatI :
                                adata.insert ( addQatI, qs )
                            else :
                                adata[ ilabels['Q-score'] ] = qs

                        C = at.coord()
                        if dmap != None :
                            C = dmap.openState.xform.inverse().apply ( at.xformCoord() )
                        adata[ ilabels['Cartn_x'] ] = "%.3f" % C.x
                        adata[ ilabels['Cartn_y'] ] = "%.3f" % C.y
                        adata[ ilabels['Cartn_z'] ] = "%.3f" % C.z

                    else :
                        #print " - atom %s in cif - not found in mol" % datId
                        deli.append ( di )

                if len(deli) > 0 :
                    deli.sort()
                    deli.reverse()
                    for di in deli :
                        del data[di]



def WriteMol ( mol, fout, dmap = None ) :

    if 0 and hasattr (mol, 'cif') :
        UpdateAtoms ( mol.cif, mol, dmap )
        WriteCif ( mol.cif, fout )

    else :
        print " ---------- making cif for %s ---------- " % mol.name
        #return
        mol.openedAs = [ fout, [] ]
        mol.cif = []
        mol.cifLoops = {}

        from os.path import splitext
        mol.cif.append ( "data_%s\n" % splitext(mol.name)[0] )
        mol.cif.append ( "#\n" )
        mol.cif.append ( "_entry.id %s\n" % splitext(mol.name)[0] )
        mol.cif.append ( "#\n" )

        AddAtoms ( mol.cif, mol.cifLoops, mol, dmap )
        if 1 :
            AddResDisplay ( mol.cif, mol.cifLoops, mol )
        WriteCif ( mol.cif, fout )



def AddAtoms ( cif, cifLoops, mol, dmap = None ) :

    name = "_atom_site"

    labels = []
    labels.append ( "group_PDB" )
    labels.append ( "id" )
    labels.append ( "type_symbol" )
    labels.append ( "label_atom_id" )
    labels.append ( "label_alt_id" )
    labels.append ( "label_comp_id" )
    labels.append ( "label_asym_id" )
    labels.append ( "label_entity_id" )
    labels.append ( "label_seq_id" )
    labels.append ( "pdbx_PDB_ins_code" )
    labels.append ( "Cartn_x" )
    labels.append ( "Cartn_y" )
    labels.append ( "Cartn_z" )
    labels.append ( "occupancy" )
    labels.append ( "B_iso_or_equiv" )
    labels.append ( "Q-score" )
    labels.append ( "pdbx_formal_charge" )
    labels.append ( "auth_seq_id" )
    labels.append ( "auth_comp_id" )
    labels.append ( "auth_asym_id" )
    labels.append ( "auth_atom_id" )
    labels.append ( "pdbx_PDB_model_num" )

    data = []

    cress = {}
    for r in mol.residues :
        if r.id.chainId in cress :
            cress[r.id.chainId].append ( [r.id.position, r] )
        else :
            cress[r.id.chainId] = [ [r.id.position, r] ]

    cids = cress.keys()
    print "- adding %d chains" % len(cids)
    cids.sort()

    chainEIds = {}
    for i, cid in enumerate ( cids ) :
        chainEIds[cid] = "%d" % (i+1)

    ati = 1

    for cid in cids :
        ress = cress[cid]
        ress.sort()
        entityId = chainEIds [cid]
        print ". %s - %d res - entity %s" % (cid, len(ress), entityId)

        for ri, r in ress :

            for at in r.atoms :

                aname = at.name
                #if '"' in aname : aname = "'" + aname + "'"
                #elif "'" in aname or " " in aname : aname = '"' + aname + '"'

                apos = at.coord()
                if dmap :
                    apos = dmap.openState.xform.inverse().apply ( at.xformCoord() )

                adata, mdata = [None] * len(labels), {}

                qstr = "?"
                if hasattr ( at, "Q" ) :
                    qstr = "%.3f" % at.Q

                mdata["group_PDB"]          = adata[0] = "ATOM"
                mdata["id"]                 = adata[1] = "%d" % ati
                mdata["type_symbol"]        = adata[2] = at.element.name
                mdata["label_atom_id"]      = adata[3] = aname
                mdata["label_alt_id"]       = adata[4] = "." if len(at.altLoc) == 0 else at.altLoc
                mdata["label_comp_id"]      = adata[5] = r.type
                mdata["label_asym_id"]      = adata[6] = r.id.chainId
                mdata["label_entity_id"]    = adata[7] = entityId
                mdata["label_seq_id"]       = adata[8] = "%d" % r.id.position
                mdata["pdbx_PDB_ins_code"]  = adata[9] = "?"
                mdata["Cartn_x"]            = adata[10] = "%.3f" % apos.x
                mdata["Cartn_y"]            = adata[11] = "%.3f" % apos.y
                mdata["Cartn_z"]            = adata[12] = "%.3f" % apos.z
                mdata["occupancy"]          = adata[13] = "%.3f" % at.occupancy
                if hasattr ( at, 'bfactor0' ) :
                    mdata["B_iso_or_equiv"]     = adata[14] = "%.3f" % at.bfactor0
                else :
                    mdata["B_iso_or_equiv"]     = adata[14] = "%.3f" % at.bfactor
                mdata["Q-score"]            = adata[15] = qstr
                mdata["pdbx_formal_charge"] = adata[16] = "?"
                mdata["auth_seq_id"]        = adata[17] = "%d" % r.id.position
                mdata["auth_comp_id"]       = adata[18] = r.type
                mdata["auth_asym_id"]       = adata[19] = r.id.chainId
                mdata["auth_atom_id"]       = adata[20] = aname
                mdata["pdbx_PDB_model_num"] = adata[21] = "1"

                data.append ( {'asArray':adata, 'asMap':mdata} )
                ati += 1

    cif.append ( [name, labels, data] )
    cif.append ( "#\n" )






def AddResDisplay ( cif, cifLoops, mol ) :

    name = "_mapq_comp_display_attributes"

    labels = []
    labels.append ( "label_asym_id" ) # chain id
    labels.append ( "label_seq_id" ) # position
    labels.append ( "ribbonDisplay" )
    labels.append ( "ribbonColor_R" )
    labels.append ( "ribbonColor_G" )
    labels.append ( "ribbonColor_B" )
    labels.append ( "ribbonColor_A" )

    data = []

    cress = {}
    for r in mol.residues :
        if r.id.chainId in cress :
            cress[r.id.chainId].append ( [r.id.position, r] )
        else :
            cress[r.id.chainId] = [ [r.id.position, r] ]

    print " - adding res display, %d residues" % len(mol.residues)

    cids = cress.keys()
    cids.sort()

    for cid in cids :
        ress = cress[cid]
        ress.sort()

        for ri, r in ress :

            if hasattr ( r, 'ribbonDisplay' ) and hasattr ( r, 'ribbonColor' ) and r.ribbonColor != None :
                adata, mdata = [None] * len(labels), {}
                C = r.ribbonColor.rgba()

                #mdata["label_comp_id"]     = adata[5] = r.type
                mdata["label_asym_id"]      = adata[0] = r.id.chainId
                mdata["label_seq_id"]       = adata[1] = "%d" % r.id.position
                mdata["ribbonDisplay"]      = adata[2] = "%d" % (1 if r.ribbonDisplay else 0)
                mdata["ribbonColor_R"]      = adata[3] = "%.3f" % C[0]
                mdata["ribbonColor_G"]      = adata[4] = "%.3f" % C[1]
                mdata["ribbonColor_B"]      = adata[5] = "%.3f" % C[2]
                mdata["ribbonColor_A"]      = adata[6] = "%.3f" % C[3]

                data.append ( {'asArray':adata, 'asMap':mdata} )

    cif.append ( [name, labels, data] )
    cif.append ( "#\n" )


#
