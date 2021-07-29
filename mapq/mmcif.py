



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
            print " - ? %d -- %s" % (li, ls)

    fp.close()

    print ( " - read %d lines, %.1fs" % (li, time.time()-start) )

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
        #print " _ %d %s" % (li, ls)

        if len(ls) == 0 :
            continue

        elif ls[0] == "_" :
            if getData :
                print " - done loop?"
                getNext = False
                break
            else :
                # loop label
                name, label = ls.split(".")
                loopLabels.append ( label )
                loopName = name

        elif ls[0] == "#" :
            #print " - done loop"
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
        if ls[0] == ';' :
            data = ls[1:]

            # look for end ;
            while 1 :
                atLine, getNext = fp.readline(), True
                if not atLine :
                    print " - ? %d - eof while getting single value starting at %d" % (li, liStart)
                    return li, name, data, lines

                li += 1
                lines += atLine
                ls = atLine.strip ()
                if ls[0] == ';' :
                    # done
                    break
                else :
                    data += " " + ls
        else :
            data = ls

    else :
        print " - ? %d - " % li, ls

    return li, name, data, lines


def GetData ( fp, atLine, ls, li, labels ) :
    # get (multiple) data for loop

    data = []
    liStart = li

    while 1 :
        if ls[0] == ';' :
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
    for t in ts :
        if addTo :
            addTo += " " + t
            if t[-1] == "'" :
                tsr.append ( addTo[0:-1] )
                addTo = None
        elif t[0] == "'" :
            if t[-1] == "'" :
                tsr.append ( t[1:-1] )
            else :
                addTo = t[1:]
        else :
            tsr.append ( t )
    if addTo :
        print " - ? %d - unmatched '" % li
        tsr += addTo
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

            cws = [0] * len(labels)
            for d in data :
                adata, mdata = d['asArray'], d['asMap']

                for i, ds in enumerate(adata) :
                    dd = "'%s'" % ds if ' ' in ds else ds
                    cws[i] = max ( cws[i], len(dd) )

            for d in data :
                adata, mdata = d['asArray'], d['asMap']
                first = True
                for di, dd in enumerate(adata) :
                    if type(dd) == dict :
                        # write original block starting with ; on start and end lines
                        if not first :
                            fp.write ( "\n" )
                        fp.write ( dd['lines'] )
                        first = True
                    else :
                        # write in columns
                        dd = "'%s'" % dd if ' ' in dd else dd
                        #dd = dd if first else ("\t" + dd)
                        padn = cws[di] - len(dd) + 1
                        fp.write ( "%s%s" % (dd, " "*padn) )
                        first = False
                fp.write ("\n")
        else :
            fp.write ( ls )

    fp.close()



def LoadMol ( fpath, log=False ) :

    import chimera
    mol = ReadMol ( fpath, log )
    chimera.openModels.add ( [mol] )
    #return mol


def ReadMol ( fpath, log=False ) :

    from random import random
    import chimera

    from chimera.resCode import nucleic3to1
    from chimera.resCode import protein3to1, protein1to3
    protein3to1['HSD'] = protein3to1['HIS']
    protein3to1['HSE'] = protein3to1['HIS']


    cif, loops = ReadCif ( fpath, log )

    descrByEntityId = GetEntityDescr ( cif, loops )

    atoms = loops['_atom_site']['data']
    print " - %d atom records" % len(atoms)

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
        nat.drawMode = nat.Sphere
        nat.color = clr
        nat.display = True # not drawRib
        nat.altLoc = altLoc
        nat.occupancy = float(occ)
        nat.bfactor = float(bfactor)

        if 'Q_score' in mp :
            try :
                Q = float ( mp['Q_score'] )
                nat.Q = Q
            except :
                pass

        #res.isHelix = res.isHelix
        #res.isHet = res.isHet
        #res.isSheet = res.isSheet
        #res.isStrand = res.isStrand
        res.ribbonDisplay = False # drawRib
        res.ribbonDrawMode = 2
        res.ribbonColor = clr

    #for bond in mol.bonds :
    #    nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
    #    nb.display = nb.Smart

    print " - created %d atoms, %.1fs" % ( len(nmol.atoms), time.time()-start )

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



def UpdateAtoms ( cif, mol ) :

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
                if addQ and not 'Q_score' in ilabels :
                    if 'B_iso_or_equiv' in ilabels :
                        addQatI = ilabels['B_iso_or_equiv'] + 1
                    else :
                        addQatI = len (labels)
                    labels.insert ( addQatI, "Q_score" )
                    print " - added Q_score column %d" % addQatI

                ilabels = {}
                for i, l in enumerate ( labels ) :
                    ilabels [ l ] = i
                    #print " -- %s - %d" % (l, i)

                for d in data :
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
                        qs = ("%.3f" % at.Q) if hasattr ( at, 'Q' ) else "?"
                        if addQatI :
                            adata.insert ( addQatI, qs )
                        else :
                            adata[ ilabels['Q_score'] ] = qs
                    else :
                        print " - atom %s in cif - not found in mol" % datId



def WriteMol ( mol, fout ) :

    if not hasattr (mol, 'cif') :
        print " - cif not found in %s" % mol.name
        return

    UpdateAtoms ( mol.cif, mol )
    WriteCif ( mol.cif, fout )













#
