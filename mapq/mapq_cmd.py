
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


import sys, os

chimeraPath = None
numProc = 1
res = 3.0
bfactor = -1
gSigma = 0.6

pdbs = []
cifs = []
maps = []

mapqVersion = "1.9.11"

#scriptPath = os.path.dirname(os.path.realpath(__file__))
#fp = os.open ( os.path.join(scriptPath, "mapq.py" ) )

print ("")
print ("Command Line Script - MapQ Version %s" % mapqVersion)
print ("")
print ("Found parameters:")

ok = True

for arg in sys.argv :

    print (": %s" % arg),

    if "mapq_cmd.py" in arg :
        print ( " -> this script" )

    elif os.path.isdir(arg) :
        #print ( " -> Chimera path" )

        arg = os.path.abspath(arg)

        cp = os.path.join ( os.path.join ( arg, 'bin' ), 'chimera' )
        if os.path.isfile(cp) :
            print ( " -> Chimera path, Unix" )
            chimeraPath = cp

        cp = os.path.join ( os.path.join ( arg, 'bin' ), 'chimera.exe' )
        if os.path.isfile(cp) :
            print ( " -> Chimera path, Windows" )
            chimeraPath = cp

        cp = os.path.join ( os.path.join ( os.path.join ( arg, 'Contents' ), 'MacOS' ), "chimera" )
        if os.path.isfile(cp) :
            print ( " -> Chimera path, Mac" )
            chimeraPath = cp


    elif os.path.isfile( arg ) :

        arg = os.path.abspath(arg)

        if ".pdb" in arg or ".ent" in arg :
            print ( " -> PDB model" )
            pdbs.append ( arg )
        elif ".cif" in arg :
            print ( " -> mmCIF model" )
            cifs.append ( arg )
        elif ".mrc" in arg or ".map" in arg :
            print ( " -> map" )
            maps.append ( arg )
        else :
            print ( " -X- extension not recognized" )
            print ( " --- use map=[] or pdb=[] or cif=[] to specify type of file" )
            ok = False

    elif arg[0:2] == "-v" or arg[0:len('--version')] == "--version" :
        print ("\n\nCommand Line Script - MapQ Version: %s\n\n" % mapqVersion)
        sys.exit(0)

    else :
        tokens = arg.split("=")
        if len(tokens) == 2 and tokens[0] == "np" :
            try :
                numProc = int ( tokens[1] )
                print ( " -> number of processes: %d" % numProc )
            except :
                print ( " -> specify an integer" )
        elif len(tokens) == 2 and tokens[0] == "res" :
            try :
                res = float ( tokens[1] )
                print ( " -> resolution: %.3f" % res )
            except :
                print ( " -> specify a number" )
        elif len(tokens) == 2 and tokens[0] == "bfactor" :
            try :
                bfactor = float ( tokens[1] )
                print ( " -> bfactor: %.0f" % bfactor )
            except :
                print ( " -> specify a number" )
        elif len(tokens) == 2 and tokens[0] == "sigma" :
            try :
                gSigma = float ( tokens[1] )
                print ( " -> sigma: %.1f" % gSigma )
            except :
                print ( " -> specify a number" )

        elif len(tokens) == 2 and tokens[0] == "map" :
            print ( " -> map" )
            fpath = os.path.abspath( tokens[1] )
            if os.path.isfile( fpath ) :
                maps.append ( fpath )
            else :
                print ( " -X- map - file not found" )
                ok = False

        elif len(tokens) == 2 and tokens[0] == "pdb" :
            print ( " -> pdb" )
            fpath = os.path.abspath( tokens[1] )
            if os.path.isfile( fpath ) :
                pdbs.append ( fpath )
            else :
                print ( " -X- pdb - file not found" )
                ok = False

        elif len(tokens) == 2 and tokens[0] == "cif" :
            print ( " -> cif" )
            fpath = os.path.abspath( tokens[1] )
            if os.path.isfile( fpath ) :
                cifs.append ( fpath )
            else :
                print ( " -X- cif - file not found" )

        else :
            print ( " -> file/path not found or parameter not recognized" )

print ("")



mods = pdbs+cifs
if len(mods) == 0 :
    print ( " -X- need at least one model" )
    ok = False
if len(maps) != 1 :
    print ( " -X- need one map (only)" )
    ok = False
if chimeraPath == None :
    print ( " -X- please specify a path to Chimera" )
    ok = False

if not ok :
    print ("")
    print ("mapq_cmd.py")
    print ("  - Calculate Q-scores from command line")
    print ("")
    print ("Parameters:")
    print ("  [path to model or map file]")
    print ("      specify a map or model - will try to autodetect type based on extension")
    print ("  [path to Chimera]")
    print ("      e.g.: ~/Chimera.app (Mac)")
    print ("          ~/Chimera (Unix)")
    print ("          C:\\Users\\name\\Chimera (Windows)")
    print ("  sigma=# (optional)")
    print ("      sigma of reference Gaussian, default is 0.6")
    print ("  res=# (optional)")
    print ("      resolution of map, e.g. res=3.2, default=3.0")
    print ("      only used in output of Q-scores/residue as comparison")
    print ("  bfactor=f (optional, f=50,100,200,...)")
    print ("      if specified, Q-scores are converted to Bfactors")
    print ("      using the formula bfactor=f*(1.0-Qscore)")
    print ("  np=# (optional, #=1,2,3,4,..., default=1")
    print ("      number of processors to use")
    print ("  pdb=<path to PDB file>")
    print ("      specify a model with PDB format")
    print ("  cif=<path to mmCIF file>")
    print ("      specify a model with mmCIF format")
    print ("  map=<path to map (ccp4 format) file>")
    print ("      specify a MAP with ccp4 format")

print ("")

if ok :

    #scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = os.path.dirname ( maps[0] )
    newScript = os.path.join ( scriptPath, "_mapqScript.py" )

    print ("Creating Chimera script in %s" % newScript)
    print ("")

    try :
        fp = open ( newScript, "w" )
    except :
        print ( " - could not write script for Chimera, check if you have write permission in %s" % scriptPath )
        exit ( 0 )

    fp.write ( "import chimera\n" )
    fp.write ( "import VolumeViewer\n" )
    fp.write ( "try:\n" )
    fp.write ( "    from mapq import qscores\n" )
    fp.write ( "    from mapq import mmcif\n" )
    fp.write ( "except:\n" )
    fp.write ( "    from Segger import qscores\n" )
    fp.write ( "    from Segger import mmcif\n" )


    print ("Running:")
    cmd = "%s --nogui --silent --nostatus " % chimeraPath

    for mod in pdbs :
        fp.write ( "\n" )
        fp.write ( "dmap = VolumeViewer.open_volume_file ( '%s', 'ccp4')[0]\n" % maps[0] )
        fp.write ( "mol = chimera.openModels.open ( '%s', type='PDB' )[0]\n" % mod )
        fp.write ( "qscores.Calc('%s', mol, %d, %f, %f, %f)\n" % (chimeraPath, numProc, res, bfactor, gSigma) )
        fp.write ( "chimera.openModels.close ( [dmap,mol] )\n" )

    for mod in cifs :
        fp.write ( "\n" )
        fp.write ( "dmap = VolumeViewer.open_volume_file ( '%s', 'ccp4')[0]\n" % maps[0] )
        fp.write ( "mol = mmcif.ReadMol ( '%s' )\n" % mod )
        fp.write ( "qscores.Calc('%s', mol, %d, %f, %f, %f)\n" % (chimeraPath, numProc, res, bfactor, gSigma) )
        fp.write ( "chimera.openModels.close ( [dmap,mol] )\n" )

    cmd += "'%s'" % newScript

    print (" : " + cmd)
    print ("")

    if 0 :
        fp.write ( "qscores.Calc('%s', mol, %d, %f, %f, %f)\n" % (chimeraPath, numProc, res, bfactor, gSigma) )

    elif 0 :
        if numProc == 1 :
            #CalcQ ( mol, None, dmap, sigma, log=True )
            fp.write ( "qscores.CalcQ ( mol, 'All', dmap, %f )\n" % ( gSigma) )
        else :
            #CalcQp ( mol, None, dmap, sigma, numProc=numProc, chimeraPath=chimeraPath )
            #CalcQpn ( mol, None, dmap, sigma, numProc=numProc, chimeraPath=chimeraPath, closeMap=True )
            fp.write ( "qscores.CalcQpn ( mol, 'All', dmap, %f, useOld=False, numProc=%d, chimeraPath='%s', closeMap=True )\n" % ( gSigma, numProc, chimeraPath) )



    fp.close()

    os.system(cmd)

    if 0 :
        print ( "Removing temp Chimera script ")
        print ( " - %s" % newScript )
        os.remove ( newScript )
        print ( " - %s" % (newScript+"c") )
        os.remove ( newScript + "c" )
