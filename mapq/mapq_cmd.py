
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

mods = []
chimeraPath = None
numProc = 1
res = 3.0
bfactor = -1
gSigma = 0.6


print ("")
print ("Found parameters:")

for arg in sys.argv :

    print (": %s" % arg),

    if "mapq_cmd.py" in arg :
        print ( " -> this script" )

    elif os.path.isdir(arg) :
        #print ( " -> Chimera path" )

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


    elif os.path.isfile(arg) :
        print ( " -> map or model" )
        if arg[0:2] == '..' :
            print ( " -X- please do not use relatives, i.e. .., in path (sorry)" )
        else :
            mods.append ( arg )

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
                print ( " -> sigma: %.0f" % gSigma )
            except :
                print ( " -> specify a number" )
        else :
            print ( " -> unknown" )

print ("")

ok = True
if len(mods) <= 1 or chimeraPath == None :
    print ("")
    print ("mapq_cmd.py")
    print ("  - Calculate Q-scores from command line")
    print ("")
    print ("Parameters:")
    print ("  [path to model or map file]")
    print ("      one map and at least one model should be specified")
    print ("  [path to Chimera]")
    print ("      e.g.: ~/Chimera.app (Mac)")
    print ("          ~/Chimera (Unix)")
    print ("          C:\\Users\\name\\Chimera (Windows)")
    print ("  sigma=# (optional)")
    print ("      sigma of reference Gaussian, default is 0.6")
    print ("  res=# (optional)")
    print ("      resolution of map, e.g. res=3.2")
    print ("      only used in output of Q-scores/residue as comparison")
    print ("  bfactor=f (optional, f=50,100,200,...)")
    print ("      if specified, Q-scores are converted to Bfactors")
    print ("      using the formula bfactor=f*(1.0-Qscore)")
    print ("  np=# (optional, #=1,2,3,4,...")
    print ("      number of processors to use")

    ok = False

print ("")

if ok :

    #scriptPath = os.path.dirname(os.path.realpath(__file__))
    scriptPath = os.path.dirname ( mods[0] )
    newScript = os.path.join ( scriptPath, "_mapqScript.py" )

    print ("Creating Chimera script in %s" % newScript)
    print ("")

    try :
        fp = open ( newScript, "w" )
    except :
        print ( " - could not write script for Chimera, check if you have write permission in %s" % scriptPath )
        exit ( 0 )

    fp.write ( "import mapq\n" )
    fp.write ( "import mapq.qscores\n" )
    fp.write ( "from mapq.mmcif import LoadMol as CifLoadMol\n" )


    print ("Running:")
    cmd = "%s --nogui --silent --nostatus " % chimeraPath
    for mod in mods :
        if os.path.splitext(mod)[1] == ".cif" :
            fp.write ( "CifLoadMol('%s')\n" % mod )
        else :
            cmd += '"%s" ' % mod

    cmd += "'%s'" % newScript

    print (" : " + cmd)
    print ("")

    fp.write ( "mapq.qscores.Calc('%s',%d,%f,%f,%f)\n" % (chimeraPath, numProc, res, bfactor, gSigma) )
    fp.close()

    os.system(cmd)

    if 1 :
        print ( "Removing temp Chimera script ")
        print ( " - %s" % newScript )
        os.remove ( newScript )
        print ( " - %s" % (newScript+"c") )
        os.remove ( newScript + "c" )
