


import sys, os

mods = []
chimeraPath = None
numProc = 1
res = 3.0
bfactor = -1


print ("")
print ("Found parameters:")

for arg in sys.argv :

    print ("  %s" % arg),

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
        mods.append ( arg )

    else :
        tokens = arg.split("=")
        if len(tokens) == 2 and tokens[0] == "np" :
            try :
                numProc = int ( tokens[1] )
                print ( " -> number of processes: %d" % numProc )
            except :
                print ( " -> unknown" )
        elif len(tokens) == 2 and tokens[0] == "res" :
            try :
                res = float ( tokens[1] )
                print ( " -> resolution: %.3f" % res )
            except :
                print ( " -> unknown" )
        elif len(tokens) == 2 and tokens[0] == "bfactor" :
            try :
                bfactor = float ( tokens[1] )
                print ( " -> bfactor: %.0f" % bfactor )
            except :
                print ( " -> unknown" )
        else :
            print ( " -> unknown" )

print ("")

ok = True
if len(mods) == 0 or chimeraPath == None :
    print ("")
    print ("mapq_cmd.py")
    print ("  - Calculate Q-scores from command line")
    print ("")
    print ("Parameters:")
    print ("  - [path to model or map file]")
    print ("    one map and at least one model should be specified")
    print ("    files must exist")
    print ("  - [path to Chimera]")
    print ("    e.g.: ~/Chimera.app (Mac)")
    print ("          ~/Chimera (Unix)")
    print ("          C:\\Users\\name\\Chimera (Windows)")
    print ("  - res=# (optional)")
    print ("    resolution of map, e.g. res=3.2")
    print ("    only used in output of Q-scores/residue as comparison")
    print ("  - bfactor=f (optional)")
    print ("    if specified, Q-scores are converted to Bfactors")
    print ("    using the formula bfactor=f*(1.0-Qscore)")
    ok = False

print ("")

if ok :

    scriptPath = os.path.dirname(os.path.realpath(__file__))
    newScript = os.path.join ( scriptPath, "mapqScript.py" )

    print ("Creating Chimera script in %s" % newScript)
    print ("")

    try :
        fp = open ( newScript, "w" )
    except :
        print ( " - could not write script for Chimera, check if you have write permission in %s" % scriptPath )
        exit ( 0 )

    fp.write ( "import mapq\n" )
    fp.write ( "import mapq.mapq\n" )
    fp.write ( "mapq.mapq.Calc('%s',%d,%f,%f)\n" % (chimeraPath, numProc, res, bfactor) )
    fp.close()


    print ("Running:")
    cmd = "%s --nogui --silent --nostatus " % chimeraPath
    for mod in mods :
        cmd += "%s " % mod
    cmd += "%s" % newScript

    print (" : " + cmd)
    print ("")

    os.system(cmd)
