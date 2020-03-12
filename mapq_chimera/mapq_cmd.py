


import sys, os

mods = []
chimeraPath = None
numProc = 1


for arg in sys.argv :

    print ("%s" % arg),

    if "mapq_cmd.py" in arg :
        print ( " -> this script" )

    elif os.path.isdir(arg) :
        #print ( " -> Chimera path" )

        cp = os.path.join ( os.path.join ( arg, 'bin' ), 'chimera' )
        if os.path.isfile(cp) :
            print ( " -> Chimera path, Unix" )
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
        else :
            print ( " -> unknown" )

print ("")

ok = True
if len(mods) < 2 :
    print (" - Please specify at least one map and one model. This script checks if the file exists and does not count it if not.")
    ok = False

if chimeraPath == None :
    print (" - Please specify path to Chimera. If specified, the script may not have found a valid path.")
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
    fp.write ( "mapq.mapq.Calc('%s',%d)\n" % (chimeraPath, numProc) )
    fp.close()


    print ("Running:")
    cmd = "%s --nogui --silent --nostatus " % chimeraPath
    for mod in mods :
        cmd += "%s " % mod
    cmd += "%s" % newScript

    print (" : " + cmd)
    print ("")

    os.system(cmd)
