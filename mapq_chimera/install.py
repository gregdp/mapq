import sys, os, shutil

if len(sys.argv) != 2 :
    print ("")
    print ("Please add the path where Chimera is installed, e.g.:")
    print ("   python install.py /home/greg/applications/Chimera")
    print ("")

    sys.exit(0)


print ("")



sharePaths = []
def FindShare (path, lev) :

    try :
        fs = os.listdir (path)
    except :
        fs = []

    for f in fs :
        atPath = os.path.join ( path, f )
        if os.path.isdir ( atPath ) :
            #for i in range(lev+1) :
            #    print "  ",
            #print f
            if f == "share" and "chimera" in atPath.lower() :
                #print " - found: ", atPath
                sharePaths.append ( atPath )
                return
            else :
                FindShare ( atPath, lev+1 )






# Mac...
opath1 = os.path.join ( sys.argv[1], "Contents" )
opath1 = os.path.join ( opath1, "Resources" )
opath1 = os.path.join ( opath1, "share" )

# Unix...
opath2 = os.path.join ( sys.argv[1], "share" )

didInstall = False

chimeraPath = sys.argv[1]
if not os.path.isdir ( chimeraPath ) :
    print ("")
    print ("The specified Chimera path '%s' doesn't seem to exist" % chimeraPath )
    print (" - please check and try again")
    print ("")
    sys.exit(0)

if not "chimera" in chimeraPath.lower() :
    print ("")
    print ("The specified path '%s' doesn't seem to be for Chimera" % chimeraPath )
    print (" - please check and try again")
    print ("")
    sys.exit(0)


#print ("finding...")
FindShare (sys.argv[1], 0)
sharePath = None
for p in sharePaths :
    if sharePath == None or len(p) < len(sharePath) :
        sharePath = p
#print (sharePath)

#exit()

if sharePath == None :
    print ("")
    print ("Could not find the 'share' folder in the specified path" )
    print (" - please check that the path points to Chimera and try again")
    print ("")
    sys.exit(0)



for opath in [sharePath] :

    if os.path.isdir( opath ) :
        opath = os.path.join ( opath, "mapq" )

        if os.path.isdir( opath ) :
            print (" - removing previous MapQ:" + opath)
            try :
                shutil.rmtree(opath)
            except :
                pass

        #print " - copying from:", os.getcwd()
        print (" - copying . ->" + opath )

        try :
            shutil.copytree ( os.getcwd(), opath )
            didInstall = True
        except :
            print ("")
            print ("-----------------------------------------------------")
            print ("Problem: Could not copy to:", opath)
            print (" 1. please check if you have write access")
            print (" 2. try with sudo python install.py <path to Chimera>")
            print ("-----------------------------------------------------")
            print ("")
            break

        didInstall = True

if didInstall :

    print ("")
    print ("------------------------")
    print ("Installation successful.")
    print ("------------------------")
    print ("")
    print ("To use:")
    print (" 1. Restart Chimera.")
    print (" 2. Select Tools -> Volume Data -> MapQ")
    print ("")

    #wh = os.path.join ( os.getcwd(), "install.html" )
    #import webbrowser
    #webbrowser.open( 'file://' + wh, new=2)


else :
    print ("")
    print ("-----------------------------------------------------------------------")
    print ("Problem: Could not find 'share' folder in '" + sys.argv[1] + "'")
    print (" 1. please check the path")
    print (" 2. remember you can auto-complete while typing the path with <tab>")
    print ("-----------------------------------------------------------------------")
    print ("")
