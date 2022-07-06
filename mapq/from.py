


import os
import shutil
import sys

print sys.argv[1]
fromPath = sys.argv[1]

files = ["mapq.py", "qscores.py", "mmcif.py", "gridm.py", "mapq_cmd.py"]

for f in files :
    print f,

    try :
        shutil.copy2 ( fromPath + "/" + f, "./" + f )
        print "  - ok"
    except :
        print "?"


for f in os.listdir ( fromPath + "/_param" ) :

    fname, fext = os.path.splitext (f)
    if fext == ".pdb" :
        print f,
        shutil.copy2 ( fromPath + "/_param/" + f, "./_param/" + f )

print ""
