


import os
import shutil
import sys

print sys.argv[1]
fromPath = sys.argv[1]

files = ["mapq.py", "qscores.py"]

for f in files :
    print f,

    try :
        shutil.copy2 ( fromPath + "/" + f, "./" + f )
        print "  - ok"
    except :
        print "?"
