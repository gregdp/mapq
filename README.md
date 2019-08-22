# mapq

Calculates Q-scores

To Install:

1. First, <a href="https://www.cgl.ucsf.edu/chimera/download.html">download</a> and install UCSF Chimera. (Run it once before installing the plugin; on some platforms, e.g. MacOS, you may see a warning message on first run which you have to accept. This may prevent further issues after adding the plugin.)
2. <a href="https://github.com/gregdp/mapq/tree/master/download">Download</a> latest version of MapQ.
3. In a terminal, navigate to where the file was downloaded, then execute the following commands:
* unzip mapq_chimera_1_2.zip
* cd mapq_chimera_1_2
* python install.py [path to Chimera]

e.g.:
* python install.py ~/Desktop/Chimera.app/


Note that on Windows, you may use the python bundled with Chimera itself, so the third command would be
* [path to Chimera]/bin/python install.py [path to Chimera]

To Run:
1. (Re)start Chimera*
2. Start MapQ: Tools -> Volume Data -> MapQ
3. See [tutorial](https://github.com/gregdp/mapq/tree/master/tutorials)
3. More details [here](https://cryoem.slac.stanford.edu/ncmi/resources/software/modelz)

\* On Mac OS, an error message may be shown on first run after installing, see [here](https://www.santoshsrinivas.com/disable-gatekeeper-in-macos-sierra/) for solution.

UCSF Chimera runs on MacOS, Linux and Windows. To install Chimera it should take ~10 minutes, and MapQ another ~5 minutes. For a tutorial on MapQ, please see the "tutorials" folder.
