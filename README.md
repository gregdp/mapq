# mapq

A plugin for <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> to calculate and visualize <strong>Q-scores</strong> in 3D cryo-EM maps.


* <a href="https://github.com/gregdp/mapq/wiki/MapQ-Install">Install the latest version</a>
* [Tutorials](https://github.com/gregdp/mapq/tree/master/tutorials)
* [Report an issue or ask a question](https://github.com/gregdp/mapq/issues)

Practical Notes:
* Q-scores are now reported in EMDB validation reports. Expected_Q scores are not shown, but can be helpful to interpret the score; see tutorial for more details. Also note that EMDB reports Q-scores calculated with sigma=0.4, which are ~0.1 lower than Q-scores calculated with the default sigma=0.6; this is mainly to allow comparison across all cryo-EM maps in the EMDB which go up to a highest of ~1Å resolution.
* When using this plugin, the Q-score is more practically compared to the expected_Q score, which will correspond correctly based on the chosen sigma. Use of sigma=0.4 is only really necessary when the map is 1.5Å resolution and higher.
* If using Google Drive: when calculating Q-scores with multiple processes, the process fails if the map and model files are on a Google Drive path. Please copy them to a non Google Drive location for this to work.
* A video that shows how to genreate a color key for displaying Q-scores on a model ribbon: https://www.youtube.com/watch?v=lxy3reAXLKI
* On Macs, new versions of Chimera are painfully laggy. Version 1.13 works much better overall, and Segger will work the same. It can be obtained from the <a href=https://www.cgl.ucsf.edu/chimera/olddownload.html>Old Releases</a> link on the Chimerea downloads page. Overdue but Q-scores are slowly being migrated to ChimeraX.

More details:
* (2020) Measurement of atom resolvability in cryo-EM maps with Q-scores <a href="https://www.nature.com/articles/s41592-020-0731-1" target="_blank">Nature Methods</a>, <a href="https://www.biorxiv.org/content/10.1101/722991v2" target="_blank">BioRXiv</a>
* (2020) Resolving individual atoms... <a href="https://www.nature.com/articles/s41422-020-00432-2">Cell Research</a>
* (2021) Validation, analysis and annotation of cryo-EM structures <a href="https://onlinelibrary.wiley.com/iucr/doi/10.1107/S2059798321006069">Acta Cryst. Sect. D.</a>
