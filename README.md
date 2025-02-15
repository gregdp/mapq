# mapq

A plugin for <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> to calculate and visualize <strong>Q-scores</strong> in 3D cryo-EM maps.


* <a href="https://github.com/gregdp/mapq/wiki/MapQ-Install">Install the latest version</a>
* [Tutorials](https://github.com/gregdp/mapq/tree/master/tutorials)
* [Report an issue or ask a question](https://github.com/gregdp/mapq/issues)

Practical Notes:
* Q-scores are now reported in EMDB validation reports. Expected_Q scores are not shown, but can be helpful to interpret the score; see tutorial for more details.
* EMDB reports Q-scores calculated with sigma=0.4, which are ~0.1 lower than Q-scores calculated with the default sigma=0.6; this is mainly to allow comparison across all cryo-EM maps in the EMDB which go up to a highest of ~1Å resolution.
* When using this plugin, Q-scores are usually compared to the expected_Q score. Sigma=0.4 should be used when the map is 1.5Å resolution and higher, otherwise either sigma is ok, and the expected_Q score will be reported for either sigma=0.4 or sigma=0.6. Using sigma=0.6 may be more beneficial when approaching ~4Å resolution maps.
* If using Google Drive: when calculating Q-scores with multiple processes, the process fails if the map and model files are on a Google Drive path. Fix is in the works.
* A video that shows how to genreate a color key for displaying Q-scores on a model ribbon: https://www.youtube.com/watch?v=lxy3reAXLKI
  
More details:
* (2020) Measurement of atom resolvability in cryo-EM maps with Q-scores <a href="https://www.nature.com/articles/s41592-020-0731-1" target="_blank">Nature Methods</a>, <a href="https://www.biorxiv.org/content/10.1101/722991v2" target="_blank">BioRXiv</a>
* (2020) Resolving individual atoms... <a href="https://www.nature.com/articles/s41422-020-00432-2">Cell Research</a>
* (2021) Validation, analysis and annotation of cryo-EM structures <a href="https://onlinelibrary.wiley.com/iucr/doi/10.1107/S2059798321006069">Acta Cryst. Sect. D.</a>
* (2022) Electron microscopy holdings of the Protein Data Bank: the impact of the resolution revolution, new validation tools, and implications for the future <a href="https://link.springer.com/article/10.1007/s12551-022-01013-w">Biophysical Reviews</a>
