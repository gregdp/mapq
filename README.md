# mapq

A plugin for <a href="https://www.cgl.ucsf.edu/chimera/">UCSF Chimera</a> to calculate and visualize <strong>Q-scores</strong> in 3D cryo-EM maps.


* <a href="https://github.com/gregdp/mapq/wiki/MapQ-Install">Install the latest version</a>
* [Tutorials](https://github.com/gregdp/mapq/tree/master/tutorials)
* [Report an issue or ask a question](https://github.com/gregdp/mapq/issues)

Practical Notes:
* Q-scores are now reported in EMDB validation reports, they are calculated using this plugin.
* The default sigma is now 0.4 to match Q-scores calculated in the EMDB.
* Q-scores can be compared to Q_peak, Q_low_95%, and Q_high_95%, which compares it to other maps and model in the EMDB at similar reported resolution. See the recent biorXiv for an explanation.
* If using Google Drive: when calculating Q-scores with multiple processes, the process fails if the map and model files are on a Google Drive path. 
* A video that shows how to genreate a color key for displaying Q-scores on a model ribbon: https://www.youtube.com/watch?v=lxy3reAXLKI

https://www.biorxiv.org/content/10.1101/2025.01.14.633006v1
  
More details:
* (2025) Q-score as a reliability measure for protein, nucleic acid, and small molecule atomic coordinate models derived from 3DEM density maps  <a href="https://www.biorxiv.org/content/10.1101/2025.01.14.633006v1" target="_blank">BioRXiv</a>
* (2020) Measurement of atom resolvability in cryo-EM maps with Q-scores <a href="https://www.nature.com/articles/s41592-020-0731-1" target="_blank">Nature Methods</a>, <a href="https://www.biorxiv.org/content/10.1101/722991v2" target="_blank">BioRXiv</a>
* (2020) Resolving individual atoms... <a href="https://www.nature.com/articles/s41422-020-00432-2">Cell Research</a>
* (2021) Validation, analysis and annotation of cryo-EM structures <a href="https://onlinelibrary.wiley.com/iucr/doi/10.1107/S2059798321006069">Acta Cryst. Sect. D.</a>
* (2022) Electron microscopy holdings of the Protein Data Bank: the impact of the resolution revolution, new validation tools, and implications for the future <a href="https://link.springer.com/article/10.1007/s12551-022-01013-w">Biophysical Reviews</a>
