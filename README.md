# lcolens
Time domain gravitational lens photometry of Las Cumbres Observatory observations

Order of things on AAS hack day, 1/8/15
- Steve: PSF-creation code within pyphot DOCUMENTED
- Curtis: distutils/setuptools configuration for pyphot DONE
- Leonidas: 
 - reading in the sequence of lcogt images and header information
 - import stellar catalog
  - Nb astropy affiliated package astroquery, for tapping into Vizier catalog
  - HIPPARCOS for astrometry
  - SDSS for photometry
  - NOMAD via Vizier, http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/297
 - WCS match of stellar catalog to positions in lcogt image
 - build PSF using pyphot pkfit routines
 - pygfit crowded-field PSF fitting setup for project, using pyphot object detection
 - convert pyphot's PSF description into galfit-format PSF for pygfit use

