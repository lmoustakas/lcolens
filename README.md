# lcolens
Time domain gravitational lens photometry of Las Cumbres Observatory observations

Order of things on AAS hack day, 1/8/15
o Steve: PSF-creation code within pyphot DOCUMENTED
o Curtis: distutils/setuptools configuration for pyphot DONE
o Leonidas: 
 - reading in the sequence of lcogt images and header information
 - import stellar catalog
  o HIPPARCOS for astrometry
  o SDSS for photometry
 - WCS match of stellar catalog to positions in lcogt image
 - build PSF using pyphot pkfit routines
 - pygfit crowded-field PSF fitting setup for project, using pyphot object detection
 - convert pyphot's PSF description into galfit-format PSF for pygfit use
