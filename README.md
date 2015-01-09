# lcolens
Time domain gravitational lens photometry of Las Cumbres Observatory observations

Order of things on AAS hack day, 1/8/15
- Steve: PSF-creation code within pyphot DOCUMENTED
- Curtis: distutils/setuptools configuration for pyphot DONE
- Leonidas: 
 - reading in the sequence of lcogt images and header information
 - import stellar catalog
  - result = Vizier.query_region("HE0435-1223", radius=Angle(0.1, "deg"), catalog='NOMAD1')
 - WCS match of stellar catalog to positions in lcogt image
  - Give list of ra/decs and the WCS, and will return x,y
  - http://astropy.readthedocs.org/en/latest/api/astropy.wcs.WCS.html#astropy.wcs.WCS.all_world2pix
 - build PSF using pyphot pkfit routines
 - pygfit crowded-field PSF fitting setup for project, using pyphot object detection
 - convert pyphot's PSF description into galfit-format PSF for pygfit use

