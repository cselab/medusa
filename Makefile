include conf.mk

all: 
	(cd netcdf && \
	 FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" ./configure --disable-f03 --disable-dependency-tracking --disable-shared --prefix="$(PREFIX)" && \
	     make install) && \
	(cd ppm && make install) && \
	(cd medusa && make)
