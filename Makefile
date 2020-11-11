include conf.mk

all:
	(cd fftw && \
	 FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" CFLAGS="$(NETCDF_CFLAGS)" ./configure --enable-fortran --enable-single --disable-dependency-tracking --disable-shared --disable-doc --prefix="$(PREFIX)" && \
	     make install) && \
	(cd netcdf && \
	 FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" CFLAGS="$(NETCDF_CFLAGS)" ./configure --disable-f03 --disable-dependency-tracking --disable-shared --disable-doxygen --prefix="$(PREFIX)" && \
	     make install) && \
	(cd ppm && make install) && \
	(cd medusa && make)
