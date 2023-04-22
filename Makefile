include conf.mk

all:
	(cd fftw && \
	 ./configure --enable-fortran --enable-single --disable-dependency-tracking --disable-shared --disable-doc --prefix="$(PREFIX)" FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" CFLAGS="$(NETCDF_CFLAGS)" && \
	     make install) && \
	(cd netcdf && \
	p=$(NETCDF_CFLAGS) && \
	 ./configure --disable-f03 --disable-dependency-tracking --disable-shared --disable-doxygen --prefix="$(PREFIX)" FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" CFLAGS="$$p" && \
	     make install) && \
	(cd ppm && make install) && \
	(cd medusa && make)
