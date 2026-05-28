include conf.mk

all:
	(cd fftw && \
	 ./configure --enable-fortran --enable-single --disable-dependency-tracking --disable-shared --disable-doc --prefix="$(PREFIX)" FC="$(FC)" FCFLAGS="$(FCFLAGS) $(FXFLAGS)" && \
	     make install) && \
	(cd ppm && make install) && \
	(cd medusa && make)
