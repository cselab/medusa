# Makefile created by mkmf $Id: Makefile,v 1.9 2006/10/17 14:17:40 pchatela Exp $ 
.SUFFIXES:

.DEFAULT:
	-touch $@
all: ppm_tv
wvic_run.o: ./wvic_run.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_run.f
wvic_iswall.o: ./wvic_iswall.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_iswall.f
wvic_readparams.o: ./wvic_readparams.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_readparams.f
wvic_dumpfield_plot3d.o: ./wvic_dumpfield_plot3d.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_dumpfield_plot3d.f
UpperCase.o: ./UpperCase.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./UpperCase.f
wvic_write_restart.o: ./wvic_write_restart.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_write_restart.f
wvic_suck_restart.o: ./wvic_suck_restart.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_suck_restart.f
wvic_ghost.o: ./wvic_ghost.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_ghost.f
wvic_dump_restart.o: ./wvic_dump_restart.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_dump_restart.f
wvic_tvdrk3.o: ./wvic_tvdrk3.f module_wvic.o wvic_module_io.o ppm_module_fft.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_tvdrk3.f
wvic_died.o: ./wvic_died.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_died.f
wvic_rhs_vort.o: ./wvic_rhs_vort.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_rhs_vort.f
wvic_compvel.o: ./wvic_compvel.f module_wvic.o wvic_les_tdmclipped.o wvic_les_hypervisc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_compvel.f
ppm_module_fft.o: ./ppm_module_fft.f
	cpp -P $(CPPDEFS) $(CPPFLAGS) $< > __$*.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./ppm_module_fft.f
wvic_init_physics_0.o: ./wvic_init_physics_0.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_physics_0.f
wvic.o: ./wvic.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic.f
wvic_map_particles.o: ./wvic_map_particles.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_map_particles.f
module_wvic.o: ./module_wvic.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./module_wvic.f
wvic_check_abort.o: ./wvic_check_abort.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_check_abort.f
wvic_enforcewall.o: ./wvic_enforcewall.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_enforcewall.f
wvic_enforcetv0.o: ./wvic_enforcetv0.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_enforcetv0.f
wvic_enforcetv1.o: ./wvic_enforcetv1.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_enforcetv1.f
wvic_enforcetv2.o: ./wvic_enforcetv2.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_enforcetv2.f
wvic_poisson.o: ./wvic_poisson.f module_wvic.o ppm_module_fft.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_poisson.f
wvic_field2netcdf.o: ./wvic_field2netcdf.f module_wvic.o
	$(FC) $(FFLAGS) -c	./wvic_field2netcdf.f
wvic_netcdf2field.o: ./wvic_netcdf2field.f module_wvic.o
	$(FC) $(FFLAGS) -c	./wvic_netcdf2field.f
wvic_module_io.o: ./wvic_module_io.f module_wvic.o
	cpp -P $(CPPDEFS) $(CPPFLAGS) $< > __$*.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_module_io.f
wvic_init_physics_4.o: ./wvic_init_physics_4.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_physics_4.f
wvic_init_physics_5.o: ./wvic_init_physics_5.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_physics_5.f
wvic_init_physics_6.o: ./wvic_init_physics_6.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       ./wvic_init_physics_6.f
wvic_init_tvphysics_0.o: ./wvic_init_tvphysics_0.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_tvphysics_0.f
wvic_init_tvphysics_1.o: ./wvic_init_tvphysics_1.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_tvphysics_1.f
wvic_init_tvphysics_2.o: ./wvic_init_tvphysics_2.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_tvphysics_2.f
sisec.o: ./sisec.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./sisec.f
rootfinders.o: ./rootfinders.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./rootfinders.f
stopwatch.o: ./stopwatch.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./stopwatch.f
wvic_init_helphysics_0.o: ./wvic_init_helphysics_0.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_helphysics_0.f
wvic_init_stdobjs.o: ./wvic_init_stdobjs.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_stdobjs.f
wvic_stepfunction.o: ./wvic_stepfunction.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_stepfunction.f
wvic_medusa.o: ./wvic_medusa.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_medusa.f
wvic_barycentric.o: ./wvic_barycentric.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_barycentric.f
wvic_stlread.o: ./wvic_stlread.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_stlread.f
wvic_init_cart.o: ./wvic_init_cart.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_init_cart.f
wvic_dgammadt.o: ./wvic_dgammadt.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_dgammadt.f
wvic_penal.o: ./wvic_penal.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_penal.f
wvic_movement.o: ./wvic_movement.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_movement.f
wvic_forcenoca.o: ./wvic_forcenoca.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_forcenoca.f
wvic_les_tdmclipped.o: ./wvic_les_tdmclipped.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_les_tdmclipped.f
wvic_les_hypervisc.o: ./wvic_les_hypervisc.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_les_hypervisc.f
wvic_alloc_field.o: ./wvic_alloc_field.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_alloc_field.f
wvic_diagnostics.o: ./wvic_diagnostics.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_diagnostics.f
wvic_rhs_loca.o: ./wvic_rhs_loca.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_rhs_loca.f
wvic_pbc.o: ./wvic_pbc.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_pbc.f
wvic_readstl.o: ./wvic_readstl.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_readstl.f
wvic_dgammadt_wall.o: ./wvic_dgammadt_wall.f module_wvic.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	./wvic_dgammadt_wall.f
SRC = ./wvic_enforcetv0.f ./wvic_enforcetv1.f ./wvic_enforcetv2.f ./wvic_run.f ./wvic_iswall.f ./wvic_readparams.f ./wvic_dumpfield_plot3d.f ./UpperCase.f ./wvic_write_restart.f ./wvic_suck_restart.f ./wvic_ghost.f ./wvic_dump_restart.f ./wvic_tvdrk3.f ./wvic_died.f ./wvic_rhs_vort.f ./wvic_compvel.f ./ppm_module_fft.f ./wvic_init_physics_0.f ./wvic_init_cart.f ./wvic.f ./wvic_map_particles.f ./module_wvic.f ./wvic_check_abort.f ./wvic_enforcewall.f ./wvic_poisson.f ./wvic_field2netcdf.f ./wvic_netcdf2field.f ./wvic_init_physics_4.f ./wvic_init_physics_5.f ./wvic_init_physics_6.f ./wvic_init_tvphysics_0.f ./wvic_init_tvphysics_1.f ./wvic_init_tvphysics_2.f ./wvic_init_helphysics_0.f ./wvic_init_stdobjs.f ./wvic_stlread.f ./wvic_medusa.f ./wvic_barycentric.f ./wvic_stepfunction.f ./sisec.f ./rootfinders.f ./stopwatch.f /wvic_dgammadt.f /wvic_penal.f /wvic_movement.f ./wvic_forcenoca.f ./wvic_les_tdmclipped.f ./wvic_les_hypervisc.f ./wvic_alloc_field.f ./wvic_diagnostics.f ./wvic_rhs_loca.f ./wvic_pbc.f ./wvic_dgammadt_wall.f ./wvic_module_io.f
OBJ = wvic_enforcetv0.o wvic_enforcetv1.o wvic_enforcetv2.o wvic_run.o wvic_iswall.o wvic_readparams.o wvic_dumpfield_plot3d.o UpperCase.o wvic_write_restart.o wvic_suck_restart.o wvic_ghost.o wvic_dump_restart.o wvic_tvdrk3.o wvic_died.o wvic_rhs_vort.o wvic_compvel.o ppm_module_fft.o wvic_init_physics_0.o wvic.o wvic_map_particles.o module_wvic.o wvic_init_cart.o wvic_check_abort.o wvic_enforcewall.o wvic_poisson.o wvic_field2netcdf.o wvic_init_physics_4.o wvic_dgammadt.o wvic_penal.o wvic_movement.o wvic_forcenoca.o wvic_les_tdmclipped.o wvic_les_hypervisc.o wvic_alloc_field.o wvic_netcdf2field.o wvic_diagnostics.o  wvic_rhs_loca.o wvic_pbc.o wvic_dgammadt_wall.o wvic_init_physics_5.o wvic_init_physics_6.o wvic_init_tvphysics_0.o wvic_init_tvphysics_1.o wvic_init_tvphysics_2.o wvic_init_helphysics_0.o wvic_init_stdobjs.o wvic_stepfunction.o  ./wvic_readstl.o ./wvic_medusa.o ./wvic_barycentric.o ./sisec.o ./rootfinders.o ./stopwatch.o wvic_module_io.o

clean: neat
	-rm -f .cppdefs $(OBJ) ppm_tv
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
ppm_tv: $(OBJ)
	$(LD) $(OBJ) -o ppm_tv $(LDFLAGS)
