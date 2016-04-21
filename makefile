objects=module_constants.f90 Model_sub.f90 Grid_sub.f90 Initialization_sub.f90 \
	IBM_module.f90 LDDRK_sub.f90 F_sub.f90 Absorbing_Coeff_sub.f90 \
		Distribution_sub.f90 Periodic_Boundary_sub.f90 \
			DRP7_sub.f90 Get_dUdV_sub.f90 Search_Index_sub.f90 \
				main.f90 lapack_LINUX.a blas_LINUX.a
scattering.exe: $(objects)
	gfortran $(objects) -o scattering.exe
