FC = ifort

FFLAGS = -O2

LDR = ifort

.F90.o:
	$(FC) $(FFLAGS) $<

test_cool: test_cool.F90 Cool_init.F90 readCloudyTable.F90 ut_hunt.F90 \
	Cool_data.F90
	$(FC) $(FFLAGS) Cool_data.F90 Cool_init.F90 readCloudyTable.F90 \
	    ut_hunt.F90  test_cool.F90 -c
	$(LDR) $(FFLAGS) test_cool.o Cool_init.o readCloudyTable.o \
	    ut_hunt.o Cool_data.o /usr/lib64/libcfitsio.so.0 -o test_cool 
