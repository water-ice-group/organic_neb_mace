#
#SUFFIX= -c 
SUFFIX= -c -O3
#LIBS = -static 
#SUFFIX= -c -fpe2
FSALGEM= reac.f poten.f ffopt.f shanno.f version.f vibra.f blas.f
OBJECT= reac.o poten.o ffopt.o shanno.o version.o vibra.o blas.o
INCLUD= cbka.blk cbkm.blk opt.blk 
#
        reac:$(UPDATE) $(OBJECT) $(INCLUD) 
	f77 -o $@ $(OBJECT) $(LIBS)
	date >> compile.log
reac.o: reac.f $(INCLUD) 
	f77 $(SUFFIX) reac.f
poten.o: poten.f $(INCLUD) 
	f77 $(SUFFIX) poten.f
ffopt.o: ffopt.f $(INCLUD) 
	f77 $(SUFFIX) ffopt.f
shanno.o: shanno.f $(INCLUD) 
	f77 $(SUFFIX) shanno.f
version.o: version.f $(INCLUD) 
	f77 $(SUFFIX) version.f
vibra.o: vibra.f $(INCLUD)
	f77 $(SUFFIX) vibra.f
blas.o: blas.f $(INCLUD)
	f77 $(SUFFIX) blas.f
