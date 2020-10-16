CC = g++
DEBUG = -g 
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
MALLOC_CHECK_ = 2
#CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11"
#QMAKE_CXXFLAGS += -std=c++11	
sc.o: sc.h global.h sc.cpp 
	$(CC) $(CFLAGS) sc.cpp
	

vars.o: vars.h global.h vars.cpp 
	$(CC) $(CFLAGS) vars.cpp

###############--NEXUS---
# run make -f makefile nex


nexmul.o: sc.h global.h vars.h nexmul.cpp 
	$(CC) $(CFLAGS) nexmul.cpp	
		
nex : nexmul.o sc.o vars.o
	$(CC) $(LFLAGS) nexmul.o vars.o sc.o -o nex

###############--NEXUS---



###############--mcklnsmooth---
# run make -f makefile nex


mcklnsmooth.o: sc.h global.h vars.h mcklnsmooth.cpp 
	$(CC) $(CFLAGS) mcklnsmooth.cpp	
		
smooth : mcklnsmooth.o sc.o vars.o
	$(CC) $(LFLAGS) mcklnsmooth.o vars.o sc.o -o smooth

###############--mcklnsmooth---

###############--trento---
# run make -f makefile tnt


trento.o: sc.h global.h vars.h trento.cpp 
	$(CC) $(CFLAGS) trento.cpp	
		
tnt : trento.o sc.o vars.o
	$(CC) $(LFLAGS) trento.o vars.o sc.o -o tnt

###############--trento---

###############--tll---
# run make -f makefile tall


tall.o: sc.h global.h vars.h tall.cpp 
	$(CC) $(CFLAGS) tall.cpp	
		
tll : tall.o sc.o vars.o
	$(CC) $(LFLAGS) tall.o vars.o sc.o -o tll

###############--tll---

###############--trento cumulants---
# run make -f makefile trt


trt2.o: sc.h global.h vars.h trt2.cpp 
	$(CC) $(CFLAGS) trt2.cpp	
		
trt : trt2.o sc.o vars.o
	$(CC) $(LFLAGS) trt2.o vars.o sc.o -o trt

###############--trento---

###############--trento ecc---
# run make -f makefile tntecc


tntecc.o: sc.h global.h vars.h tntecc.cpp 
	$(CC) $(CFLAGS) tntecc.cpp	
		
tntecc : tntecc.o sc.o vars.o
	$(CC) $(LFLAGS) tntecc.o vars.o sc.o -o tntecc

###############--trento---

###############--jazma ecc---
# run make -f makefile jazecc


jazecc.o: sc.h global.h vars.h jazecc.cpp 
	$(CC) $(CFLAGS) jazecc.cpp	
		
jazecc : jazecc.o sc.o vars.o
	$(CC) $(LFLAGS) jazecc.o vars.o sc.o -o jazecc

###############--jazma---


###############--cgcrhic---
# run make -f makefile nex


rcgc.o: sc.h global.h vars.h rcgc.cpp 
	$(CC) $(CFLAGS) rcgc.cpp	
		
rhic : rcgc.o sc.o vars.o
	$(CC) $(LFLAGS) rcgc.o vars.o sc.o -o rhic

###############--cgcrhic---

###############--Glasma---
# run make -f makefile glas

glasmul.o: sc.h global.h  vars.h  glasmul.cpp 
	$(CC) $(CFLAGS) glasmul.cpp
			
glas : glasmul.o sc.o vars.o 
	$(CC) $(LFLAGS) glasmul.o vars.o sc.o -o glas

###############--Glasma---

###############--eccentricities---
# run make -f makefile ecc

eccmul.o: sc.h global.h vars.h eccmul.cpp 
	$(CC) $(CFLAGS) eccmul.cpp
			
ecc : eccmul.o sc.o vars.o 
	$(CC) $(LFLAGS) eccmul.o vars.o sc.o -o ecc

###############--eccentricities---

###############--LHC---
# run make -f makefile lhc

lhcmul.o: sc.h global.h lhcmul.cpp 
	$(CC) $(CFLAGS) lhcmul.cpp
			
lhc : lhcmul.o sc.o
	$(CC) $(LFLAGS) lhcmul.o sc.o -o lhc

###############--LHC---

###############--LHC2---
# run make -f makefile lhc

lhcmul2.o: sc.h global.h lhc2mul.cpp 
	$(CC) $(CFLAGS) lhc2mul.cpp
			
lhc2 : lhcmul2.o sc.o
	$(CC) $(LFLAGS) lhc2mul.o sc.o -o lhc2

###############--LHC---

###############--Glauber---
# run make -f makefile glb

glbmul.o: sc.h global.h glbmul.cpp 
	$(CC) $(CFLAGS) glbmul.cpp
			
glb : glbmul.o sc.o
	$(CC) $(LFLAGS) glbmul.o sc.o -o glb

###############--Glauber---


###############--Glauber LHC---
# run make -f makefile glblhcv2

glblhc.o: sc.h global.h vars.h glblhc.cpp 
	$(CC) $(CFLAGS) glblhc.cpp
			
glblhcv2 : glblhc.o vars.o  sc.o
	$(CC) $(LFLAGS) glblhc.o vars.o sc.o -o glblhcv2

###############--Glauber---


