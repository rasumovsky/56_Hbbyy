# -*- mode: makefile -*-
#
# Makefile containing platform dependencies different projects.
MAKEFLAGS = --no-print-directory -r -s -j2


ARCH_LOC_1 := $(wildcard $(shell root-config --prefix)/test/Makefile.arch)
ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/etc/Makefile.arch)
ARCH_LOC_3 := $(wildcard $(shell root-config --prefix)/share/doc/root/test/Makefile.arch)
ifneq ($(strip $(ARCH_LOC_1)),)
  $(info Using $(ARCH_LOC_1))
  include $(ARCH_LOC_1)
else
  ifneq ($(strip $(ARCH_LOC_2)),)
    $(info Using $(ARCH_LOC_2))
    include $(ARCH_LOC_2)
  else
    ifneq ($(strip $(ARCH_LOC_3)),)
      $(info Using $(ARCH_LOC_3))
      include $(ARCH_LOC_3)
    else
      $(error Could not find Makefile.arch!)
    endif
  endif
endif


CXXFLAGS += -Wall -Wno-overloaded-virtual -Wno-unused

INCLUDES += -I./inc

VPATH	= ./src ./inc ./ws

GLIBS	+= -lTMVA -lMLP
GLIBS	+= -lTreePlayer -lProof -lProofPlayer -lutil -lRooFit -lRooFitCore  -lRooStats -lFoam -lMinuit -lHistFactory -lXMLParser -lXMLIO -lCore -lGpad -lMathCore  -lPhysics
.PHONY:

OBJS_Template		= obj/template.o
DEPS_Template		:= $(OBJS_Template:.o=.d) 

bin/%	: obj/%.o obj/statistics.o obj/statisticsDict.o obj/RooBernsteinM.o obj/RooBernsteinMDict.o obj/HggTwoSidedCBPdf.o obj/HggTwoSidedCBPdfDict.o obj/FlexibleInterpVarMkII.o obj/FlexibleInterpVarMkIIDict.o obj/CommonFunc.o obj/Config.o obj/DHCheckJobs.o obj/DHAnalysis.o obj/DHDataReader.o obj/SigParam.o obj/DHWorkspace.o obj/DHTestStat.o obj/DHToyTree.o obj/DHToyAnalysis.o obj/DHSystematics.o 

	@echo "Linking " $@
	echo $(LD) $(LDFLAGS) $^ $(GLIBS) -o $@	
	@$(LD) $(LDFLAGS) $^ $(GLIBS) -o $@ 

obj/%.o : %.cxx
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)

obj/%.o : %.cc
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)

obj/%.o : %.C
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -O2 -c $< -MD -o $@ $(INCLUDES)


obj/%Dict.o:	%.h
		@echo "Compiling $@"
		echo rootcint -f $*Dict.cxx -c $*.h $*LinkDef.h
		rootcint -f $*Dict.cxx -c inc/$*.h inc/$*LinkDef.h
		echo $(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o $*Dict.o
		$(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o obj/$*Dict.o
		rm -f $*Dict.cxx $*Dict.h 

obj/%Dict.o:	%.hh
		@echo "Compiling $@"
		echo rootcint -f $*Dict.cxx -c $*.hh $*LinkDef.h
		rootcint -f $*Dict.cxx -c inc/$*.hh inc/$*LinkDef.h
		echo $(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o $*Dict.o
		$(CXX) $(INCLUDES) $(CXXFLAGS) -c $*Dict.cxx -o obj/$*Dict.o
		rm -f $*Dict.cxx $*Dict.h 

clean:
	@echo "Cleaning $<..."
	rm -fr *~ obj/*.o */*~ *_Dict.* *.a 
	@echo "Done"
