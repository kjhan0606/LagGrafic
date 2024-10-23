#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  meshcorr.mod1.exe

#--- Objects to build

#OBJS		=  readparticle.mod2.o meshcorr.o isctgt.o
OBJS		=  meshcorr.mod1.o 

OS = LINUX_GCC



#--- Targets

SUBDIRS =


all: $(OBJS)
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) all);\
	done 
	rm -f $(PROGS)
	$(MAKE) $(PROGS)

new: 
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) new);\
	done 
	rm -f $(PROGS)
	rm -f $(OBJS)
	$(MAKE) $(PROGS)

clean: 
	rm -f $(OBJS)
	rm -f $(PROGS)
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) clean);\
	done

default: all

$(PROGS): $(OBJS)

include Rules.make
