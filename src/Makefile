#
# Copyright (c) 2021, Istvan Bondar,
# Written by Istvan Bondar, ibondar2014@gmail.com
#
# BSD Open Source License.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
# sciloc library Makefile
#
#
SHELL:=/bin/bash
HOME:=${HOME}
MAKE:=make
ISLA:=1
TARGETLIB := $(HOME)/lib

#
# color definitions
#
blue := $(shell tput setaf 6)
sgr0 := $(shell tput sgr0)
#
#
# Determine OS specifics
#
#
OS:=$(shell uname -s | tr [:upper:] [:lower:])
#
# detect bit size
#
UNAME_M:=$(shell uname -m)
ifeq ($(UNAME_M),$(filter $(UNAME_M),x86_64 arm64 aarch64))
        export ARCH := 64
else
        export ARCH := 32
endif
#
# get linux distribution
#
ifeq ($(OS),linux)
    LSB:=$(shell lsb_release -i | tr [:upper:] [:lower:])
    RL:=$(word 3, $(LSB))
    ifeq ($(RL),ubuntu)
        REL:=Ubuntu
    else ifeq ($(RL),debian)
        REL:=Ubuntu
    else ifeq ($(RL),centos)
        REL:=CentOS
    else ifeq ($(RL),redhat)
        REL:=CentOS
    else ifeq ($(RL),suse)
        REL:=CentOS
    else
        REL:=CentOS
    endif
endif
#
#
# Set OS-specific compiler variables
#
#
ifeq ($(OS),darwin)
#
# Mac OS X 10.6 and later:
#    - LAPACK is part of the Accelerate framework
#    - LD_LIBRARY_PATH is disabled; create ~/lib instead
    LAPACK := -framework Accelerate
    MACOS := -DMACOSX=1
    LDOPTS := -dynamiclib
    LIBEXT := dylib
    CCFLAGS = -Wall -g -O0 -m$(ARCH) $(MACOS)
else
#
#   Linux (tested on CentOS and Ubuntu):
#
    LAPACK := -L$(TARGETLIB) -llapack -lblas
    LDOPTS := -shared -z defs
    LIBEXT := so
    CCFLAGS = -Wall -g -O0 -m$(ARCH) -fPIC
#
#   Ubuntu, Debian specific
#
    ifeq ($(REL),Ubuntu)
#
#       check for the existence of lapack and blas libraries
#
        ifeq (,$(wildcard $(TARGETLIB)/liblapack.$(LIBEXT)))
            ifneq (,$(wildcard /usr/lib/liblapack.$(LIBEXT)))
                LLIB:=/usr/lib/liblapack.$(LIBEXT)
                BLIB:=/usr/lib/libblas.$(LIBEXT)
            else ifneq (,$(wildcard /usr/lib/liblapack.$(LIBEXT).3))
                LLIB:=/usr/lib/liblapack.$(LIBEXT).3
                BLIB:=/usr/lib/libblas.$(LIBEXT).3
            else ifneq (,$(wildcard /usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT)))
                LLIB:=/usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT)
                BLIB:=/usr/lib/x86_64-linux-gnu/libblas.$(LIBEXT)
            else ifneq (,$(wildcard /usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT).3))
                LLIB:=/usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT).3
                BLIB:=/usr/lib/x86_64-linux-gnu/libblas.$(LIBEXT).3
            else
                ISLA:=0
            endif
        endif
#
#   CentOS, Red Hat, SuSe
#
    else
#
#       check for the existence of lapack and blas libraries
#
        ifeq (,$(wildcard $(TARGETLIB)/liblapack.$(LIBEXT)))
            ifneq (,$(wildcard /usr/lib64/liblapack.$(LIBEXT)))
                LLIB:=/usr/lib64/liblapack.$(LIBEXT)
                BLIB:=/usr/lib64/libblas.$(LIBEXT)
            else ifneq (,$(wildcard /usr/lib64/liblapack.$(LIBEXT).3))
                LLIB:=/usr/lib64/liblapack.$(LIBEXT).3
                BLIB:=/usr/lib64/libblas.$(LIBEXT).3
            else ifneq (,$(wildcard /usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT)))
                LLIB:=/usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT)
                BLIB:=/usr/lib/x86_64-linux-gnu/libblas.$(LIBEXT)
            else ifneq (,$(wildcard /usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT).3))
                LLIB:=/usr/lib/x86_64-linux-gnu/liblapack.$(LIBEXT).3
                BLIB:=/usr/lib/x86_64-linux-gnu/libblas.$(LIBEXT).3
            else
                ISLA:=0
            endif
        endif
    endif
endif
#
# compiler and linker flags
#
CC = gcc
OUTPUT = $(TARGETLIB)/libsciLoc.$(LIBEXT)
SCILOCINC = -I../include -I../rstt/SLBM/include -I../rstt/SLBM_C_shell/include
ILOCLIBS = $(LAPACK) -L$(TARGETLIB) -lslbm -lslbmCshell -lgeotesscpp -lm
COMPILE = ${CC} ${CCFLAGS} $(SCILOCINC) -o $@ -c
LINK = ${CC} ${CCFLAGS} $(LDOPTS) -o $(OUTPUT) ${OBJS} $(ILOCLIBS)

#
#
# Source files
#
#
SRCS =  sciLocCluster.c \
	sciLocDataCovariance.c \
	sciLocDepthPhases.c \
	sciLocDistAzimuth.c \
	sciLocFreeMemory.c \
	sciLocGregion.c \
	sciLocInitializations.c \
	sciLocInterpolate.c \
	sciLocLocalTT.c \
	sciLocLocationQuality.c \
	sciLocLocator.c \
	sciLocNA.c \
	sciLocPhaseIdentification.c \
	sciLocPrintEvent.c \
	sciLocReadAuxDataFiles.c \
	sciLocSVD.c \
	sciLocTravelTimes.c \
	sciLocTravelTimeAPI.c

OBJS = $(SRCS:.c=.o)

################################################################################
# recipes
################################################################################


all: checks sciLoc clean
.PHONY: all checks sciLoc clean

checks: clean
	@ echo "$(blue)----------------------------------------$(sgr0)"
	@ echo "$(blue)Checking libraries... $(sgr0)"
	@ echo "$(blue)----------------------------------------$(sgr0)"
#
# exit with error if lapack and blas libraries are missing
#
	@echo "$(blue)liblapack and libblas $(sgr0)"
ifeq ($(ISLA),0)
	$(error(lapack and blas libraries are not found!))
else
#
# create symbolic links of lapack and blas libraries to ~/lib
#
ifeq ($(OS),linux)
ifeq (,$(wildcard $(TARGETLIB)/liblapack.$(LIBEXT)))
	ln -srf $(LLIB) $(TARGETLIB)
	ln -srf $(BLIB) $(TARGETLIB)
endif
endif
endif

%.o: $($*:.o=.c)
	${COMPILE} $*.c

sciLoc: ${OBJS}
	@echo "$(blue)--------------------$(sgr0)"
	@echo "$(blue)building $(OUTPUT)  $(sgr0)"
	@echo "$(blue)--------------------$(sgr0)"
	${LINK}
	rm -f *.o
	@echo "$(blue)---------$(sgr0)"
	@echo "$(blue)done.    $(sgr0)"
	@echo "$(blue)---------$(sgr0)"

clean:
	@echo "$(blue)Cleaning up object files $(sgr0)"
	rm -f *.o
