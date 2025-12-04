#
# Copyright (c) 2018, Istvan Bondar,
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
SHELL:=/bin/bash
HOME:=${HOME}
MAKE:=make
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
ifeq ($(UNAME_M),x86_64)
    ARCH:=64
else ifeq ($(UNAME_M),arm64)
    ARCH:=64
else
    ARCH:=32
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
    LIBEXT := dylib
else
#
#   Linux (tested on CentOS and Ubuntu):
#
    LAPACK := -L$(TARGETLIB) -llapack -lblas
    LIBEXT := so
endif

CFLAGS = -Wall -g -O0 -m$(ARCH) $(MACOS)
CC = gcc
LIBS = $(LAPACK) -L$(TARGETLIB) -lsciLoc -lslbm -lslbmCshell -lgeotesscpp -lm
INCS = -I../include -I../rstt/SLBM/include -I../rstt/SLBM_C_shell/include
COMPILE = ${CC} ${CFLAGS} $(INCS) -o $@ -c

#-------------------------------------------------------------------------------
# targets
#-------------------------------------------------------------------------------

all: ttimes Philippines Ukraine ARA1 ARA2 VacA VacTA VeryBig Jan Deep cleanobjs

%.o: $($*:.o=.c)
	${COMPILE} $*.c

ttimes: sciLocTTimes.o
	$(CC) ${CFLAGS} $(INCS) -o sciLocTTimes sciLocTTimes.o ${LIBS}

Philippines: sciLocMainPhilippines.o
	$(CC) ${CFLAGS} $(INCS) -o PhilippinesTest sciLocMainPhilippines.o ${LIBS}

Ukraine: sciLocMainUkraine.o
	$(CC) ${CFLAGS} $(INCS) -o UkraineTest sciLocMainUkraine.o ${LIBS}

ARA1: sciLocMainARA1.o
	$(CC) ${CFLAGS} $(INCS) -o ARA1Test sciLocMainARA1.o ${LIBS}

ARA2: sciLocMainARA2.o
	$(CC) ${CFLAGS} $(INCS) -o ARA2Test sciLocMainARA2.o ${LIBS}

VacTA: sciLocMainVacTA.o
	$(CC) ${CFLAGS} $(INCS) -o VacTestTA sciLocMainVacTA.o ${LIBS}

VacA: sciLocMainVacA.o
	$(CC) ${CFLAGS} $(INCS) -o VacTestA sciLocMainVacA.o ${LIBS}

VeryBig: sciLocMainVeryBig.o
	$(CC) ${CFLAGS} $(INCS) -o VeryBigTest sciLocMainVeryBig.o ${LIBS}

Jan: sciLocMainJan.o
	$(CC) ${CFLAGS} $(INCS) -o JanTest sciLocMainJan.o ${LIBS}

Deep: sciLocMainDeep.o
	$(CC) ${CFLAGS} $(INCS) -o DeepTest sciLocMainDeep.o ${LIBS}
	@echo "$(blue)---------$(sgr0)"
	@echo "$(blue)done.    $(sgr0)"
	@echo "$(blue)---------$(sgr0)"

cleanobjs:
	-rm -f *.o


.PHONY: all
