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
ifeq ($(UNAME_M),$(filter $(UNAME_M),x86_64 arm64 aarch64))
        export ARCH := 64
else
        export ARCH := 32
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

all: ttimes cleanobjs

%.o: $($*:.o=.c)
	${COMPILE} $*.c

ttimes: sciLocTTimes.o
	$(CC) ${CFLAGS} $(INCS) -o sciLocTTimes sciLocTTimes.o ${LIBS}


cleanobjs:
	-rm -f *.o


.PHONY: all
