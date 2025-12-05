#
# Copyright (c) 2021-2024, Istvan Bondar,
# Written by Istvan Bondar, ibondar2014@gmail.com
#
# BSD Open Source License.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the <organization> nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
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
# Makefile to create RSTT and sciLoc libraries
#
#
################################################################################
# initializations
################################################################################
SCILOC_VERSION := 4.3
HOME:= ${HOME}
MAKE:=make
export SCILOCROOT = $(CURDIR)
TARGETLIB := $(HOME)/lib
#
# Determine OS and set OS-specific compiler variables
#    if the command below does not work in zsh add the line below to .zshrc
#    unsetopt nomatch
#
OS:=$(shell uname -s | tr [:upper:] [:lower:])
ifeq ($(OS),darwin)
#
#   Mac OS X 10.6 and later:
#    - LAPACK is part of the Accelerate framework
#    - LD_LIBRARY_PATH is disabled; create ~/lib instead
#
	LIBEXT := dylib
	LDOPTS := -dynamiclib
else
#
#   Linux:
#
	LIBEXT := so
	LDOPTS := -shared -z defs
endif
# color definitions
blue := $(shell tput setaf 6)
sgr0 := $(shell tput sgr0)

#-------------------------------------------------------------------------------
# targets
#-------------------------------------------------------------------------------
all: setup rstt iLoc ttimes

# command that tells us which recipes don't build source files
.PHONY: all setup rstt iLoc ttimes

setup:
#
#   create target library if necessary
#
	@echo "$(blue)-----------------------$(sgr0)"
	@echo "$(blue)Building sciLoc v$(ILOC_VERSION)$(sgr0)"
	@echo "$(blue)Creating directories.. $(sgr0)"
	@echo "$(blue)-----------------------$(sgr0)"
	@ [ -d $(TARGETLIB) ] || mkdir $(TARGETLIB)
	@ [ -d $(HOME)/bin ] || mkdir $(HOME)/bin
#
#   compile RSTT libraries: geotess, C++ and C interface
#
rstt:
	@echo "$(blue)--------------------------------------$(sgr0)"
	@echo "$(blue)Creating RSTT C interface libraries.. $(sgr0)"
	@echo "$(blue)--------------------------------------$(sgr0)"
	$(MAKE) -C rstt cleanall
	$(MAKE) -C rstt geotess
	$(MAKE) -C rstt slbm
	$(MAKE) -C rstt slbmc
	$(MAKE) -C rstt clean
	@echo "$(blue)---------$(sgr0)"
	@echo "$(blue)done     $(sgr0)"
	@echo "$(blue)---------$(sgr0)"
#
#   compile sciLoc library and ttimes API
#
iLoc:
	@echo "$(blue)-------------------------$(sgr0)"
	@echo "$(blue)Compiling iLoc library.. $(sgr0)"
	@echo "$(blue)-------------------------$(sgr0)"
	$(MAKE) -C src all
	@echo "$(blue)---------$(sgr0)"
	@echo "$(blue)done     $(sgr0)"
	@echo "$(blue)---------$(sgr0)"
#
#   compile ttimes API
#
ttimes:
	@echo "$(blue)-----------------------------------$(sgr0)"
	@echo "$(blue)Compiling tests and sciLocTTimes.. $(sgr0)"
	@echo "$(blue)-----------------------------------$(sgr0)"
	$(MAKE) -C TTimes all
	mv TTimes/sciLocTTimes $(HOME)/bin
	@echo "$(blue)---------$(sgr0)"
	@echo "$(blue)done.    $(sgr0)"
	@echo "$(blue)---------$(sgr0)"


