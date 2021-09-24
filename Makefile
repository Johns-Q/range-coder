#
#	@file Makefile	@brief	make file for range-coder.
#
#	Copyright (c) 2015 by Johns.  All Rights Reserved.
#
#	Contributor(s):
#
#	License: AGPLv3
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU Affero General Public License as
#	published by the Free Software Foundation, either version 3 of the
#	License.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU Affero General Public License for more details.
#
#	$Id$
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#	Config

#	enable debug
CONFIG += -DDEBUG

#	64-bit version
#CONFIG += -DRANGE_CODER_64

#	32-bit version
CONFIG += -DRANGE_CODER_32

DEFS += $(CONFIG)

#----------------------------------------------------------------------------

VERSION	=	"1.00"
GIT_REV =	$(shell git describe --always 2>/dev/null)

CC=	gcc

#MARCH=	-march=armv6j -mtune=arm1136jf-s -mfpu=vfp -mfloat-abi=softfp
MARCH=	-march=native
#MARCH=	-muclibc
#
OPTIM=	-U_FORTIFY_SOURCE -D__OPTIMIZE__ -O3 -fomit-frame-pointer
CFLAGS= $(MARCH) $(OPTIM) -W -Wall -W -g -pipe \
	-I. $(DEFS) -DVERSION='$(VERSION)' \
	$(if $(GIT_REV), -DGIT_REV='"$(GIT_REV)"')
LIBS	=

OBJS	= range-coder.o
SRCS	= $(OBJS:.o=.c)
HDRS	= range-coder.h

FILES=	Makefile contrib/range-coder.doxyfile

all:	range-coder # range-coder32

$(OBJS):Makefile $(HDRS)

range-coder:	$(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

range-coder32:	range-coder.c
	$(CC) -m32 $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

#----------------------------------------------------------------------------
#	Developer tools

.PHONY: doc indent clean clobber dist

doc:	$(SRCS) $(HDRS) contrib/range-coder.doxyfile
	(cat contrib/range-coder.doxyfile; \
	echo 'PROJECT_NUMBER=${VERSION} $(if $(GIT_REV), (GIT-$(GIT_REV)))'; \
	echo 'INPUT=$(SRCS) $(HDRS) README.md'; \
	) | doxygen -

indent:
	for i in $(OBJS:.o=.c) $(HDRS); do \
		indent $$i; unexpand -a $$i > $$i.up; mv $$i.up $$i; \
	done

clean:
	-rm *.o *~

clobber:	clean
	-rm -rf doc/html range-coder range-coder32

dist:
	tar cjCf .. range-coder-`date +%F-%H`.tar.bz2 \
		$(addprefix range-coder/, $(FILES) $(HDRS) $(OBJS:.o=.c))

dist-git:
	git tag v$(VERSION)
	git archive --format=tar --prefix="range-coder-$(VERSION)/" v$(VERSION) | \
		gzip > range-coder-$(VERSION).tar.gz

install: all
	mkdir -p $(DESTDIR)/usr/local/bin
	install -s range-coder $(DESTDIR)/usr/local/bin/
