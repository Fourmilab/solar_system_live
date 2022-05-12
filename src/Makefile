
VERSION = "\"2.4 (2022-05-12)\""

CFLAGS = -O3 -Wall -DVersion=$(VERSION)
#CFLAGS = -g -Wall -DMALLOC_CHECK_=1
CC	= gcc

OBJS = vplanet.o astro.o gifout.o orbits.o vsop87.o \
	   planetp.o pluto.o asteroid.o astrelem.o cometel.o sglib.o \
	   csvelements.o jplelements.o htmlutil.o uncgi.o strlcpy.o

#	Build test version with public host key
solar:	$(OBJS) packargs_PUB.o
	$(CC) -o solar $(OBJS) packargs_PUB.o -lm
	@rm -f core.*

#	Build production version with Fourmilab private host key
fourmilab:	$(OBJS) packargs_FOURMILAB.o
	$(CC) -o solar $(OBJS) packargs_FOURMILAB.o -lm
	@rm -f core.*

#	Generate BMP planet images from assembled PPM file
images: images.ppm
	ppmtobmp -bpp 8 images.ppm >solar_images.bmp

#	Run regression test against reference output from production application
DDI = | sed '/di=/d'
check: FORCE
	rm -rf regress/cache
	mkdir regress/cache
	#	Test image generation
	echo "1P/Halley 49400 0.58597811 0.96714291 162.26269 111.33249 58.42008 19860205.89532 JPL J863/77" | \
	./solar -oregress/cache/solar1.gif -b0 -t -i640 -e22476948.34564 -dA90 -DO0 47 7
	cmp regress/cache/solar1.gif regress/expected/solar1.gif
	echo "20461 Dioretsa 51388 24.4459301 0.9027402 160.19391 102.19213 296.97833 359.0956660 13.80 0.15 JPL" | \
	./solar -oregress/cache/solar2.gif -b1 -t -i640 -e22476948.34564 -dA45 -DO30 47 7
	cmp regress/cache/solar2.gif regress/expected/solar2.gif
	./solar -oregress/cache/solar3.gif -si -b0 -i320 -v0 -e22476948.34564 -dA33 -DO0 47 7
	cmp regress/cache/solar3.gif regress/expected/solar3.gif
	./solar -oregress/cache/solar4.gif -b2 -i480 -e22476948.34564 -dA90 -DO0 47 7
	cmp regress/cache/solar4.gif regress/expected/solar4.gif
	#	Test HTML output generation
	echo "1P/Halley 49400 0.58597811 0.96714291 162.26269 111.33249 58.42008 19860205.89532 JPL J863/77" | \
	./solar -f -h -b0 -t -i640 -e22476948.34564 -dA90 -DO0 47 7 $(DDI) >regress/cache/solar1.html
	diff regress/cache/solar1.html regress/expected/solar1.html
	echo "20461 Dioretsa 51388 24.4459301 0.9027402 160.19391 102.19213 296.97833 359.0956660 13.80 0.15 JPL" | \
	./solar -f -h -b1 -t -i640 -e22476948.34564 -dA45 -DO30 47 7 $(DDI) >regress/cache/solar2.html
	diff regress/cache/solar2.html regress/expected/solar2.html
	./solar -f -h -si -b0 -i320 -v0 -e22476948.34564 -dA33 -DO0 47 7 $(DDI) >regress/cache/solar3.html
	diff regress/cache/solar3.html regress/expected/solar3.html
	./solar -f -h -b2 -i480 -e22476948.34564 -dA90 -DO0 47 7 $(DDI) >regress/cache/solar4.html
	diff regress/cache/solar4.html regress/expected/solar4.html
	rm -rf regress/cache
	@echo All regression tests passed.

vplanet.o:  vplanet.c vplanet.h

astro.o:    astro.c vplanet.h

astrelem.o: astrelem.c vplanet.h

cometel.o:  cometel.c vplanet.h

csvelements.o:	csvelements.c vplanet.h

orbits.o:   orbits.c vplanet.h

packargs_FOURMILAB.o: packargs.c vplanet.h
	$(CC) -c -o packargs_FOURMILAB.o packargs.c -DFOURMILAB

packargs_PUB.o: packargs.c vplanet.h
	$(CC) -c -o packargs_PUB.o packargs.c

gifout.o:   gifout.c vplanet.h

jplelements.o: jplelements.c vplanet.h

sglib.o:    sglib.c sglib.h

vsop87.o:   vsop87.c vplanet.h

planetp.o: planetp.c vplanet.h

pluto.o:    pluto.c vplanet.h

asteroid.o: asteroid.c vplanet.h

uncgi.o:    uncgi.c vplanet.h

strlcpy.o:  strlcpy.c

htmlutil.o:  htmlutil.c vplanet.h

clean:
	rm -f *.o solar *.bak core *.out *.gif
	rm -rf regress/cache

FORCE:
