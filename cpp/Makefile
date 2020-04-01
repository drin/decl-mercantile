FLATBUF_HEADERS=src/flatbuffers
GEN_HEADERS=src/generated
FB_META_SOURCES=src/dataformats/fb_meta.cc

INCLUDES=-I${FLATBUF_HEADERS} -I${GEN_HEADERS} -I${HOME}/toolbox/arrow/include -I/usr/include
LD_LIBRARY_PATH=-L${HOME}/toolbox/arrow/lib
LIBS=-lboost_program_options -larrow

CPPFLAGS=${INCLUDES} ${LD_LIBRARY_PATH} ${LIBS}
CC=g++

FB_META_BIN=build/bin/fb_meta

${FB_META_BIN}: ${FB_META_SOURCES}
	mkdir -p build/bin/
	${CC} ${CPPFLAGS} $^ -o $@

clean:
	rm -rf build
