#
#  makefile
#  fastp_lite
#
#  Created by Jordan Berg on 3/12/20.
#  Copyright © 2020 Jordan Berg. All rights reserved.
#

CXX = c++

SRCS =  src/main.cpp src/fastqreader.cpp src/read.cpp src/sequence.cpp
OBJS = ${SRCS:.cpp=.o}

MAIN = fastp_lite 

all: ${MAIN} clean
	@echo fastp_lite compiled successfully

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	rm src/*.o