SRCS=util.cc random.cc pri_queue.cc m_srp.cc kl_srp.cc lcsb.cc interface.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o lcsb ${OBJS}

util.o: util.h

random.o: random.h

pri_queue.o: pri_queue.h

m_srp.o: m_srp.h

kl_srp.o: kl_srp.h

lcsb.o: lcsb.h

interface.o: interface.h

main.o:

clean:
	-rm ${OBJS} lcsb
