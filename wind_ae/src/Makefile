# Makefile for relaxed_ae, a program to simulate atmospheric escape
# to list defaults, type make -p

BIN=${EXEDIR}relaxed_ae
EXEDIR=../bin/

# use gcc
CC=gcc
CFLAGS=-O3 -Wall
LDFLAGS=-lm
LDFLAGS=-lm -fcommon

NRDIR=../ext/nrecipes
NRSRCS=${NRDIR}/bksub.c \
       ${NRDIR}/bsstep.c \
       ${NRDIR}/rkck.c \
       ${NRDIR}/rkqs.c \
       ${NRDIR}/rk4.c \
       ${NRDIR}/rkdumb.c \
       ${NRDIR}/mmid.c \
       ${NRDIR}/nr_util.c \
       ${NRDIR}/odeint.c \
       ${NRDIR}/pinvs.c \
       ${NRDIR}/polint.c \
       ${NRDIR}/pzextr.c \
       ${NRDIR}/red.c \
       ${NRDIR}/solvde.c

SRCS=difeq.c \
     glq_rates.c \
     intode.c \
     io.c \
     main.c \
     relax.c \
     soe.c \
     utils.c

compile: ${SRCS}
	${CC} ${CFLAGS} -I ${NRDIR} ${SRCS} ${NRSRCS} ${LDFLAGS} -o ${BIN}
