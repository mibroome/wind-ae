# Makefile for relaxed_ae, a program to calculate atmospheric escape

DIRS = bin

all: dirs compile

dirs:
	-@for i in ${DIRS}; do \
	(if [ -d $$i ]; then \
		echo DIR $$i already exists; \
	else \
		echo DIR $$i created with input_files; \
		mkdir $$i; \
		cp -r input_files $$i; \
	fi); done

compile:
	(cd src; ${MAKE} compile)

clean:
	(rm bin/relaxed_ae)
