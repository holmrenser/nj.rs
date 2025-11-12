# top-level Makefile
.PHONY: all python wasm cli lib clean

all: lib wasm cli python

lib:
	$(MAKE) -C nj_lib

wasm:
	$(MAKE) -C wasm

cli:
	$(MAKE) -C cli

python:
	$(MAKE) -C python

clean:
	$(MAKE) -C nj_lib clean
	$(MAKE) -C wasm clean
	$(MAKE) -C cli clean
	$(MAKE) -C python clean
