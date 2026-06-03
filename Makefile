VERSION:=$(toml get Cargo.toml workspace.package.version)

.PHONY: all python wasm nj bench clean

all: lib cli wasm python

lib:
	$(MAKE) -C nj lib

cli:
	$(MAKE) -C nj cli

wasm:
	$(MAKE) -C wasm

python:
	$(MAKE) -C python

test:
	$(MAKE) -C nj test
	$(MAKE) -C wasm test
	$(MAKE) -C python test

bench:
	$(MAKE) -C nj bench

bump-%:
	./release.sh $*

clean:
	$(MAKE) -C nj clean
	$(MAKE) -C wasm clean
	$(MAKE) -C python clean
