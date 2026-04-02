VERSION:=$(toml get Cargo.toml workspace.package.version)

.PHONY: all python wasm nj clean

all: nj wasm python

nj:
	$(MAKE) -C nj

wasm:
	$(MAKE) -C wasm

python:
	$(MAKE) -C python

test:
	$(MAKE) -C nj test
	$(MAKE) -C wasm test
	$(MAKE) -C python test

bump-%:
	./release.sh $*

clean:
	$(MAKE) -C nj clean
	$(MAKE) -C wasm clean
	$(MAKE) -C python clean
