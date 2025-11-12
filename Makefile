VERSION=$(toml get Cargo.toml workspace.package.version)

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

bump:
	cargo workspaces version patch --yes
	@VERSION=$$(toml get Cargo.toml workspace.package.version -r)
	@echo "New version: $$VERSION"
	@make sync-version

sync-version:
	@toml set python/pyproject.toml project.version $(VERSION)
	@jq --arg v "$VERSION" '.version = $v' package.json | sponge package.json

clean:
	$(MAKE) -C nj_lib clean
	$(MAKE) -C wasm clean
	$(MAKE) -C cli clean
	$(MAKE) -C python clean
