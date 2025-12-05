VERSION:=$(toml get Cargo.toml workspace.package.version)

.PHONY: all python wasm nj clean

all: nj wasm python

nj:
	$(MAKE) -C nj

wasm:
	$(MAKE) -C wasm

python:
	$(MAKE) -C python

bump:
	cargo workspaces version patch --yes --no-git-commit
	$(MAKE) sync-version VERSION=$$(toml get Cargo.toml workspace.package.version)

sync-version:
	echo "syncing version to $(VERSION)"
	toml set python/pyproject.toml project.version "$(VERSION)" | sponge python/pyproject.toml
	jq --arg v $(VERSION) '.version = $$v' wasm/package.json | sponge wasm/package.json
	git tag -a "v$(VERSION)" -m "Version $(VERSION)"

clean:
	$(MAKE) -C nj clean
	$(MAKE) -C wasm clean
	$(MAKE) -C python clean
