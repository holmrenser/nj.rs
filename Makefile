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
	cargo test
	$(MAKE) -C wasm test
	$(MAKE) -C python test

bump-%:
	@if [ "$*" != "patch" ] && [ "$*" != "minor" ] && [ "$*" != "major" ]; then \
		echo "Usage: make bump-{patch|minor|major}"; \
		exit 1; \
	fi

	cargo workspaces version $* --yes --no-git-tag
	$(eval VERSION := $(shell toml get Cargo.toml workspace.package.version | tr -d '"'))
	toml set python/pyproject.toml project.version "$$VERSION" | sponge python/pyproject.toml
	jq --arg v "$$VERSION" '.version = $$v' wasm/package.json | sponge wasm/package.json
	git commit -am "Release v$$VERSION"
	git tag v$$VERSION

sync-version:
	echo "syncing version to $(VERSION)"
	toml set python/pyproject.toml project.version "$(VERSION)" | sponge python/pyproject.toml
	jq --arg v $(VERSION) '.version = $$v' wasm/package.json | sponge wasm/package.json
	git tag -a "v$(VERSION)" -m "Version $(VERSION)"

clean:
	$(MAKE) -C nj clean
	$(MAKE) -C wasm clean
	$(MAKE) -C python clean
