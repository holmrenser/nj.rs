PYTHON_LIBDIR := $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYTHON_VERSION := $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LDVERSION'))")
RUSTFLAGS := -L$(PYTHON_LIBDIR) -lpython$(PYTHON_VERSION)

all:
	cp ../README.md README.md
	CARGO_TARGET_DIR=target RUSTFLAGS="${RUSTFLAGS}" cargo build
	CARGO_TARGET_DIR=target env -u CONDA_PREFIX uv run maturin develop --release

test:
	cp ../README.md README.md
	CARGO_TARGET_DIR=target env -u CONDA_PREFIX uv run maturin develop
	uv run pytest

clean:
	CARGO_TARGET_DIR=target cargo clean