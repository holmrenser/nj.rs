RUSTFLAGS="-L/opt/homebrew/opt/python@3.14/Frameworks/Python.framework/Versions/3.14/lib -lpython3.14"

all:
	CARGO_TARGET_DIR=target RUSTFLAGS=${RUSTFLAGS} cargo build
	CARGO_TARGET_DIR=target maturin develop --release --uv

clean:
	CARGO_TARGET_DIR=target cargo clean