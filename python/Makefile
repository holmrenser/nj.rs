RUSTFLAGS="-L/opt/homebrew/opt/python@3.14/Frameworks/Python.framework/Versions/3.14/lib -lpython3.14"

all:
	RUSTFLAGS=${RUSTFLAGS} cargo build
	maturin develop --release --uv

clean:
	cargo clean