all: lib cli

lib:
	CARGO_TARGET_DIR=target cargo build --release

cli:
	CARGO_TARGET_DIR=target cargo build --release --features cli --bin nj

test:
	CARGO_TARGET_DIR=target cargo test
	CARGO_TARGET_DIR=target cargo test --features cli

clean:
	CARGO_TARGET_DIR=target cargo clean