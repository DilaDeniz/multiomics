#!/bin/bash
set -euxo pipefail

export CARGO_HOME="${SRC_DIR}/.cargo-home"

# Override the repo's .cargo/config.toml which targets x86-64-v3 (AVX2).
# conda-forge/bioconda CI spans x86-64-v1, aarch64, and osx-arm64 — let
# rustc pick the best target for each host instead.
export RUSTFLAGS="-C target-cpu=native -C strip=symbols"

cargo install \
    --path cli \
    --root "${PREFIX}" \
    --locked \
    --features cnv,atac
