# ── Build stage ───────────────────────────────────────────────────────────────
# Full Rust toolchain image used only for compilation; kept out of the final
# image to minimise attack surface and layer size.
FROM rust:1.77-slim AS builder
WORKDIR /build

# Install C compiler required by mimalloc and any crate with a build.rs that
# calls a C compiler (e.g. ring, openssl).
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc libc6-dev pkg-config \
    && rm -rf /var/lib/apt/lists/*

COPY . .

# Target x86-64-v3 (AVX2) — safe for most cloud instances (AWS c5/m5/r5,
# GCP n2, Azure Dsv3).  Override at build time with:
#   docker build --build-arg RUSTFLAGS="-C target-cpu=native" .
ARG RUSTFLAGS="-C target-cpu=x86-64-v3"
ENV RUSTFLAGS=${RUSTFLAGS}

# Build only the CLI binary; exclude pybioomics (cdylib requires Python headers).
RUN cargo build --release --bin bioomics --workspace --exclude pybioomics

# ── Runtime stage ─────────────────────────────────────────────────────────────
# Minimal Debian image — only ca-certificates is needed for HTTPS access in
# potential future network features.
FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="BioMultiOmics" \
      org.opencontainers.image.description="Fast multi-omics analysis: genomics, transcriptomics, epigenomics" \
      org.opencontainers.image.url="https://github.com/diladeniz/multiomics" \
      org.opencontainers.image.source="https://github.com/diladeniz/multiomics" \
      org.opencontainers.image.licenses="Apache-2.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/bioomics /usr/local/bin/bioomics

# Smoke-test: verify the binary executes and prints a version string.
RUN bioomics --version

ENTRYPOINT ["bioomics"]
CMD ["--help"]
