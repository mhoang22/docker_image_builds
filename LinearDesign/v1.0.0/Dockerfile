# ─────────────────────────────────────────────────────────────────────────────
# LinearDesign v1.0.0 — mRNA Design Tool
# Image : myhoang04/LinearDesign:v1.0.0
# Repo  : https://github.com/LinearDesignSoftware/LinearDesign
# ─────────────────────────────────────────────────────────────────────────────

# --- Stage 1: Build ──────────────────────────────────────────────────────────
FROM ubuntu:22.04 AS builder

# Avoid interactive prompts during apt installs
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        g++ \
        make \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Clone the official LinearDesign repository (pinned to v1.0.0)
WORKDIR /opt
RUN git clone --depth 1 --branch v1.0.0 \
        https://github.com/LinearDesignSoftware/LinearDesign.git lineardesign

# Compile
WORKDIR /opt/lineardesign
RUN make

# ─────────────────────────────────────────────────────────────────────────────
# --- Stage 2: Runtime ────────────────────────────────────────────────────────
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Minimal runtime dependencies (g++ runtime libs for the compiled binary + .so)
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
        libgomp1 \
        python3 \
    && rm -rf /var/lib/apt/lists/*

# Copy compiled artefacts from builder
WORKDIR /opt/lineardesign
COPY --from=builder /opt/lineardesign /opt/lineardesign

# Make sure the shared library is found at runtime
ENV LD_LIBRARY_PATH=/opt/lineardesign:${LD_LIBRARY_PATH}

# Expose the tool on PATH
ENV PATH=/opt/lineardesign:${PATH}

# Default working directory for user data / codon tables
WORKDIR /data

# Health-check: verify the binary runs
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD echo "M" | /opt/lineardesign/lineardesign > /dev/null 2>&1 || exit 1

# Entrypoint — users can pass any lineardesign flags directly
ENTRYPOINT ["/opt/lineardesign/lineardesign"]
CMD ["--help"]
