FROM gcc:15.2.0

# Install software dependencies for your support libraries
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    unzip \
    git \
    cmake \
    python3 \
    python3-pip \
    python3-dev \
    nano \
    libgmp-dev \
    libmpfr-dev \
    libmpc-dev \
    libisl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install code-server
RUN curl -fsSL https://code-server.dev | sh

# Configure environment to leverage the official extensions store (for Copilot Business)
ENV EXTENSIONS_GALLERY='{"serviceUrl":"https://visualstudio.com","itemUrl":"https://visualstudio.com"}'

WORKDIR /workspace
EXPOSE 8080

CMD ["code-server", "--bind-addr", "0.0.0.0:8080", "/workspace"]
