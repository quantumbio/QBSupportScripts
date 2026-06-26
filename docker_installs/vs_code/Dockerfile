FROM gcc:15.2.0

# Install software dependencies for your support libraries
RUN echo 'APT::Update::Post-Invoke {"true";};' > /etc/apt/apt.conf.get && \
    rm -f /etc/apt/apt.conf.d/docker-clean && \
    apt-get update && apt-get install -y \
    curl \
    wget \
    unzip \
    git \
    python3 \
    python3-pip \
    python3-dev \
    libgmp-dev \
    libmpfr-dev \
    libmpc-dev \
    libisl-dev \
    ripgrep \
    nano \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/gcc && ln -sf /usr/local /opt/gcc/current

RUN useradd -m -s /bin/bash developer && echo "developer ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Inject environment variables globally for code-server terminals
ENV PATH="/opt/gcc/current/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/gcc/current/lib64:${LD_LIBRARY_PATH}"
ENV CC="/opt/gcc/current/bin/gcc"
ENV CXX="/opt/gcc/current/bin/g++"

# Install code-server
RUN curl -fsSL https://code-server.dev | sh


# Configure environment to leverage the official extensions store (for Copilot Business)
ENV EXTENSIONS_GALLERY='{"serviceUrl":"https://visualstudio.com","itemUrl":"https://visualstudio.com"}'

USER developer
RUN git config --global pull.rebase false

WORKDIR /workspace
EXPOSE 8080

CMD ["code-server", "--bind-addr", "0.0.0.0:8080", "/workspace"]
