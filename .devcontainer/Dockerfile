FROM debian:bookworm-slim

# Atualiza pacotes e instala ferramentas essenciais
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    coreutils \
    curl \
    wget \
    vim \
    git \
    && rm -rf /var/lib/apt/lists/*

# Define o shell padrão como bash
CMD ["/bin/bash"]

