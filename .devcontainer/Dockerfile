FROM debian:bookworm-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    coreutils \
    curl \
    wget \
    nano \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory (change this to your preferred default directory)
WORKDIR /workspaces/bioinformatics_course

# Change the prompt (PS1) permanently for all bash sessions
RUN echo 'export PS1="(bioinfo) \u@\h:\w$ "' >> /etc/bash.bashrc

# Start with bash
CMD ["/bin/bash"]
