FROM ubuntu:24.04

# ANTLR 4 tool (Java) 4.10.1 + C++ runtime 4.10 from apt
RUN apt-get update && apt-get install -y --no-install-recommends \
    libantlr4-runtime-dev \
    g++ \
    make \
    wget \
    default-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Download ANTLR 4.10.1 complete jar (matches libantlr4-runtime-dev 4.10)
RUN wget -q https://www.antlr.org/download/antlr-4.10.1-complete.jar \
    -O /usr/local/lib/antlr-4.10.1-complete.jar

# Wrapper script
RUN printf '#!/bin/bash\nexec java -jar /usr/local/lib/antlr-4.10.1-complete.jar "$@"\n' \
    > /usr/local/bin/antlr4 && chmod +x /usr/local/bin/antlr4

WORKDIR /workspace
