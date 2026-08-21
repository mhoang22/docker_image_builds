FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    openjdk-21-jdk \
    bash \
    python3 \
    python3-pip \
    procps \
    samtools \
    curl \
    wget \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# WORKDIR /app
# Copy the JAR file into the container
# COPY PepQueryMHC.v1.0.8.jar /app/PepQueryMHC.jar
# ENTRYPOINT ["java", "-jar", "/app/PepQueryMHC.jar"]

WORKDIR /opt/pepquerymhc
COPY PepQueryMHC.v1.0.8.jar /opt/pepquerymhc/PepQueryMHC.jar
RUN chmod 644 /opt/pepquerymhc/PepQueryMHC.jar
ENTRYPOINT ["java", "-jar", "/opt/pepquerymhc/PepQueryMHC.jar"]

CMD ["--help"]
