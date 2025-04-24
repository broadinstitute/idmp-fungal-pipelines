FROM marcoteix/bracken:1.0.0

# Set the working directory
WORKDIR /app

# Download and extract Kraken2 database
RUN mkdir -p db && \
    wget -O kraken2_db.tar.gz "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz" && \
    tar -C ./db/ -xzvf kraken2_db.tar.gz && \
    rm kraken2_db.tar.gz
