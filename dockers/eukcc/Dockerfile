FROM --platform=linux/amd64 us.gcr.io/broad-gotc-prod/eukcc:2.1.3

# Set the working directory
WORKDIR /app

# Download and extract Kraken2 database
RUN mkdir -p db && \
    wget -O eukcc2_db.tar.gz "http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz" && \
    tar -C ./db/ -xzvf eukcc2_db.tar.gz && \
    rm eukcc2_db.tar.gz


ENV EUKCC2_DB=/app/db/eukcc2_db_ver_1.1