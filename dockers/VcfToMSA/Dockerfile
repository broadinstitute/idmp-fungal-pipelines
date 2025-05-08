# Use an official Python base image
FROM python:3.10-slim

LABEL description="Converts VCF to PHYLIP, FASTA, NEXUS, or binary NEXUS formats for phylogenetic analysis"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Copy the script into the container
COPY vcf2matrix.py .

# Make script executable
RUN chmod +x vcf2matrix.py

# Define the default command
ENTRYPOINT ["./vcf2matrix.py"]

