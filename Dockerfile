# Enterprise DeltaChem Environment
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libsm6 \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements and install
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY src/ src/
COPY main.py .
COPY run_production.py .

# Environment Variables
ENV PYTHONPATH="/app"
ENV DATA_DIR="/app/data"
ENV ORCA_PATH="/usr/local/bin/orca"
ENV XTB_PATH="/usr/local/bin/xtb"

# Default command
CMD ["python", "run_production.py"]
