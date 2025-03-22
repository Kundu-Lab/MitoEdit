#!/bin/bash

# Check if running as root
if [ "$EUID" -ne 0 ]; then
  echo "Please run as root (use sudo)"
  exit 1
fi

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
  echo "Docker is not installed. Please install Docker first."
  exit 1
fi

# Stop the service if it's running
systemctl stop mitoedit.service || true

# Wait a moment for the container to stop
sleep 2

# Build the Docker image
echo "Building Docker image..."
docker build -t mitoedit .

# Restart the service
echo "Starting service..."
systemctl start mitoedit.service

# Check status
systemctl status mitoedit.service --no-pager

echo "MitoEdit has been rebuilt and restarted."
echo "You can access it at https://localhost:443"
echo "To view logs: sudo journalctl -u mitoedit.service -f"