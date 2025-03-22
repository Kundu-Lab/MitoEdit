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

# Check if the Docker image exists, build it if not
if ! docker images | grep -q "mitoedit"; then
  echo "Building Docker image..."
  docker build -t mitoedit .
fi

# Copy service file to systemd directory
cp mitoedit.service /etc/systemd/system/

# Reload systemd to recognize the new service
systemctl daemon-reload

# Stop the service if it's already running
systemctl stop mitoedit.service || true

# Enable the service to start on boot
systemctl enable mitoedit.service

# Start the service
systemctl start mitoedit.service

# Check status
systemctl status mitoedit.service --no-pager

echo "MitoEdit service has been installed and started."
echo "You can manage it with the following commands:"
echo "  sudo systemctl start mitoedit.service"
echo "  sudo systemctl stop mitoedit.service"
echo "  sudo systemctl restart mitoedit.service"
echo "  sudo systemctl status mitoedit.service"