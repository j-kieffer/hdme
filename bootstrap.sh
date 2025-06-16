#!/bin/sh

# This script generates the configure script from configure.ac using autoreconf

echo "Generating configure script using autoreconf..."

# Check for autoreconf (part of autoconf package)
if ! command -v autoreconf >/dev/null 2>&1; then
    echo "Error: autoreconf is not installed"
    echo "On macOS: brew install autoconf"
    echo "On Ubuntu/Debian: sudo apt-get install autoconf"
    exit 1
fi

# Check for libtoolize (part of libtool package)
if ! command -v libtoolize >/dev/null 2>&1 && ! command -v glibtoolize >/dev/null 2>&1; then
    echo "Error: libtoolize is not installed"
    echo "On macOS: brew install libtool"
    echo "On Ubuntu/Debian: sudo apt-get install libtool"
    exit 1
fi

# Run autoreconf to generate configure script and all necessary files
autoreconf -fiv

if [ $? -eq 0 ]; then
    echo "Configure script generated successfully!"
    echo "You can now run ./configure to configure the build"
else
    echo "Error: Failed to generate configure script"
    exit 1
fi 