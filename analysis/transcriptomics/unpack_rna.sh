#!/bin/bash

# Directory containing the .tar files
INPUT_DIR="rna_data"

# Check if the directory exists
if [ ! -d "$INPUT_DIR" ]; then
  echo "Directory $INPUT_DIR does not exist."
  exit 1
fi

# Function to unpack a .tar file if no matching base name files exist (excluding .tar files)
unpack_tar() {
  local tar_file="$1"
  local base_name="${tar_file%.tar}"

  # Check if any file or directory with the base name exists, excluding .tar files
  if ls "${base_name}"* 1> /dev/null 2>&1 | grep -v "\.tar$" >/dev/null; then
    echo "A file or directory starting with $base_name already exists (excluding .tar), skipping."
  else
    echo "Unpacking $tar_file..."
    tar -xf "$tar_file"
  fi
}

export -f unpack_tar

# Find all .tar files in the directory and unpack them if not already unpacked
# Run up to 10 processes in parallel using xargs
find "$INPUT_DIR" -type f -name '*.tar' | xargs -n 1 -P 10 -I {} bash -c 'unpack_tar "$@"' _ {}

echo "Unpacking complete."
