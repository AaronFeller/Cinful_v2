#!/bin/bash

# Directory containing the .gz files
INPUT_DIR="dna_data"

# Check if the directory exists
if [ ! -d "$INPUT_DIR" ]; then
  echo "Directory $INPUT_DIR does not exist."
  exit 1
fi

# Function to unpack a .gz file if the decompressed file does not already exist
unpack_file() {
  local gz_file="$1"
  local uncompressed_file="${gz_file%.gz}"

  if [ -f "$uncompressed_file" ]; then
    echo "File $uncompressed_file already exists, skipping."
  else
    echo "Unpacking $gz_file..."
    gunzip "$gz_file"
  fi
}

export -f unpack_file

# Find all .gz files in the directory and unpack them if not already unpacked
# Run up to 10 processes in parallel using xargs
find "$INPUT_DIR" -type f -name '*.gz' | xargs -n 1 -P 40 -I {} bash -c 'unpack_file "$@"' _ {}

echo "Unpacking complete."
