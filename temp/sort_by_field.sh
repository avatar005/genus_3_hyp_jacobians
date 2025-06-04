#!/usr/bin/env bash
set -euo pipefail

# Enable nullglob so that the loop skips if no files match
shopt -s nullglob

for file in data_*.txt; do
  # Strip the .txt extension to work on the basename
  name="${file%.txt}"

  # Temporarily split on underscores
  IFS_OLD="$IFS"
  IFS="_"
  read -ra parts <<< "$name"
  IFS="$IFS_OLD"

  # Need at least two fields (tuple + last) to classify
  if (( ${#parts[@]} < 2 )); then
    echo "Skipping '$file' (unexpected format)"
    continue
  fi

  # Second‐to‐last is the tuple field, last is the final number
  tuple="${parts[$(( ${#parts[@]} - 2 ))]}"
  last="${parts[$(( ${#parts[@]} - 1 ))]}"

  # Remove spaces from the tuple so directory names don’t contain spaces
  dir_tuple="${tuple// /}"

  # Build a single directory name: "<tuple_without_spaces>_<last>"
  target_dir="./${dir_tuple}_${last}"

  mkdir -p "$target_dir"
  mv -- "$file" "$target_dir/"
done
