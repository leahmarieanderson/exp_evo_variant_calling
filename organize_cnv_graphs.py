#!/usr/bin/env python3

import os
import shutil
import argparse
import glob

def organize_cnv_graphs(cnv_graphs_dir):
    """
    Organize CNV graphs by moving all PDFs from subdirectories to the main CNV_graphs folder.

    Current structure: CNV_graphs/sample1/sample1.pdf
    Target structure: CNV_graphs/sample1.pdf
    """

    if not os.path.exists(cnv_graphs_dir):
        print(f"Error: Directory {cnv_graphs_dir} does not exist")
        return

    print(f"Organizing CNV graphs in: {cnv_graphs_dir}")

    # Find all PDF files in subdirectories
    pdf_pattern = os.path.join(cnv_graphs_dir, "*", "*.pdf")
    pdf_files = glob.glob(pdf_pattern)

    if not pdf_files:
        print("No PDF files found in subdirectories")
        return

    moved_count = 0

    for pdf_file in pdf_files:
        # Get the filename
        filename = os.path.basename(pdf_file)

        # Target location in the main CNV_graphs directory
        target_path = os.path.join(cnv_graphs_dir, filename)

        # Check if target already exists
        if os.path.exists(target_path):
            print(f"Warning: {filename} already exists in target directory, skipping")
            continue

        # Move the file
        try:
            shutil.move(pdf_file, target_path)
            print(f"Moved: {filename}")
            moved_count += 1
        except Exception as e:
            print(f"Error moving {filename}: {e}")

    print(f"\nSummary: Moved {moved_count} PDF files")

    # Optionally remove empty subdirectories
    for item in os.listdir(cnv_graphs_dir):
        item_path = os.path.join(cnv_graphs_dir, item)
        if os.path.isdir(item_path):
            try:
                os.rmdir(item_path)  # Only removes if empty
                print(f"Removed empty directory: {item}")
            except OSError:
                # Directory not empty, that's fine
                pass

def main():
    parser = argparse.ArgumentParser(description="Organize CNV graphs by moving all PDFs to the main CNV_graphs folder")
    parser.add_argument("cnv_graphs_dir", help="Path to CNV_graphs directory")

    args = parser.parse_args()

    organize_cnv_graphs(args.cnv_graphs_dir)

if __name__ == "__main__":
    main()