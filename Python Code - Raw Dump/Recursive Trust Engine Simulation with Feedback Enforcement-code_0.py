import os

def merge_markdown_files(input_folder, output_file, sort_by='name'):
    files = [f for f in os.listdir(input_folder) if f.endswith('.md')]

    # Optional sorting
    if sort_by == 'name':
        files.sort()
    elif sort_by == 'modified':
        files.sort(key=lambda x: os.path.getmtime(os.path.join(input_folder, x)))

    with open(output_file, 'w', encoding='utf-8') as outfile:
        for filename in files:
            filepath = os.path.join(input_folder, filename)
            try:
                with open(filepath, 'r', encoding='utf-8') as infile:
                    content = infile.read()
            except UnicodeDecodeError:
                # Fallback for incompatible encodings (like Windows-1252 or Latin-1)
                with open(filepath, 'r', encoding='latin1') as infile:
                    content = infile.read()

            outfile.write(f"\n\n<!-- BEGIN {filename} -->\n\n")
            outfile.write(content)
            outfile.write(f"\n\n<!-- END {filename} -->\n\n")

    print(f"Merged {len(files)} files into {output_file}")

# Usage
merge_markdown_files("D:\\Nexus\\Obsidian Valut\\Nexus3\\", "merged_output.md", sort_by='name')
