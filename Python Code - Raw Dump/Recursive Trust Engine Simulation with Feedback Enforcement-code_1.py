import os
import hashlib
import math

CHUNK_SIZE_MB = 8
CHUNK_SIZE_BYTES = CHUNK_SIZE_MB * 1024 * 1024

def sha256_of_text(text):
    return hashlib.sha256(text.encode('utf-8')).hexdigest()

def merge_markdown_files(input_folder, output_base, sort_by='name'):
    files = [f for f in os.listdir(input_folder) if f.endswith('.md')]

    if sort_by == 'name':
        files.sort()
    elif sort_by == 'modified':
        files.sort(key=lambda x: os.path.getmtime(os.path.join(input_folder, x)))

    contents = []
    toc = ["# 🧭 Table of Contents\n"]

    for idx, filename in enumerate(files):
        filepath = os.path.join(input_folder, filename)
        try:
            with open(filepath, 'r', encoding='utf-8') as infile:
                content = infile.read()
        except UnicodeDecodeError:
            with open(filepath, 'r', encoding='latin1') as infile:
                content = infile.read()

        hash_digest = sha256_of_text(content)
        anchor = f"### ⬇ BBP_ANCHOR_{idx}_{filename.replace('.md','')}\n"
        block = (
            f"\n\n<!-- BEGIN {filename} -->\n"
            f"{anchor}\n"
            f"**SHA256**: `{hash_digest}`\n\n"
            f"{content}\n"
            f"\n<!-- END {filename} -->\n"
        )

        contents.append(block)
        toc.append(f"- [{filename}](#bbp_anchor_{idx}_{filename.replace('.md','').lower().replace(' ','-')})")

    full_text = "\n".join(toc) + "\n\n" + "\n".join(contents)

    num_chunks = math.ceil(len(full_text.encode('utf-8')) / CHUNK_SIZE_BYTES)
    chunk_outputs = []

    current_chunk = []
    current_size = 0
    part_number = 1

    for block in contents:
        block_bytes = block.encode('utf-8')
        if current_size + len(block_bytes) > CHUNK_SIZE_BYTES and current_chunk:
            output_path = f"{output_base}_part{part_number}.md"
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write("\n".join(toc) + "\n\n" + "".join(current_chunk))
            print(f"Written: {output_path}")
            chunk_outputs.append(output_path)
            part_number += 1
            current_chunk = []
            current_size = 0

        current_chunk.append(block)
        current_size += len(block_bytes)

    if current_chunk:
        output_path = f"{output_base}_part{part_number}.md"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(toc) + "\n\n" + "".join(current_chunk))
        print(f"Written: {output_path}")
        chunk_outputs.append(output_path)

    print(f"✅ Completed: {len(files)} files merged into {len(chunk_outputs)} part(s).")

# Usage Example
merge_markdown_files(
    input_folder="D:\\Nexus\\Nexus1\\@GPT Downloads\\Chat History\\MD Files",
    output_base="D:\\Nexus\\Nexus1\\@GPT Downloads\\Chat History\\MD Files",
    sort_by='name'
)
