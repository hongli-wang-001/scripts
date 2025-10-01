#!/bin/bash

base_dir="ens_b"

for i in $(seq -w 1 15); do
    src_index=$(printf "%03d" $((10#$i))) 
    src_dir="${base_dir}/mem${src_index}"
    tgt_index=$(printf "%03d" $((10#$i + 15)))  # Ensure base-10 arithmetic
    tgt_dir="${base_dir}/mem${tgt_index}"

    echo "🔗 Linking mem${i} → mem${tgt_index}"

    # Check source exists
    if [ ! -d "$src_dir" ]; then
        echo "⚠️ Source directory not found: $src_dir"
        continue
    fi

    # Create target directory if not exists
    mkdir -p "$tgt_dir"

    # Create symlinks for all files inside mem${i}
    for file in "$src_dir"/*; do
        filename=$(basename "$file")
        ln -s "$PWD/$file" "$tgt_dir/$filename"
    done

    echo "✅ mem${tgt_index} linked"
    echo "-----------------------------"
done

echo "🎉 All symbolic links created for mem031 to mem060."

