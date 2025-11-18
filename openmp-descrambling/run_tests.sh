#!/bin/bash

# Thread counts to test
THREADS=(1 2 4 8)

# Input file to encode
INPUT_FILE="Frankenstein.txt"

# Output directory
OUTDIR="results"
mkdir -p "$OUTDIR"

echo "Running encoder..."
for t in "${THREADS[@]}"; do
    echo "  Encoding with $t threads..."
    
    # Run encoder (assumes it writes encoded_<t>_threads.txt internally)
    ./thread_file_encoder "$t" "$INPUT_FILE"
    
    ENCODED_FILE="encoded_${t}_threads.txt"
    
    # Move the generated file into results/
    if [ -f "$ENCODED_FILE" ]; then
        mv "$ENCODED_FILE" "$OUTDIR/$ENCODED_FILE"
    else
        echo "Error: encoder did not produce $ENCODED_FILE"
    fi
done

echo "Running decoder..."
for t in "${THREADS[@]}"; do
    echo "  Decoding with $t threads..."
    
    ENCODED_FILE="$OUTDIR/encoded_${t}_threads.txt"
    
    # Run decoder (assumes it writes decoded_<t>_threads.txt internally)
    ./thread_decoder "$t" "$ENCODED_FILE"
    
    DECODED_FILE="decoded_${t}_threads.txt"
    
    # Move the generated file into results/
    if [ -f "$DECODED_FILE" ]; then
        mv "$DECODED_FILE" "$OUTDIR/$DECODED_FILE"
    else
        echo "Error: decoder did not produce $DECODED_FILE"
    fi
done

echo "Done! All files are in $OUTDIR/"