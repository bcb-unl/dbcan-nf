#!/bin/bash
# Generate SVG diagrams from Mermaid files
# Usage: ./generate_diagrams.sh [format]
# Format options: svg (default), png, pdf

set -e

FORMAT=${1:-svg}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DOCS_DIR="$PROJECT_ROOT/docs"
IMAGES_DIR="$DOCS_DIR/images"

# Create images directory if it doesn't exist
mkdir -p "$IMAGES_DIR"

# Check if mmdc is installed
if ! command -v mmdc &> /dev/null; then
    echo "Error: Mermaid CLI (mmdc) is not installed."
    echo "Please install it using one of the following methods:"
    echo "  npm install -g @mermaid-js/mermaid-cli"
    echo "  conda install -c conda-forge mermaid-cli"
    exit 1
fi

echo "Generating $FORMAT diagrams..."
echo "Output directory: $IMAGES_DIR"
echo ""

# Generate diagrams with appropriate dimensions
WIDTH=2400
HEIGHT=1800

case $FORMAT in
    svg)
        EXT="svg"
        ;;
    png)
        EXT="png"
        ;;
    pdf)
        EXT="pdf"
        ;;
    *)
        echo "Error: Unsupported format '$FORMAT'. Use svg, png, or pdf."
        exit 1
        ;;
esac

# Generate short-read workflow diagram
echo "Generating short-read workflow diagram..."
mmdc -i "$DOCS_DIR/pipeline_shortreads.mmd" \
     -o "$IMAGES_DIR/pipeline_shortreads.$EXT" \
     -w $WIDTH -H $HEIGHT -b white

# Generate long-read workflow diagram
echo "Generating long-read workflow diagram..."
mmdc -i "$DOCS_DIR/pipeline_longreads.mmd" \
     -o "$IMAGES_DIR/pipeline_longreads.$EXT" \
     -w $WIDTH -H $HEIGHT -b white

# Generate assembly-free workflow diagram
echo "Generating assembly-free workflow diagram..."
mmdc -i "$DOCS_DIR/pipeline_assemfree.mmd" \
     -o "$IMAGES_DIR/pipeline_assemfree.$EXT" \
     -w $WIDTH -H $HEIGHT -b white

echo ""
echo "Successfully generated all $FORMAT diagrams in $IMAGES_DIR"
echo ""
echo "Generated files:"
ls -lh "$IMAGES_DIR"/*.$EXT 2>/dev/null || echo "No $EXT files found"

