# Pipeline Diagrams - How to Generate SVG Images

This document explains how to generate SVG images from the Mermaid diagram files for the three workflows.

## Files

- `pipeline_shortreads.mmd` - Short-read workflow diagram
- `pipeline_longreads.mmd` - Long-read workflow diagram  
- `pipeline_assemfree.mmd` - Assembly-free workflow diagram

## Method 1: Using Mermaid CLI (Recommended)

### Installation

Install Mermaid CLI using npm:

```bash
npm install -g @mermaid-js/mermaid-cli
```

Or using conda:

```bash
conda install -c conda-forge mermaid-cli
```

### Generate SVG

Generate SVG files from the Mermaid files:

```bash
# Generate SVG for short-read workflow
mmdc -i docs/pipeline_shortreads.mmd -o docs/images/pipeline_shortreads.svg

# Generate SVG for long-read workflow
mmdc -i docs/pipeline_longreads.mmd -o docs/images/pipeline_longreads.svg

# Generate SVG for assembly-free workflow
mmdc -i docs/pipeline_assemfree.mmd -o docs/images/pipeline_assemfree.svg
```

### Options

You can customize the output with additional options:

```bash
# Generate with custom width and background color
mmdc -i docs/pipeline_shortreads.mmd -o docs/images/pipeline_shortreads.svg \
  -w 2400 -H 1800 -b white

# Generate PNG instead of SVG
mmdc -i docs/pipeline_shortreads.mmd -o docs/images/pipeline_shortreads.png \
  -w 2400 -H 1800
```

## Method 2: Using Online Mermaid Editor

1. Go to [Mermaid Live Editor](https://mermaid.live/)
2. Copy the content from one of the `.mmd` files
3. Paste it into the editor
4. Click "Actions" → "Download SVG" or "Download PNG"

## Method 3: Using GitHub/GitLab

If your repository is hosted on GitHub or GitLab, you can:

1. Create a markdown file (e.g., `docs/pipeline_diagrams.md`)
2. Include the Mermaid diagram using code blocks:

````markdown
```mermaid
%%{init: {'theme':'base', 'themeVariables': { 'primaryColor':'#1ba1e2', 'primaryTextColor':'#fff', 'primaryBorderColor':'#006EAF', 'lineColor':'#1ba1e2', 'secondaryColor':'#6c8ebf', 'tertiaryColor':'#d5e8d4'}}}%%
flowchart LR
    ...
```
````

3. GitHub/GitLab will automatically render the diagram
4. Right-click on the rendered diagram and "Save image as..." to download as SVG

## Method 4: Using Docker

If you have Docker installed:

```bash
# Pull the Mermaid CLI Docker image
docker pull minlag/mermaid-cli

# Generate SVG
docker run --rm -v $(pwd):/data minlag/mermaid-cli \
  -i /data/docs/pipeline_shortreads.mmd \
  -o /data/docs/images/pipeline_shortreads.svg
```

## Method 5: Using Python (mermaid-py)

Install mermaid-py:

```bash
pip install mermaid-py
```

Generate SVG:

```python
from mermaid import Mermaid

# Read the Mermaid file
with open('docs/pipeline_shortreads.mmd', 'r') as f:
    mermaid_code = f.read()

# Generate SVG
mermaid = Mermaid(mermaid_code)
svg = mermaid.to_svg()

# Save to file
with open('docs/images/pipeline_shortreads.svg', 'w') as f:
    f.write(svg)
```

## Batch Generation Script

Create a simple bash script to generate all diagrams at once:

```bash
#!/bin/bash
# generate_diagrams.sh

# Create output directory if it doesn't exist
mkdir -p docs/images

# Generate all SVG diagrams
mmdc -i docs/pipeline_shortreads.mmd -o docs/images/pipeline_shortreads.svg -w 2400 -H 1800
mmdc -i docs/pipeline_longreads.mmd -o docs/images/pipeline_longreads.svg -w 2400 -H 1800
mmdc -i docs/pipeline_assemfree.mmd -o docs/images/pipeline_assemfree.svg -w 2400 -H 1800

echo "All diagrams generated successfully!"
```

Make it executable and run:

```bash
chmod +x generate_diagrams.sh
./generate_diagrams.sh
```

## Troubleshooting

### Issue: mmdc command not found

**Solution**: Make sure Mermaid CLI is installed and in your PATH:
```bash
npm install -g @mermaid-js/mermaid-cli
```

### Issue: Diagrams are too small or cut off

**Solution**: Increase the width and height parameters:
```bash
mmdc -i docs/pipeline_shortreads.mmd -o docs/images/pipeline_shortreads.svg -w 3200 -H 2400
```

### Issue: Font rendering issues

**Solution**: Install required fonts or use a different theme. You can modify the theme in the `.mmd` files.

## Notes

- SVG format is recommended for documentation as it scales without quality loss
- PNG format can be used if SVG is not supported
- The diagrams use nf-core color scheme (blue theme)
- All diagrams are in landscape orientation (LR - left to right)

