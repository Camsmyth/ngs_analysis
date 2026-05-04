#!/usr/bin/env bash
# =============================================================================
# setup.sh — ngs_analysis environment installer
# Usage:  bash setup.sh
# =============================================================================
# Creates a virtual environment named 'ngs' inside this directory,
# installs all Python dependencies, and verifies critical imports.
# Tested on: macOS (arm64/x86_64), Ubuntu 22.04/24.04
# Requires:  Python 3.10+, pip, (optionally) hmmer for ANARCI full functionality
# =============================================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_NAME="ngs"
VENV_PATH="$REPO_DIR/$VENV_NAME"
PYTHON_MIN="3.10"

# ── Colours ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
info()  { echo -e "${GREEN}[setup]${NC} $*"; }
warn()  { echo -e "${YELLOW}[warn] ${NC} $*"; }
error() { echo -e "${RED}[error]${NC} $*"; exit 1; }

echo ""
echo "============================================"
echo "  ngs_analysis — environment setup"
echo "============================================"
echo ""

# ── 1. Locate Python ──────────────────────────────────────────────────────────
PYTHON=""
for candidate in python3.12 python3.11 python3.10 python3 python; do
    if command -v "$candidate" &>/dev/null; then
        ver=$("$candidate" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>/dev/null || echo "0.0")
        major=${ver%%.*}; minor=${ver##*.}
        req_major=${PYTHON_MIN%%.*}; req_minor=${PYTHON_MIN##*.}
        if [ "$major" -gt "$req_major" ] || { [ "$major" -eq "$req_major" ] && [ "$minor" -ge "$req_minor" ]; }; then
            PYTHON="$candidate"
            break
        fi
    fi
done

[ -z "$PYTHON" ] && error "Python $PYTHON_MIN or higher not found. Install from https://python.org"
info "Using Python: $($PYTHON --version)  →  $(command -v $PYTHON)"

# ── 2. Create venv ────────────────────────────────────────────────────────────
if [ -d "$VENV_PATH" ]; then
    warn "venv '$VENV_NAME' already exists at $VENV_PATH — skipping creation."
    warn "To rebuild: rm -rf $VENV_PATH && bash setup.sh"
else
    info "Creating virtual environment: $VENV_PATH"
    "$PYTHON" -m venv "$VENV_PATH"
fi

PIP="$VENV_PATH/bin/pip"
PYTHON_VENV="$VENV_PATH/bin/python"

# ── 3. Upgrade pip ────────────────────────────────────────────────────────────
info "Upgrading pip..."
"$PIP" install --quiet --upgrade pip

# ── 4. Install requirements ───────────────────────────────────────────────────
REQ="$REPO_DIR/requirements.txt"
[ -f "$REQ" ] || error "requirements.txt not found at $REQ"
info "Installing dependencies from requirements.txt..."
"$PIP" install --quiet -r "$REQ"
info "Dependencies installed."

# ── 5. Optional system dependency: HMMER (for ANARCI full functionality) ──────
echo ""
if command -v hmmbuild &>/dev/null; then
    info "HMMER detected: $(hmmbuild -h 2>&1 | head -2 | tail -1 | xargs)"
else
    warn "HMMER not found. ANARCI will fall back to offset-based CDR extraction."
    warn "Install HMMER for full IMGT numbering support:"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        warn "  brew install hmmer"
    else
        warn "  sudo apt install hmmer   # Ubuntu/Debian"
        warn "  conda install -c bioconda hmmer   # Conda"
    fi
fi

# ── 6. Verify critical imports ────────────────────────────────────────────────
echo ""
info "Verifying imports..."
"$PYTHON_VENV" - <<'PYCHECK'
import importlib, sys
required = [
    ("Bio",        "biopython"),
    ("pysam",      "pysam"),
    ("anarci",     "anarci"),
    ("regex",      "regex"),
    ("scipy",      "scipy"),
    ("statsmodels","statsmodels"),
    ("numpy",      "numpy"),
    ("pandas",     "pandas"),
    ("matplotlib", "matplotlib"),
    ("openpyxl",   "openpyxl"),
    ("rich",       "rich"),
    ("tqdm",       "tqdm"),
    ("rapidfuzz",  "rapidfuzz"),
    ("networkx",   "networkx"),
]
failed = []
for mod, pkg in required:
    try:
        importlib.import_module(mod)
        print(f"  \033[0;32m✓\033[0m  {pkg}")
    except ImportError:
        print(f"  \033[0;31m✗\033[0m  {pkg}  ← MISSING")
        failed.append(pkg)
if failed:
    print(f"\n\033[0;31m[error] Missing packages: {', '.join(failed)}\033[0m")
    sys.exit(1)
PYCHECK

# ── 7. Done ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo -e "  ${GREEN}Setup complete!${NC}"
echo "============================================"
echo ""
echo "  Activate the environment:"
echo "    source $VENV_PATH/bin/activate"
echo ""
echo "  Run the full pipeline:"
echo "    python run_pipeline.py --help"
echo ""
echo "  Run individual steps:"
echo "    python bam_extract.py <bam_dir> --min-q 12 --cdr-method anarci"
echo "    python cluster_levenshtein.py --input <csv> --threshold 0.85"
echo "    python cluster_enrichment.py <R1_consensus.csv> <R2_consensus.csv>"
echo ""
