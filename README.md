## Installation

### Using conda

Enforcing python 3.11 is recommended for ensuring working antiSMASH installation on all platforms.
```bash
conda create -n poliwag python=3.11 bioconda::antismash
conda install -c bioconda ncbi-datasets-cli
conda activate poliwag
pip install paras --no-deps
pip install poliwag
```