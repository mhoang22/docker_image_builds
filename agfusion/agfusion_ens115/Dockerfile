FROM ubuntu:22.04

# AGFusion for Ensembl release 115
# Uses murphycj/AGFusion and openvax/pyensembl installed from source

ENV DEBIAN_FRONTEND=noninteractive
ENV PYENSEMBL_CACHE_DIR=/opt
ENV ensembl_version=115

# Install ONLY minimal system packages — no python3-pandas, python3-biopython etc.
# to avoid apt/pip mixing which causes the bottleneck/_ARRAY_API conflict
RUN apt-get update -y && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    curl \
    less \
    vim \
    default-libmysqlclient-dev \
    pkg-config \
    git \
    zip \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip first
RUN pip3 install --upgrade pip

# Install ALL Python dependencies via pip only (no apt python packages)
# Pin versions known to be compatible
RUN pip3 install \
    "numpy<2.0" \
    "pandas>=1.5,<3.0" \
    "matplotlib" \
    "biopython==1.79" \
    "future" \
    "mysqlclient" \
    "mysql-connector-python"

# Install pyensembl from GitHub main (to get latest MAX_ENSEMBL_RELEASE)
RUN pip3 install "git+https://github.com/openvax/pyensembl.git@main"

# Patch MAX_ENSEMBL_RELEASE in pyensembl if needed
# Find the file dynamically rather than assuming its name
RUN python3 - <<'EOF'
import os, re, pyensembl

pkg_dir = os.path.dirname(pyensembl.__file__)
print(f"pyensembl package dir: {pkg_dir}")

# Search for whichever file contains MAX_ENSEMBL_RELEASE
for fname in os.listdir(pkg_dir):
    if not fname.endswith('.py'):
        continue
    fpath = os.path.join(pkg_dir, fname)
    content = open(fpath).read()
    if 'MAX_ENSEMBL_RELEASE' in content and re.search(r'MAX_ENSEMBL_RELEASE\s*=\s*\d+', content):
        current = int(re.search(r'MAX_ENSEMBL_RELEASE\s*=\s*(\d+)', content).group(1))
        print(f"Found MAX_ENSEMBL_RELEASE = {current} in {fname}")
        if current < 115:
            patched = re.sub(r'MAX_ENSEMBL_RELEASE\s*=\s*\d+', 'MAX_ENSEMBL_RELEASE = 115', content)
            open(fpath, 'w').write(patched)
            print(f"Patched to 115 in {fname}")
        else:
            print(f"Already >= 115, no patch needed")
        break
else:
    print("WARNING: MAX_ENSEMBL_RELEASE not found in any .py file — pyensembl may not need patching")
EOF

# Install AGFusion from murphycj upstream
WORKDIR /usr/local/bin
RUN git clone https://github.com/murphycj/AGFusion.git

WORKDIR /usr/local/bin/AGFusion

# Patch AGFusion utils.py to add release 115 if not already listed
RUN python3 - <<'EOF'
import re

path = 'agfusion/utils.py'
content = open(path).read()

if '115' in content:
    print("Release 115 already present in utils.py")
else:
    # Add 115 to the end of each species release list
    patched = re.sub(
        r'((?:homo_sapiens|mus_musculus)[^\[]*\[)([\d\s,]+)(\])',
        lambda m: m.group(1) + m.group(2).rstrip().rstrip(',') + ', 115' + m.group(3),
        content
    )
    open(path, 'w').write(patched)
    print("Patched utils.py: added release 115 to species lists")

# Verify
content2 = open(path).read()
print("115 present after patch:", '115' in content2)
EOF

RUN pip3 install .

# Build AGFusion databases
RUN mkdir -p /opt/agfusiondb/
WORKDIR /opt/agfusiondb

# Download Pfam (current release from EBI)
RUN curl -fsSL https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/pfamA.txt.gz -o pfamA.txt.gz && gunzip pfamA.txt.gz

# Build databases — these steps are slow (30-60 min each), network-dependent
RUN agfusion build --dir . --species homo_sapiens --release ${ensembl_version} --pfam pfamA.txt
RUN agfusion build --dir . --species mus_musculus --release ${ensembl_version} --pfam pfamA.txt

RUN rm pfamA.txt

# Install pyensembl reference data (downloads from Ensembl FTP)
RUN pyensembl install --species homo_sapiens --release ${ensembl_version}
RUN pyensembl install --species mus_musculus --release ${ensembl_version}