#!/usr/bin/env bash
#
# Build PepQueryMHC JAR locally
# Prerequisites: git, Maven (mvn), Java 17+
#
# Usage:
#   ./build_pepquerymhc.sh                 # Builds v1.0.8 (default)
#   ./build_pepquerymhc.sh v1.0.7          # Builds specific version
#

set -euo pipefail

# Configuration
PEPQUERYMHC_VERSION="${1:-v1.0.8}"
REPO_URL="https://github.com/bzhanglab/PepQueryMHC.git"
BUILD_DIR="${2:-.}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Building PepQueryMHC ${PEPQUERYMHC_VERSION}${NC}"

# Check prerequisites
for cmd in git mvn java; do
  if ! command -v "$cmd" &> /dev/null; then
    echo -e "${RED}Error: $cmd is not installed${NC}" >&2
    exit 1
  fi
done

# Verify Java version (need 17+)
JAVA_VERSION=$(java -version 2>&1 | grep -oP 'version "\K[0-9]+' | head -1)
if [[ $JAVA_VERSION -lt 17 ]]; then
  echo -e "${RED}Error: Java 17+ required (found version $JAVA_VERSION)${NC}" >&2
  exit 1
fi

# Create temp directory for build
WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT

echo -e "${YELLOW}Cloning repository...${NC}"
git clone --quiet "$REPO_URL" "$WORK_DIR"

echo -e "${YELLOW}Checking out ${PEPQUERYMHC_VERSION}...${NC}"
cd "$WORK_DIR"
git checkout --quiet "$PEPQUERYMHC_VERSION"

echo -e "${YELLOW}Building JAR with Maven...${NC}"
mvn -q clean package -DskipTests

# Find the built JAR
BUILT_JAR=$(find "$WORK_DIR/target" -name "PepQueryMHC.*.jar" -o -name "PepQueryMHC.jar" | head -1)

if [[ ! -f "$BUILT_JAR" ]]; then
  echo -e "${RED}Error: Build failed - JAR not found${NC}" >&2
  exit 1
fi

# Copy JAR to build directory
OUTPUT_JAR="$BUILD_DIR/PepQueryMHC.jar"
cp "$BUILT_JAR" "$OUTPUT_JAR"

echo -e "${GREEN}✓ Successfully built: $OUTPUT_JAR${NC}"
ls -lh "$OUTPUT_JAR"
