#!/usr/bin/env bash

# --- 1. Strict Error Trap Settings ---
# Exit immediately if a command exits with a non-zero status
set -e
# Ensure variables are defined before use; stop execution on a pipeline failure
set -u
set -o pipefail

function pause(){
   read -n 1 -s -p "Press any key to continue..."
   echo ""
}

GIT_USER="${1:-}"

# If no argument was provided, ask for it interactively
if [ -z "$GIT_USER" ]; then
    read -p "Please enter your GitHub Username: " GIT_USER
fi

if [ -z "$GIT_USER" ]; then
    echo "❌ Error: Username cannot be blank." >&2
    exit 1
fi

echo "Starting workspace configuration for user: ${GIT_USER}"

# --- 2. Pull Support Libraries ---
declare -a arr=("zlib" "bzip2" "openssl" "hdf5" "cython" "numpy" "h5py" "lemon" "lmx" "rapidjson" "antlr" "cif" "icu4c" "log4cplus" "boost" "pion" "fmt" "gemmi" "mmdb2" "libccp4" "clipper" )

for i in "${arr[@]}"
do
   if [ -d "${i}" ] ; then
      echo "ℹ️  Support library [${i}] already exists. Skipping download."
      continue
   fi

   echo "🚀 Fetching support library: ${i}"
   
   # Download to a targeted zip name to avoid collisions
   ZIP_FILE="${i}_download.zip"
   
   # Using a subshell loop block to capture failures elegantly
   if ! wget --no-check-certificate -q --show-progress \
      "http://ci.quantumbioinc.com:8081/view/rhel-build64/job/${i}/lastSuccessfulBuild/artifact/${i}/*zip*/${i}.zip" \
      -O "${ZIP_FILE}"; then
         echo "❌ Error: Failed to download archive for ${i}." >&2
         rm -f "${ZIP_FILE}"
         exit 1
   fi

   # Force unzip and instantly clean up the file chunk
   if ! unzip -q -o "${ZIP_FILE}"; then
         echo "❌ Error: Failed to extract zip archive for ${i}." >&2
         rm -f "${ZIP_FILE}"
         exit 1
   fi
   rm -f "${ZIP_FILE}"

      # 🛠️ JENKINS SYMLINK REPAIR WORKAROUND 🛠️
      if [ -d "${i}" ]; then
      echo "⚙️  Checking ${i} for unzipped shared libraries needing symlink repair..."
         
      # Find files matching *.so.* (like libclipper-core.so.2.0.1) Safely handles spaces
      find "${i}" -type f -name "*.so.*" -print0 | while read -r -d '' real_so_file; do
            so_dir=$(dirname "$real_so_file")
            so_filename=$(basename "$real_so_file")
            
         # Isolate the base name (e.g., cuts down to libclipper-core.so)
            base_so_name=$(echo "$so_filename" | sed -E 's/(\.so)(\.[0-9]+)+$/\1/')
            
            # Check if the clean base .so file name is missing in that directory
            if [ ! -f "${so_dir}/${base_so_name}" ] && [ ! -L "${so_dir}/${base_so_name}" ]; then
               echo "🔗 Restoring symlink: ${base_so_name} -> ${so_filename}"
            # Run the symlink creation inside its matching library directory safely
               (cd "$so_dir" && ln -sf "$so_filename" "$base_so_name")
            fi
         done
      fi
done

pause

# --- 3. Configure Git + Dynamic GHP Token Mapping ---
git config --global user.name "${GIT_USER}"

if [ -n "${GITHUB_TOKEN:-}" ]; then
    echo "🔒 Applying secure GitHub credential intercept maps..."
    git config --global url."https://x-access-token:${GITHUB_TOKEN}@github.com/".insteadOf "https://github.com/"
    git config --global url."https://x-access-token:${GITHUB_TOKEN}@github.com/quantumbio/".insteadOf "https://${GIT_USER}@github.com/quantumbio/"
fi

# --- 4. Pull GitHub Repositories (Safe Multi-Value Array Map) ---
# We bundle the project and its branch explicitly inside a single string loop.
# Format: "ProjectName:BranchName"
declare -a projects=(
    "FockianIntegrals:develop"
    "Persistence:develop"
    "OOBackbone:develop"
    "JPRoothaan:develop"
    "OBMM:develop"
    "ExpDensity:develop"
    "Solvation:develop"
    "SpeciesService:develop"
    "libMovableType:develop"
    "qbdiff:develop"
    "QMechanic:develop"
)

for project_entry in "${projects[@]}"
do
   # Parse the project name and branch using bash pattern matching
   repo="${project_entry%%:*}"
   branch="${project_entry#*:}"

   if [ -d "${repo}" ]; then
      echo "ℹ️  Repository [${repo}] already cloned. Skipping."
      continue
   fi

   echo "📥 Cloning repo: ${repo} [Branch: ${branch}]"
   if ! git clone -b "${branch}" "https://${GIT_USER}@github.com/quantumbio/${repo}.git" "${repo}"; then
      echo "❌ Error: Failed to clone repository ${repo}." >&2
      exit 1
   fi
done

echo "✅ Workspace initialization completely successful and verified!"
