#!/bin/bash

function pause(){
   read -n 1 -p "Press any key to continue..."
}

GIT_USER="$1"

# If no argument was provided, ask for it interactively
if [ -z "$GIT_USER" ]; then
    read -p "Please enter your GitHub Username: " GIT_USER
fi

if [ -z "$GIT_USER" ]; then
    echo "❌ Error: Username cannot be blank."
    exit 1
fi

echo "Starting workspace configuration for user: ${GIT_USER}"

# --- 1. Pull Support Libraries ---
declare -a arr=("zlib" "bzip2" "openssl" "hdf5" "cython" "numpy" "h5py" "lemon" "lmx" "rapidjson" "antlr" "fmt" "icu4c" "log4cplus" "boost" "pion" "saxon-c" "gemmi" "mmdb2" "libccp4" "clipper" )

for i in "${arr[@]}"
do
   echo "Fetching support library: $i"
   if ! [ -d "${i}" ] ; then
      # Note: 'localhost' automatically works here because of our compose 'extra_hosts' mapping!
      #wget --no-check-certificate http://ci.quantumbioinc.com:8081/view/darwin/job/${i}-darwin/lastSuccessfulBuild/artifact/${i}/*zip*/${i}.zip
      wget --no-check-certificate http://ci.quantumbioinc.com:8081/view/rhel-build64/job/${i}/lastSuccessfulBuild/artifact/${i}/*zip*/${i}.zip
      #wget --no-check-certificate http://roger:plum@localhost:8080/job/${i}/lastSuccessfulBuild/artifact/${i}/*zip*/${i}.zip
      #curl -o ${i}.zip -0 http://roger:plum@localhost:8080/job/${i}/lastSuccessfulBuild/artifact/${i}/*zip*/${i}.zip
      unzip -o ${i}.zip
      # 🛠️ JENKINS SYMLINK REPAIR WORKAROUND 🛠️
      # Scan the unzipped folder for any concrete shared libraries with version numbers
      # and dynamically link them to the base .so name NetBeans needs.
      if [ -d "${i}" ]; then
        # Rename the directory if the target folder doesn't exist yet
        if [ ! -d "antlr4" ]; then
            mv antlr antlr4
        fi

        # Step into the correct antlr4 directory to fix lib64
        cd antlr4
        if [ -d "lib64" ]; then
            echo "🗂️ Flattening lib64 into standard lib folder..."
            mkdir -p lib
            cp -rp lib64/* lib/
            rm -rf lib64
        fi
         echo "Checking ${i} for unzipped shared libraries needing symlink repair..."
         
         # Find files matching *.so.* (like libclipper-core.so.2.0.1)
         find "${i}" -type f -name "*.so.*" | while read -r real_so_file; do
            # Extract the folder path and the file name
            so_dir=$(dirname "$real_so_file")
            so_filename=$(basename "$real_so_file")
            
            # Isolate the base name (e.g., cut libclipper-core.so.2.0.1 down to libclipper-core.so)
            base_so_name=$(echo "$so_filename" | sed -E 's/(\.so)(\.[0-9]+)+$/\1/')
            
            # Check if the clean base .so file name is missing in that directory
            if [ ! -f "${so_dir}/${base_so_name}" ] && [ ! -L "${so_dir}/${base_so_name}" ]; then
               echo "🔗 Restoring symlink: ${base_so_name} -> ${so_filename}"
               # Run the symlink creation inside its matching library directory
               (cd "$so_dir" && ln -sf "$so_filename" "$base_so_name")
            fi

            # 2. NEW: Generate the intermediate SONAME link (e.g., libclipper-minimol.so.2)
            # This captures just the first digit after the .so
            major_so_name=$(echo "$so_filename" | sed -E 's/(\.so\.[0-9]+).*/\1/')
            if [ "$major_so_name" != "$so_filename" ]; then
               if [ ! -f "${so_dir}/${major_so_name}" ] && [ ! -L "${so_dir}/${major_so_name}" ]; then
                  echo "🔗 Restoring intermediate SONAME link: ${major_so_name} -> ${so_filename}"
                  (cd "$so_dir" && ln -sf "$so_filename" "$major_so_name")
               fi
            fi
         done
      fi
   fi
done

pause
if ! [ -d "eigen" ] ; then
mv eigeneigen/include/eigen3 eigen
rm eigeneigen
fi
pause
# --- 2. Configure Git + Dynamic GHP Token Mapping ---
# Use the dynamic username variable here
git config --global user.name "${GIT_USER}"

# Automatically intercept GitHub authentication headers using your ghp_ variable
if [ -n "$GITHUB_TOKEN" ]; then
    git config --global url."https://x-access-token:${GITHUB_TOKEN}@github.com/".insteadOf "https://github.com/"
    # If your original URL format includes ${GIT_USER}@, explicitly match it:
    git config --global url."https://x-access-token:${GITHUB_TOKEN}@github.com/quantumbio/".insteadOf "https://${GIT_USER}@github.com/quantumbio/"
fi

# --- 3. Pull GitHub Repositories ---
declare -a projectarr=("FockianIntegrals" "Persistence" "OOBackbone" "ExpDensity" "JPRoothaan" "OBMM" "Solvation" "SpeciesService" "libMovableType" "qbdiff" "QMechanic")
declare -a projectbranch=("develop" "develop" "develop" "develop" "develop" "develop" "develop" "develop" "develop" "develop" "develop")

ii=0
for i in "${projectarr[@]}"
do
   echo "Cloning repo: $i [Branch: ${projectbranch[$ii]}]"
   if ! [ -d "${i}" ] ; then
      git clone -b ${projectbranch[$ii]} https://${GIT_USER}@github.com/quantumbio/${i}.git ${i}
      #cd ${i}
      #git reset --hard $(git rev-list -1 $(git rev-parse --until=2024-04-11) ${projectbranch[$ii]})
      #cd ..
   fi
   ii=$((ii+1))
done

echo "✅ Workspace initialization complete!"
