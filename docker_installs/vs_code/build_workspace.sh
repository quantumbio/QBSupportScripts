#!/bin/bash

declare -a projects=("FockianIntegrals" "Persistence" "OOBackbone" "JPRoothaan" "OBMM" "ExpDensity" "Solvation" "SpeciesService" "libMovableType" "QMechanic")

for i in "${projects[@]}"
do
   echo "=========================================="
   echo "Building Project: $i"
   echo "=========================================="
   
   # Check if the folder actually exists before trying to enter it
   if [ ! -d "${i}" ]; then
      echo "❌ Error: Directory '${i}' does not exist!"
      exit 1
   fi

   cd "${i}"

   # Run the NetBeans Makefile targets. 
   # Using 'make clean' followed by 'make' ensures a reliable build order.
   make -f Makefile CONF=Release clean
   make -f Makefile CONF=Release
   
   # 🚨 CRITICAL SAFETY CHECK 🚨
   # Capture the exit code of the 'make' execution
   COMPILE_STATUS=$?
   
   if [ $COMPILE_STATUS -ne 0 ]; then
      echo "❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌"
      echo "Compilation FAILED in project: $i"
      echo "Stopping workspace build execution immediately."
      echo "❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌ ❌"
      exit 1
   fi

   cd ..
done

echo "✅ Success! All sibling projects compiled flawlessly."
