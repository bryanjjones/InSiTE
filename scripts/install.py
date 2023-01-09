# First, make sure that pip is installed
import subprocess
import os
import sys

subprocess.run(["pip3", "install", "pip"])
subprocess.run(["pip3", "install", "--upgrade", "pip"])

# Now, use pip to install the packages
subprocess.run(["pip3", "install", "biopython"])
subprocess.run(["pip3", "install", "nrel-pysam"])
subprocess.run(["pip3", "install", "twobitreader"])
subprocess.run(["pip3", "install", "pybedtools"])

# Add the path to the scripts to the system PATH
# Get the path to the scripts folder
scripts_path = os.path.join(sys.prefix, "Scripts")

# Add the scripts folder to the system PATH
os.environ["PATH"] += os.pathsep + scripts_path

# Print the updated PATH
print(os.environ["PATH"])
