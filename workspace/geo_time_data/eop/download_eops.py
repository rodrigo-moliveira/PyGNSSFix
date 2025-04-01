import wget
import os

# In this script, the EOP files from https://datacenter.iers.org/data are downloaded
# The TAI-UTC dat file is downloaded from https://maia.usno.navy.mil/ser7/

# 1 - Download finals.all IAU1980 file from https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt
url = "https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt"
# Check if the file exists, and if so, remove it
output_file = "finals1980.all"
if os.path.exists(output_file):
    os.remove(output_file)
wget.download(url, out=output_file)

# 2 - Download finals.all IAU2000 file from https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt
url = "https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt"
output_file = "finals2000.all"
if os.path.exists(output_file):
    os.remove(output_file)
wget.download(url, out=output_file)

# 3 - Download the tai-utc.dat file from https://maia.usno.navy.mil/ser7/tai-utc.dat
url = "https://maia.usno.navy.mil/ser7/tai-utc.dat"
output_file = "tai-utc.dat"
if os.path.exists(output_file):
    os.remove(output_file)
wget.download(url, out=output_file)
