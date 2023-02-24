import sys
import json
# do Imports as
from src.io.config.gnss_config import ConfigGNSS


def main():
    # Opening JSON file
    filename = "C:\\Users\\rooo\\Documents\\other_projects\\PyGNSSFix\\src\\io\\config\\resources\\gnss_sample.json"
    with open(filename) as json_file:
        data = json.load(json_file)
    config = ConfigGNSS(**data)


print("#--------------------------------------------------#")
print("#           Welcome to PyGNSSFix Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
