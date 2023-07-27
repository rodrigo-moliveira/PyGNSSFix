import csv
import numpy as np


class CSVReader:
    @staticmethod
    def read_csv_file(filepath, **options):
        print(f"reading file {filepath}...")
        data = []

        with open(filepath) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')

            if options.get("ignore_header", False):
                next(reader, False)

            for row in reader:
                data.append(row)

        return np.array(data)
