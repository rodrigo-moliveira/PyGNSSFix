import json
from jsonschema import validate, ValidationError


def validate_config(json_data, schema):
    try:
        validate(json_data, schema)
        print("Validation successful.")
    except ValidationError as e:
        print(f"Validation failed: {e}")


# Read the schema from the file
with open('./resources/schema.json') as schema_file:
    schema = json.load(schema_file)

# Example usage
with open('./resources/gnss_sample.json') as config_file:
    config_data = json.load(config_file)

validate_config(config_data, schema)
