from json_schema_for_humans.generate import generate_from_filename
from json_schema_for_humans.generation_configuration import GenerationConfiguration

config = GenerationConfiguration(copy_css=False, expand_buttons=True)

generate_from_filename("./src/io/config/resources/performance_schema.json", "configs/performance_schema_doc.html", config=config)