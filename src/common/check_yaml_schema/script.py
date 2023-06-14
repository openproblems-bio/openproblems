import jsonschema
import yaml
from pathlib import Path

## VIASH START
par = {
  'input': 'src/tasks/batch_integration/methods/bbknn/config.vsh.yaml',
  'schema': 'src/common/api/schema_task_method.yaml'
}
meta = {
  'functionality_name': 'foo',
}
## VIASH END

def yaml_to_dict(file_path):
    with open(file_path, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def load_schemas(schema_dir):
    schema_files = list(schema_dir.glob("./**/schema_*.yaml"))
    
    schemas = {}
    for file in schema_files:
        schema = yaml_to_dict(file)
        schemas[file.absolute()] = schema
    
    return schemas

def create_validator(schema_name, schemas):
    schema_store = {}
    for name, value in schemas.items():
        schema_store[f"file://{name}"] = value

    # Setting the first schema as the main schema
    
    main_schema = schemas[schema_name]
    resolver = jsonschema.RefResolver(
        base_uri=f"file://{schema_name}",
        referrer=main_schema,
        store=schema_store
    )

    return jsonschema.Draft7Validator(main_schema, resolver=resolver)

print(">> Read input yaml", flush=True)
input_yaml_file = Path(par["input"])
with open(input_yaml_file, 'r') as f:
  input_yaml = yaml.safe_load(f)

print(">> Read schema(s)", flush=True)
schema_yaml_file = Path(par["schema"])
schemas = load_schemas(schema_yaml_file.parent)

print(">> Validate input yaml against schema", flush=True)
validator = create_validator(schema_yaml_file.absolute(), schemas)
validator.validate(input_yaml)
