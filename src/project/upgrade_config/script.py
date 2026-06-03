import ruamel.yaml
# import re

## VIASH START
par = {
  "input": "input.vsh.yaml",
  "output": "output.vsh.yaml",
}
## VIASH END

yaml = ruamel.yaml.YAML()

# Set indentation rules
yaml.indent(mapping=2, sequence=4, offset=2) 

# Load input config
with open(par["input"], "r") as file:
    data = yaml.load(file)

transformed_yaml_content = ruamel.yaml.CommentedMap()

# Add __merge__, if necessary
if "__merge__" in data:
    transformed_yaml_content["__merge__"] = data["__merge__"]

# Remove .functionality
if "functionality" in data:
    if "info" in data["functionality"]:
        info_content = data["functionality"]["info"]
        label = info_content.pop("label", None)
        summary = info_content.pop("summary", None)
        description = info_content.pop("description", None)
        reference = info_content.pop("reference", None)
        repository = info_content.pop("repository_url", None)
        documentation = info_content.pop("documentation_url", None)

        # Remove 'info' if it becomes empty
        if not info_content:
            data["functionality"].pop("info")
        
    updated_functionality = ruamel.yaml.CommentedMap()
    updated_functionality["name"] = data["functionality"].pop("name")
    
    # Move out of info
    if label is not None:
        updated_functionality["label"] = label
    if summary is not None:
        updated_functionality["summary"] = summary
    if description is not None:
        updated_functionality["description"] = description
    
    # Fetch doi using reference key
    if reference is not None:
        updated_functionality["references"] = {}
        # with open(f"library.bib", "r") as file:
        #     bib = file.read()
        # entry_pattern =  r"(@\w+{[^}]*" + reference + r"[^}]*}(.|\n)*?)(?=@)"
        # bib_entry = re.search(entry_pattern, bib)
        # if bib_entry:
        #     doi_pattern = r"(?=[Dd][Oo][Ii]\s*=\s*{([^,}]+)})"
        #     entry_doi = re.search(doi_pattern, bib_entry.group(1))
        #     updated_functionality["references"]["doi"] = entry_doi.group(1)
        updated_functionality["references"]["bibtex"] = reference
    
    # Add links
    updated_functionality["links"] = {}
    if repository is not None:
        updated_functionality["links"]["repository"] = repository
    if documentation is not None:
        updated_functionality["links"]["documentation"] = documentation

    # Add remaining contents from .functionality
    updated_functionality.update(data["functionality"])
    
    transformed_yaml_content.update(updated_functionality)

# Mapping platforms to engines and runners
transformed_yaml_content["engines"] = []
transformed_yaml_content["runners"] = []
for platform in data["platforms"]:
    if platform["type"] == "docker":
        transformed_yaml_content["engines"].append(platform)
    elif platform["type"] == "nextflow":
        transformed_yaml_content["runners"].append(platform)

# Insert `type: executable` into runners
transformed_yaml_content["runners"].insert(0, {"type": "executable"})

# Write the transformed YAML to a new file
with open(par["output"], 'w') as file:
    yaml.dump(transformed_yaml_content, file)