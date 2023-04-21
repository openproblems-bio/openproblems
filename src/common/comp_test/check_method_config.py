import yaml


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

print("Load config data", flush=True)
with open(meta["config"], "r") as file:
                config = yaml.safe_load(file)


print("check general fields", flush=True)
assert "name" in config["functionality"] is not None, "Name not a field or is empty"
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"

print("Check info fields", flush=True)
info = config['functionality']['info']
assert "type" in info, "type not an info field"
info_types = ["method", "control_method"]
assert info["type"] in info_types , f"got {info['type']} expected one of {info_types}"
assert "pretty_name" in info is not None, "pretty_name not an info field or is empty"
assert "summary" in info is not None, "summary not an info field or is empty"
assert "description" in info is not None, "description not an info field or is empty"
if ("control" not in info["type"]):
    assert "reference" in info, "reference not an info field"
    assert "documentation_url" in info is not None, "documentation_url not an info field or is empty"
    assert "repository_url" in info is not None, "repository_url not an info field or is empty"
assert "variants" in info,  "variants not an info field"
assert "preferred_normalization" in info, "preferred_normalization not an info field"




print("All checks succeeded!", flush=True)
