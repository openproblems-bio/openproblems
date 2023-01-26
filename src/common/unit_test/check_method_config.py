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
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"

print("Check info fields", flush=True)
info = config['functionality']['info']
assert "type" in info, "type not an info field"
info_types = ["method", "negative_control", "positive_control"]
assert info["type"] in info_types , f"got {info['type']} expected one of {info_types}"
assert "method_name" in info, "method_name not an info field"
assert "variants" in info,  "variants not an info field"
assert "preferred_normalization" in info, "preferred_normalization not an info field"
if ("control" not in info["type"]):
    assert "paper_reference" in info, "paper_reference not an info field"



print("All checks succeeded!", flush=True)
