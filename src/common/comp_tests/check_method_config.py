import yaml


## VIASH START

meta = {
    "config" : "foo"
}

## VIASH END

NAME_MAXLEN = 50

SUMMARY_MAXLEN = 400

DESCRIPTION_MAXLEN = 1000


def assert_dict(dict, functionality):

    arg_names = []
    args = functionality["arguments"]

    for i in args:
        arg_names.append(i["name"].replace("--",""))
    
    info = functionality["info"]
    if dict:
        for key in dict:
            assert key in arg_names or info, f"{key} is not a defined argument or .functionality.info field"



print("Load config data", flush=True)
with open(meta["config"], "r") as file:
    config = yaml.safe_load(file)


print("Check general fields", flush=True)
assert len(config["functionality"]["name"]) <= NAME_MAXLEN, f"Component id (.functionality.name) should not exceed {NAME_MAXLEN} characters."
assert "namespace" in config["functionality"] is not None, "namespace not a field or is empty"

print("Check info fields", flush=True)
info = config['functionality']['info']
assert "type" in info, "type not an info field"
info_types = ["method", "control_method"]
assert info["type"] in info_types , f"got {info['type']} expected one of {info_types}"
assert "pretty_name" in info is not None, "pretty_name not an info field or is empty"
assert "summary" in info is not None, "summary not an info field or is empty"
assert "FILL IN:" not in info["summary"], "Summary not filled in"
assert len(info["summary"]) <= SUMMARY_MAXLEN, f"Component id (.functionality.info.summary) should not exceed {SUMMARY_MAXLEN} characters."
assert "description" in info is not None, "description not an info field or is empty"
assert "FILL IN:" not in info["description"], "description not filled in"
assert len(info["description"]) <= DESCRIPTION_MAXLEN, f"Component id (.functionality.info.description) should not exceed {DESCRIPTION_MAXLEN} characters."
if ("control" not in info["type"]):
    assert "reference" in info, "reference not an info field"
    assert "documentation_url" in info is not None, "documentation_url not an info field or is empty"
    assert "repository_url" in info is not None, "repository_url not an info field or is empty"

if "variants" in info:
    for key in info["variants"]:
        assert_dict(info["variants"][key], config['functionality'])

assert "preferred_normalization" in info, "preferred_normalization not an info field"
norm_methods = ["log_cpm", "counts", "log_scran_pooling", "sqrt_cpm", "l1_sqrt"]
assert info["preferred_normalization"] in norm_methods, "info['preferred_normalization'] not one of '" + "', '".join(norm_methods) + "'."



print("All checks succeeded!", flush=True)
