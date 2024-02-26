import yaml

## VIASH START
meta = {
    "config" : "foo"
}
## VIASH END


NAME_MAXLEN = 50

SUMMARY_MAXLEN = 400

DESCRIPTION_MAXLEN = 5000

_MISSING_DOIS = ["vandermaaten2008visualizing", "hosmer2013applied"]

TIME_LABELS = ["lowtime", "midtime", "hightime"]
MEM_LABELS = ["lowmem", "midmem", "highmem"]
CPU_LABELS = ["lowcpu", "midcpu", "highcpu"]

def _load_bib():
    with open(f"{meta['resources_dir']}/library.bib", "r") as file:
        return file.read()

def check_url(url):
    import requests
    from urllib3.util.retry import Retry
    from requests.adapters import HTTPAdapter

    # configure retry strategy
    session = requests.Session()
    retry = Retry(connect=3, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)

    get = session.head(url)

    if get.ok or get.status_code == 429: # 429 rejected, too many requests
        return True
    else:
        return False

def search_ref_bib(reference):
    import re
    bib = _load_bib()
    
    entry_pattern =  r"(@\w+{[^}]*" + reference + r"[^}]*}(.|\n)*?)(?=@)"

    bib_entry = re.search(entry_pattern, bib)

    if bib_entry:

        type_pattern = r"@(.*){" + reference
        doi_pattern = r"(?=[Dd][Oo][Ii]\s*=\s*{([^,}]+)})"

        entry_type = re.search(type_pattern, bib_entry.group(1))

        if not (entry_type.group(1) == "misc" or reference in _MISSING_DOIS):
            entry_doi = re.search(doi_pattern, bib_entry.group(1))
            assert entry_doi.group(1), "doi not found in bibtex reference"
            url = f"https://doi.org/{entry_doi.group(1)}"
            assert check_url(url), f"{url} is not reachable, ref= {reference}."

        return True

    else:
        return False

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
assert "label" in info is not None, "label not an info field or is empty"
assert "summary" in info is not None, "summary not an info field or is empty"
assert "FILL IN:" not in info["summary"], "Summary not filled in"
assert len(info["summary"]) <= SUMMARY_MAXLEN, f"Component id (.functionality.info.summary) should not exceed {SUMMARY_MAXLEN} characters."
assert "description" in info is not None, "description not an info field or is empty"
assert "FILL IN:" not in info["description"], "description not filled in"
assert len(info["description"]) <= DESCRIPTION_MAXLEN, f"Component id (.functionality.info.description) should not exceed {DESCRIPTION_MAXLEN} characters."
if info["type"] == "method":
    assert "reference" in info, "reference not an info field"
    bib = _load_bib()
    if info["reference"]:
        reference = info["reference"]
        if not isinstance(reference, list):
            reference = [reference]
        for ref in reference:
            assert search_ref_bib(ref), f"reference {ref} not added to library.bib"
    assert "documentation_url" in info is not None, "documentation_url not an info field or is empty"
    assert "repository_url" in info is not None, "repository_url not an info field or is empty"
    assert check_url(info["documentation_url"]), f"{info['documentation_url']} is not reachable"
    assert check_url(info["repository_url"]), f"{info['repository_url']} is not reachable"

if "variants" in info:
    arg_names = [arg["name"].replace("--", "") for arg in config["functionality"]["arguments"]] + ["preferred_normalization"]

    for paramset_id, paramset in info["variants"].items():
        if paramset:
            for arg_id in paramset:
                assert arg_id in arg_names, f"Argument '{arg_id}' in `.functionality.info.variants['{paramset_id}']` is not an argument in `.functionality.arguments`."

assert "preferred_normalization" in info, "preferred_normalization not an info field"
norm_methods = ["log_cpm", "log_cp10k", "counts", "log_scran_pooling", "sqrt_cpm", "sqrt_cp10k", "l1_sqrt"]
assert info["preferred_normalization"] in norm_methods, "info['preferred_normalization'] not one of '" + "', '".join(norm_methods) + "'."

print("Check platform fields", flush=True)
platforms = config['platforms']
for platform in platforms:
    if not platform["type"] == "nextflow":
        continue
    nextflow= platform

assert nextflow, "nextflow not a platform"
assert nextflow["directives"], "directives not a field in nextflow platform"
assert nextflow["directives"]["label"], "label not a field in nextflow platform directives"

assert [i for i in nextflow["directives"]["label"] if i in TIME_LABELS], "time label not filled in"
assert [i for i in nextflow["directives"]["label"] if i in MEM_LABELS], "mem label not filled in"
assert [i for i in nextflow["directives"]["label"] if i in CPU_LABELS], "cpu label not filled in"

print("All checks succeeded!", flush=True)
