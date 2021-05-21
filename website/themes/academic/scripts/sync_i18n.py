#!/usr/bin/env python3

"""Sync Language Packs.

Script to synchronize each language pack's items against
Academic's main pack (English). https://sourcethemes.com/academic/

Prerequisites: pip3 install PyYAML

TODO: Switch from PyYAML to Ruamel in order to load/dump comments - see
stackoverflow.com/questions/47382227/python-yaml-update-preserving-order-and-comments
"""

from pathlib import Path

import copy
import yaml

I18N_PATH = Path(__file__).resolve().parent.parent.joinpath("i18n")
MASTER_PACK = I18N_PATH.joinpath("en.yaml")


# Load main language pack (English).
with open(MASTER_PACK) as f:
    main_map = yaml.safe_load(f)
    # if (DEBUG)
    #   print(main_map)

# Iterate over each child language pack.
cnt = 0
for filename in Path(I18N_PATH).glob("*.yaml"):
    if filename.stem != "en":
        i18n_file = I18N_PATH.joinpath(filename)
        print(f"Processing {i18n_file} ...")

        # Load a child language pack.
        with open(i18n_file) as f:
            child_map = yaml.safe_load(f)

        # Synchronize the language pack's structure against the main language pack.
        tmp_map = copy.deepcopy(
            main_map
        )  # Make a temporary deep copy of the main map (list of objects).
        main_index = 0
        for main_item in main_map:
            translation = next(
                (
                    item["translation"]
                    for item in child_map
                    if item["id"] == main_item["id"]
                ),
                main_item["translation"],
            )
            tmp_map[main_index]["translation"] = translation
            main_index += 1

        # Write the synced language pack to file.
        with open(i18n_file, "w") as f:
            yaml.dump(
                tmp_map, f, allow_unicode=True, width=float("inf")
            )  # PyYAML will break lines unless a large column width is set.
        cnt += 1

# Print results.
print(f"{cnt} child language packs successfully synchronized!")
