#!/usr/bin/env python3

# List Available Language Packs
# Used for updating the Languages page of the documentation.
# https://sourcethemes.com/academic/
#
# Prerequisites: pip3 install PyYAML

from pathlib import Path

import yaml

LANG_PATH = Path(__file__).resolve().parent.parent.joinpath("data").joinpath("i18n")
LANG_YAML = LANG_PATH.joinpath("languages.yaml")


# Iterate over languages.
with open(LANG_YAML) as f:
    main_map = yaml.safe_load(f)

    # Print languages as a plaintext list.
    # lang_list = [main_map[main_item] for main_item in main_map]
    # print(', '.join(lang_list))

    # Print languages as a Markdown list.
    i = 0
    for main_item in main_map:
        print(f"- **{main_map[main_item]}** ({main_item})")
        i += 1

print("\n")
print(f"{i} language packs found!")
