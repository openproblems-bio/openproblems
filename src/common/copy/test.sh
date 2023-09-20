#!/usr/bin/env bash
set -ex

touch test_file.txt

echo ">>> Testing if publish in local dir works"
"$meta_executable" \
  --input test_file.txt \
  --output another_file.txt

[[ ! -f another_file.txt ]] && echo "It seems no output is generated" && exit 1

echo ">>> Testing if publish in local dir works"
"$meta_executable" \
  --input test_file.txt \
  --output adir/yadir/another_file.txt

[[ ! -d adir ]] && echo "It seems no output is generated" && exit 1
[[ ! -f adir/yadir/another_file.txt ]] && echo "It seems no output is generated" && exit 1

echo ">>> Test finished successfully"
