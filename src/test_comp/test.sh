#!/bin/bash

set -ex

echo ">>> Creating dummy input file"
cat > input.txt << HERE
one
two
three
HERE

echo ">>> Running executable"
./test_comp --input input.txt --output output.txt --option FOO

echo ">>> Checking whether output file exists"
[[ ! -f output.txt ]] && echo "Output file could not be found!" && exit 1

# create expected output file
cat > expected_output.txt << HERE
FOO-one
FOO-two
FOO-three
HERE

echo ">>> Checking whether content matches expected content"
diff output.txt expected_output.txt
[ $? -ne 0 ] && echo "Output file did not equal expected output" && exit 1

# print final message
echo ">>> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
