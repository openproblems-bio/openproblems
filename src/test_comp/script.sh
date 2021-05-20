#!/bin/bash

echo "This is a skeleton component"
echo "The arguments are:"
echo " - input:  $par_input"
echo " - output: $par_output"
echo " - option: $par_option"
echo

echo "Writing output file"
cat "$par_input" | sed "s#.*#$par_option-&#" > "$par_output"
