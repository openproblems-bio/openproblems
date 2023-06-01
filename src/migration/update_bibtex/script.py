import bibtexparser
from tempfile import NamedTemporaryFile
import urllib.request

## VIASH START
par = {
  'library': 'src/common/library.bib',
  'library_v1': 'https://raw.githubusercontent.com/openproblems-bio/openproblems/main/main.bib'
}
## VIASH END

# Load the BibTeX file
print(">> Read input bibtex file", flush=True)
bib_input = bibtexparser.parse_file(par["library"])

print("  Library keys: " + ', '.join(bib_input.entries_dict.keys()), flush=True)

# Merge with v1 library
if par["library_v1"]:
  print(">> Merge with v1 library", flush=True)
  with NamedTemporaryFile("w", suffix=".bib") as tempfile:
    _ = urllib.request.urlretrieve(par["library_v1"], tempfile.name)
    bib_v1 = bibtexparser.parse_file(tempfile.name)

  print("  Library v1 keys: " + ', '.join(bib_v1.entries_dict.keys()), flush=True)
  blocks = bib_input.blocks + bib_v1.blocks
else:
  blocks = bib_input.blocks

# Remove duplicates
print(">> Remove duplicates", flush=True)
unique_blocks = {block.key : block for block in blocks if not hasattr(block, "error")}
bib_new = bibtexparser.Library(unique_blocks.values())

print("  New keys: " + ', '.join(bib_new.entries_dict.keys()), flush=True)

# Save to a new BibTeX file
print(">> Write to file", flush=True)
bibtexparser.write_file(par["library"], bib_new)
