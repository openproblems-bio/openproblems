
from openproblems.project import read_viash_config, check_config, run_and_check_output

## VIASH START
meta = {
  # ...
}
## VIASH END

# read viash config
config = read_viash_config(meta["config"])

# check whether the config is valid and contains all required fields
check_config(config)

# run the component with test arguments and check the output
run_and_check_output(meta, config)
