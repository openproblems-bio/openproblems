
from openproblems.project import read_viash_config, check_config

## VIASH START
meta = {
  # ...
}
## VIASH END

# read viash config
config = read_viash_config(meta["config"])

# check whether the config is valid and contains all required fields
check_config(config)
