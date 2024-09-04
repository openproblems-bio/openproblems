TASK=dimensionality_reduction
viash run src/common/create_component/config.vsh.yaml -- --task $TASK --type metric --name foor --language r
viash run src/common/create_component/config.vsh.yaml -- --task $TASK --type method --name foor --language r
viash run src/common/create_component/config.vsh.yaml -- --task $TASK --type method --name foopy
viash run src/common/create_component/config.vsh.yaml -- --task $TASK --type metric --name foopy