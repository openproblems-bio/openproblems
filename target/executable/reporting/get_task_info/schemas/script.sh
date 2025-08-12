#!/bin/bash

set -e

VIASH_VERSION=0.9.0 viash export json_schema -f yaml > schemas/schema_viash.yaml