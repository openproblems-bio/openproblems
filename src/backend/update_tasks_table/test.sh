#!/usr/bin/bash

set -e

viash run src/backend/update_tasks_table/config.vsh.yaml -- \
  --supabase_url "https://bleficzaoyltozjjndha.supabase.co" \
  --supabase_key $SERVICE_ROLE \
  --github_auth $GH_TOKEN \
  --github_org "openproblems-bio"