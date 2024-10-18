# import system packages
import datetime
import os
# import the Github module from the PyGithub library
from github import Github
from github import Auth
# import supabase module
from supabase import create_client, Client


## VIASH START
par = {
  "supabase_url": "https://bleficzaoyltozjjndha.supabase.co",
  "supabase_key": os.environ.get("SERVICE_ROLE"),
  "github_auth": os.environ.get("GH_TOKEN"),
  "github_repo": "task_perturbation_prediction",
}

## VIASH END

# authenticate with the Supabase API
print(">> Authenticating with Supabase", flush=True)
url: str = par["supabase_url"]
key: str = par["supabase_key"]
supabase: Client = create_client(url, key)

print(">> Authenticating with Github", flush=True)
auth = Auth.Token(par["github_auth"])
g = Github(auth=auth)

print(">> Getting repo info github", flush=True)
repo = g.get_repo(f"openproblems-bio/{par['github_repo']}")

releases = repo.get_releases()

if releases.totalCount > 0:
  rl_date = releases[0].published_at
  for rl in releases:
    if rl.publshed_at < rl_date:
      rl_date = rl.published_at
    
