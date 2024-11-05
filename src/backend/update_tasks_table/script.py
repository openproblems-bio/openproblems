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
  "github_org": "openproblems-bio",
}

## VIASH END



# authenticate with the Github API
print(">> Authenticating with Github", flush=True)
auth = Auth.Token(par["github_auth"])
g = Github(auth=auth)


# authenticate with the Supabase API
print(">> Authenticating with Supabase", flush=True)
url: str = par["supabase_url"]
key: str = par["supabase_key"]
supabase: Client = create_client(url, key)

# Get tasks already in db
print(">> Getting tasks from db", flush=True)
tasks = supabase.table("tasks").select("task_name").execute()

task_names=[]
for task in tasks.data:
    task_names.append(task["task_name"])

# get repos from Openproblems-bio organization
print(">> Getting repos from Openproblems-bio organization", flush=True)
repos = g.get_organization(par["github_org"]).get_repos()

# Create rows for tasks that are not in db
print(">> Creating rows for tasks that are not in db", flush=True)
rows= []
for repo in repos:
  if "task_" in repo.name:
    name = repo.name
    created = repo.created_at
    if name not in task_names:
      created_at = created.strftime('%Y-%m-%d %H:%M:%S%z')
      rows.append({"task_name": name, "created_at": created_at})

# Add rows to db tasks table
if len(rows) > 0:
  print(">> Adding rows to db tasks table", flush=True)
  response = (
    supabase.table("tasks")
    .insert(rows)
    .execute()
  )
else:
  print(">> No new tasks to add to db", flush=True)