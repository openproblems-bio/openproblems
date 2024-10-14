# import system packages
import datetime
import os

import boto3
# import supabase module
from supabase import create_client, Client


## VIASH START
par = {
  "supabase_url": "https://bleficzaoyltozjjndha.supabase.co",
  "supabase_key": os.environ.get("SERVICE_ROLE"),
  "s3_bucket": "openproblems-data",
}

## VIASH END

# authenticate with the Supabase API
print(">> Authenticating with Supabase")
url: str = par["supabase_url"]
key: str = par["supabase_key"]
supabase: Client = create_client(url, key)

print(">> Getting tasks from db")
tasks = supabase.table("tasks").select("task_name").execute()
task_names=[]
for task in tasks.data:
    task_names.append(task["task_name"])

prefix = f"resources/denoising/results/"

# Let's use Amazon S3
s3 = boto3.client('s3')
response = s3.list_buckets()

objects = s3.list_objects_v2(Bucket=par["s3_bucket"], Prefix=prefix)