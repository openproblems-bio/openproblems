# import system packages
import datetime
import os

import boto3
# import supabase module
from supabase import create_client, Client


# Let's use Amazon S3
s3 = boto3.resource('s3')