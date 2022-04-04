# From https://github.community/t/github-actions-new-pulling-from-private-docker-repositories/16089/28
# The goal is to retrieve the ecr password every 6 hours and put it as a secret
from base64 import b64decode

import boto3

def get_ecr_password() -> str:
    """Retrieve ECR password, it comes b64 encoded, in the format user:password
    From https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ecr.html#ECR.Client.get_authorization_token
    """
    ecr = boto3.client("ecr")
    response = ecr.get_authorization_token()
    encoded_login_password = response["authorizationData"][0]["authorizationToken"]

    decoded_login_password = b64decode(encoded_login_password).decode("UTF-8")
    return decoded_login_password.split(":")[1]


if __name__ == "__main__":
    print(get_ecr_password())

