# From https://github.community/t/github-actions-new-pulling-from-private-docker-repositories/16089/28
# The goal is to retrieve the ecr password every 6 hours and put it as a secret
from base64 import b64encode, b64decode
from nacl import encoding, public
import requests
import os
import json
import boto3


def encrypt(raw_public_key: str, secret_value: str) -> str:
    """Encrypt a Unicode string using the public key."""
    public_key = public.PublicKey(raw_public_key.encode("utf-8"), encoding.Base64Encoder())
    sealed_box = public.SealedBox(public_key)
    encrypted = sealed_box.encrypt(secret_value.encode("utf-8"))
    return b64encode(encrypted).decode("utf-8")


def get_ecr_password() -> str:
    """Retrieve ECR password, it comes b64 encoded, in the format user:password
       From https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ecr.html#ECR.Client.get_authorization_token
    """
    ecr = boto3.client('ecr')
    response = ecr.get_authorization_token()
    encoded_login_password = response['authorizationData'][0]['authorizationToken']

    decoded_login_password = b64decode(encoded_login_password).decode('UTF-8')
    return decoded_login_password.split(':')[1]


if __name__ == '__main__':

    get_public_key = requests.get('https://api.github.com/repos/ORG/REPOSITORY/actions/secrets/public-key',
                                  headers={'Accept': 'application/vnd.github.v3+json',
                                           'Authorization': 'token ' + os.environ['GH_API_ACCESS_TOKEN']})
    if get_public_key.ok is False:
        print('could not retrieve public key')
        print(get_public_key.text)
        exit(1)
    get_public_key_response = get_public_key.json()
    public_key_value = get_public_key_response['key']
    public_key_id = get_public_key_response['key_id']

    password = get_ecr_password()
    encrypted_password = encrypt(public_key_value, password)
    update_password = requests.put('https://api.github.com/repos/ORG/REPOSITORY/actions/secrets/ECR_PASSWORD',
                                   headers={'Accept': 'application/vnd.github.v3+json',
                                            'Authorization': 'token ' + os.environ['GH_API_ACCESS_TOKEN']},
                                   data=json.dumps({'encrypted_value': encrypted_password, 'key_id': public_key_id,
                                                    'visibility': 'all'}))
    if update_password.ok is False:
        print('could not update password')
        print(update_password.text)
        exit(1)
