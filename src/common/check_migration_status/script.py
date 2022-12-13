import json

## VIASH START

par = {
    'git_sha': 'temp/openproblems-v1.json',
    'comp_info': 'temp/denoising_metrics.yaml',
    'output': 'temp/migration_status.yaml'
}

## VIASH END

output = {}

with open(par['git_sha'], 'r') as f1:
    git = json.load(f1)

    with open(par['comp_info'], 'r') as f2:
        comp = json.load(f2)
        for comp_item in comp:
            if comp_item['namespace'] not in output:
                output[comp_item['namespace']] = {}
            if comp_item['v1_url']:
                for obj in git:
                    if obj['path'] in comp_item['v1_url']:
                        if obj['sha'] != comp_item['v1_commit']:
                            output[comp_item['namespace']][comp_item['id']] = {'v1_url': comp_item['v1_url'], 'status': "not latest commit"}
                        else :
                            output[comp_item['namespace']][comp_item['id']] = {'v1_url': comp_item['v1_url'],'status': "up to date"}
            else:
                output[comp_item['namespace']][comp_item['id']] ={'v1_url': "v1_url missing"}


with open(par['output'], 'w') as outf:
    json.dump(output, outf)






