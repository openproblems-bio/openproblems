set -x
pip install -U ..[evaluate]
snakemake -j $(grep -c processor /proc/cpuinfo) all || travis_terminate 1
if [ "$TRAVIS_PYTHON_VERSION" = "3.7" ]; then
    git remote set-url origin https://singlecellopenproblems:${GITHUB_PASSWORD}@github.com/${TRAVIS_REPO_SLUG}.git
    git checkout -B master
    git pull origin master
    git add ../results.md
    git add ../website/data/results
    git config user.name "Travis CI"
    git config user.email "singlecellopenproblems@protonmail.com"
    git commit -m "Travis CI: update results for ${TRAVIS_TAG}"
    git push origin master
fi
