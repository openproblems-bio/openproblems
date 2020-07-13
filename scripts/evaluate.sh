set -x
pip install -U ..[evaluate]
python evaluate.py
if [ "$TRAVIS_PYTHON_VERSION" = "3.6" ]; then
    git remote set-url origin https://scottgigante:${GITHUB_PASSWORD}@github.com/${TRAVIS_REPO_SLUG}.git
    git fetch --depth=1 origin master
    git checkout master
    git add ../results.md
    git config user.name "Travis CI"
    git config user.email "scottgigante@gmail.com"
    git commit -m "Travis CI: update results for ${TRAVIS_TAG}"
    git push origin master
fi
