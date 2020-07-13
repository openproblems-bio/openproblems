if [ "$TRAVIS_PYTHON_VERSION" = "3.6" ]; then
    set -x
    pip install -U ..[evaluate]
    python evaluate.py
    git checkout master
    git set-url origin https://scottgigante:${GITHUB_PASSWORD}@github.com/${TRAVIS_REPO_SLUG}.git
    git pull origin master
    git add ../results.md
    git commit -m "Travis CI: update results for ${TRAVIS_TAG}"
    git push origin master
fi
