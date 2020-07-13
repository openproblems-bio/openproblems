if [ "$TRAVIS_PYTHON_VERSION" = "3.6" ]; then
    python evaluate.py
    git add ../results.md
    git commit -m "update results for ${TRAVIS_TAG}"
    git push https://scottgigante:${GITHUB_PASSWORD}@github.com/${TRAVIS_REPO_SLUG}.git master
fi
