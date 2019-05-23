In order to create a tagged version of the software, perform the following:

1. Push all uncommitted changes, so that the master branch is up to date
2. Update version number in __init__, meta.yaml, and docs
3. ```git add .```
4. ```git tag``` to see last tag
5. ```git tag <tag_id```
6. ```git commit -m "message"```
7. ```git push origin <tag_id>```
8. Go to GitHub versions, click the tag, edit, add release id to title and add notes
