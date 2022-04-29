# How to contribute
Any changes in this GitHub repository should be made by [creating a new branch](#creating-a-new-branch), committing and pushing changes, [creating a pull request](#creating-a-pull-request), and [merging the pull request](#merging-a-pull-request). [Releases](#creating-a-release) should be generated periodically along with the associated Docker containers.

## Creating a new branch
Before making any changes, [create a new branch](https://www.howtogeek.com/714112/how-to-create-a-new-branch-in-github/). To make a new branch, enter the following on terminal,
```bash
git checkout -b {new-branch-name}
```

Replace `{new-branch-name}` with something short that is self-descriptive. All changes should be made in this GitHub repository. Changes should be [committed](https://github.com/git-guides/git-commit) and [pushed](https://github.com/git-guides/git-push) to the repository frequently.

## Creating a pull request 
Navigate to the [branches](https://github.com/vanallenlab/moalmanac/branches) page of this GitHub repository and click `New pull request`, found on the right side of the `Your branches` section. Revise the pull request title to be self-descriptive and copy the [pull request template](template-pull-request.md) to the comment body. Fill out the pull request template and click `Create pull request`. Comment on the pull request as the pull request is reviewed and the pre-pull request check list is completed.

## Merging a pull request
Upon completing review of a pull request and all items on the pre-pull request checklist are checked, click `Squash and merge` to merge the pull request. The branch can then be [deleted](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-branches-in-your-repository/deleting-and-restoring-branches-in-a-pull-request).

## Creating a release
Releases follow [semantic versioning](https://semver.org/) and should be created, at minimum, for every minor patch. Navigate to the [releases](https://github.com/vanallenlab/moalmanac/releases) page and click `Draft new release`. Enter the major, minor, and patch numbers delimited by periods as the tag number, such as `0.4.2`. The release title should be `Release {major}.{minor}.{patch}`, such as `Release 0.4.2`. For the release description, follow the [release template](template-release.md) to describe the merged pull requests since the previous release.
