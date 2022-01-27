# General Flow/Commands for Version Control

## Checking branches
```
git status (shows your current branch)
git branch (shows local branches)
git branch -r (shows remote branches)
git branch -a (shows all branches)
```

## Creating + Switching to a branch

```
git branch <branch name> (Creates the branch)
git checkout <branch name> (Switches to the branch)
```

## Pushing the branch to the git repository
While on the branch:
```
git add <changes>
git commit -m "<commit message>"
git push origin <branch name>
```

## Making a Pull Request
Go to Github -> Pull Requests -> New Pull Request -> Base: Main, Compare: Your Branch

## Deleting a branch
```
git branch -d <branch name> (Delete branch locally)
git push origin --delete <branch name> (Delete branch remotely)
```


