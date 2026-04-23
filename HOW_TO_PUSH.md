# How To Push This Project

This folder is meant to be managed as its own Git repository.

## Important rule

From now on, do Git work **inside this folder**, not in the larger `MACLEAN_LAB` root:

```bash
cd /Users/moyingli/USC/SP26/MACLEAN_LAB/github_ready_maclean_lab
```

## First push

Run:

```bash
git push -u origin main
```

If GitHub asks for credentials:

- Username: `mungbean138516-jpg`
- Password: use a GitHub personal access token, not your normal GitHub password

## Normal update workflow

Whenever you change files and want GitHub to reflect those changes:

```bash
cd /Users/moyingli/USC/SP26/MACLEAN_LAB/github_ready_maclean_lab
git status
git add .
git commit -m "Describe what changed"
git push
```

## What "push" means

- `git add`: stage the files you want to include
- `git commit`: save a local snapshot
- `git push`: upload your local commits to GitHub

You do **not** need to push every time you save a file.
Only push when you want the remote GitHub repository updated.

## Recommended habit

- Small edits: commit after a meaningful chunk of work
- Good commit messages:
  - `Refine beta6 sweep summary`
  - `Add post-meeting robustness plots`
  - `Clean README and reorganize results`

## If you want to see what changed

```bash
git status
git diff
```

## If Git says nothing to commit

That means your working tree is already clean, and there is nothing new to push.
