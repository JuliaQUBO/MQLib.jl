# Dependabot automatic merges

All same-repository Dependabot-authored PRs to `main`, including major and grouped
updates, are eligible when their current head passes the checks below and every
additional reported check has finished successfully (neutral/skipped optional
checks are allowed). Required checks cannot be skipped. Human PRs, drafts,
forks, changed heads, conflicts, and behind branches are left to maintainers.
Identity follows the PR author: maintainer commits on a Dependabot branch are
eligible only after that head passes the same gates.

The workflow checks out trusted `main`, never PR code or artifacts, and uses the
built-in Actions token. It re-reads identity and clean mergeability, verifies
enforced required checks and current-main ancestry, and squash-merges with `--match-head-commit`, without
admin bypass or a deferred merge request. Existing review/conversation gates
still apply. Errors are isolated per PR and reported after both merge and
publication passes; one failure cannot starve the other PRs.

## Required successful checks

- `Native C ABI smoke test`
- `Julia 1.10 - ubuntu-latest - x64 - pull_request`
- `Julia 1 - ubuntu-latest - x64 - pull_request`
- `Julia 1 - windows-latest - x64 - pull_request`
- `Dependabot policy`

## Branch protection and activation

Keep the existing strict, up-to-date branch protection and required CI checks.
The workflow verifies the public branch summary: protection is enabled and these
contexts are enforced for non-admins or everyone:

- `Julia 1 - ubuntu-latest - x64 - pull_request`
- `Julia 1.10 - ubuntu-latest - x64 - pull_request`
- `Julia 1 - windows-latest - x64 - pull_request`

The public branch summary does not expose the admin-only `strict` setting. The
workflow separately compares the current main SHA with the verified PR head and
requires main to be its ancestor. GitHub's strict rule also checks up-to-date
status atomically when merging; the workflow never bypasses it.

When `main` advances, update
existing Dependabot branches with "Update branch" or `@dependabot rebase` and
wait for new CI. The policy runs after configured workflow completions, every
15 minutes, and through manual dispatch. Repository native auto-merge settings
only control manually queued requests; this workflow merges immediately after
verifying all gates.

## Publication and live verification

`GITHUB_TOKEN` merges suppress ordinary push/PR-close workflows. This repository has no documentation publisher to dispatch after token merges.
Main push CI is also suppressed. The up-to-date PR CI remains required.
Package releases remain manual.

After activation, verify the first real token merge at its merge SHA.
GitHub documents Contents write access
for the [PR merge API](https://docs.github.com/en/rest/pulls/pulls#merge-a-pull-request).
Actions-file updates use the same path; token permissions and completion events
still need live verification. If GitHub rejects such a merge, a maintainer merges
that PR manually; its visible error does not block other PRs. No PAT or new secret
is required. Do not force-push publication history or hide failures to obtain green CI.
