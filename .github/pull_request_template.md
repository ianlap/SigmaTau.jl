## Summary

Describe the change and the reason for it.

## Checklist

- [ ] Tests pass locally (`julia --project=. -e 'using Pkg; Pkg.test()'`).
- [ ] Added or updated tests covering the change.
- [ ] Exported functions have docstrings.
- [ ] `CHANGELOG.md` has an entry under `## [Unreleased]` (or the body notes
      why none is warranted, e.g. pure docs / typo).
- [ ] Removed the matching item from `TODO.md` if one existed.
- [ ] Refreshed `project_overview.md` if the public surface changed.
- [ ] For a core-kernel change: the Stable32 / allantools / legacy-parity
      testsets in `test/stab/` still pass.

## Notes

Call out boundary policies, breaking changes, reference citations, or anything
that needs extra review.
