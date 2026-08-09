# AGENTS.md - Agent Workflow Instructions

This project uses **beads** (`bd`) for git-backed issue tracking. See https://github.com/steveyegge/beads

## Essential Commands

| Command | Purpose |
|---------|---------|
| `bd ready` | List tasks without blockers (your next work) |
| `bd create "title" -p 1` | Create task (P0=critical, P1=high, P2=medium, P3=low) |
| `bd show <id>` | View issue details and history |
| `bd update <id> --status in_progress` | Mark task as in progress |
| `bd close <id> --reason "text"` | Close completed task |
| `bd dep add <child> <parent>` | Add dependency |
| `bd list --json` | List all open issues |
| `bd export --no-memories -o .beads/issues.jsonl` | Export embedded tracker state to Git |

## Critical Rules for Agents

1. **NEVER use `bd edit`** - it opens an interactive editor. Use flag-based updates:
   ```bash
   bd update <id> --description "new description"
   bd update <id> --title "new title"
   ```

2. **Always use `--json` flag** for programmatic access

3. **Export after tracker changes** with
   `bd export --no-memories -o .beads/issues.jsonl`. The installed embedded-Dolt
   `bd` no longer provides `bd sync`.

## Landing the Plane Protocol

When ending a work session, you MUST complete these steps in order:

1. **File remaining work** as new issues for anything not completed
2. **Run quality gates** (tests, linting, builds as appropriate)
3. **Update issue statuses** - close completed, update in-progress
4. **Sync and push**:
   ```bash
   bd export --no-memories -o .beads/issues.jsonl
   git add .beads/issues.jsonl .beads/interactions.jsonl
   git pull --rebase
   git push
   ```
5. **Verify clean state**: `git status` shows nothing pending
6. **Provide handoff context** for next session

**Work is NOT complete until `git push` succeeds.**

## Finding Work

```bash
bd ready --json          # Tasks without blockers
bd list --status open    # All open tasks
bd stale --days 7        # Neglected tasks
```

## Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --status in_progress  # Claim work
bd close <id>         # Complete work
bd export --no-memories -o .beads/issues.jsonl  # Export tracker state
```

## Critical Rules

- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
