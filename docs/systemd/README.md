# systemd examples

These units are examples, not drop-in production files.

Adjust at least:

- `User` and `Group`
- `WorkingDirectory`
- `ExecStart`
- project root paths

Recommended split:

- `rnaa-download.service`: long-running independent download worker
- `rnaa-run.service` plus `rnaa-run.timer`: periodic orchestration without owning downloads
