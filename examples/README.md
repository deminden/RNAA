# Offline installation check

From the workspace root, with Python 3.11+ and a built or installed executable:

```bash
python3 examples/offline_reference.py --binary target/release/rnaa
python3 examples/offline_reference.py --binary target/release/rnaa --automatic-resources
# Or pass /your/install/prefix/bin/rnaa after installing from source.
```

The example uses a temporary project and a private resource pool. It initializes
the catalog, prepares two synthetic transcripts from local FASTA/GTF, verifies
cache reuse, and checks status and engine diagnostics. Temporary data is removed
on exit. No network access, R installation or existing project is needed.
The automatic-resources variant also exercises operating-system CPU/memory
detection while retaining a private pool, including the Linux detection path.

This is an installation/reference smoke check. It does not validate a public
cohort, quantify reads or establish differential-expression accuracy; the native
workflow integrations and independent fixtures cover those scientific paths.

For a small end-to-end public cohort, see the [yeast workflow](public_yeast.md).
