# ketcher-shell

One-shot Vite shell that bundles `ketcher-react` + `ketcher-standalone` for
ChemAgent Console embedding (iframe under `/vendor/ketcher/`).

## Requirements

- Node 18+
- npm

## Build vendor into `web/vendor/ketcher/`

From repo root:

```bash
bash tools/build_ketcher_vendor.sh
```

Pinned version: **3.7.0** (see `web/vendor/ketcher/VERSION` after build).

Do not commit `node_modules/` or this package's `dist/` — only the synced
`web/vendor/ketcher/` output is committed for offline Console use.
