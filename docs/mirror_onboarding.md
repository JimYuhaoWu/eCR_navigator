# Mirror onboarding playbook (reusable for every model)

Each genomic foundation model (ChromBERT, ChromFound, GET, EpiAgent, …) runs in
its **own server mirror** with incompatible deps/GPU runtime. Every mirror is
onboarded the same way — only the per-model provisioning at the end differs. Follow
this once per new mirror; it takes ~5 min plus any model download.

> Platform note: the mirrors run on a **university-built managed GPU platform**.
> Its defining quirk: the **data disk persists** across restarts (`~/.cache`,
> `/root/bin`, `/root/.bashrc`), but **`/root/.ssh/authorized_keys` is regenerated
> on every boot** — so the SSH key never survives a restart on its own. The
> `.bashrc` self-heal (step 4) is what makes access durable.

## What the user provides

- **IP** (fixed per mirror) and **PORT** (may change on restart), e.g.
  `root@172.16.78.10 -p 35963`.
- One-time **password** (or web-terminal access) to install the key the first time.

## Steps

### 1. Reuse the automation keypair
One keypair is shared across all mirrors: `~/.ssh/ecr_navigator` (ed25519).
Public key:
```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDgrjZA4n+pYrj+RIQlCmDuOofS+wfaUU+sAHec6nI1b ecr_navigator@claude
```
If it does not exist yet: `ssh-keygen -t ed25519 -f ~/.ssh/ecr_navigator -N "" -C ecr_navigator@claude`.

### 2. Bootstrap access (once, via the platform web terminal)
The user pastes the public key into `authorized_keys` on the mirror:
```bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
echo 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDgrjZA4n+pYrj+RIQlCmDuOofS+wfaUU+sAHec6nI1b ecr_navigator@claude' >> ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
```
Do NOT paste the password into chat — the user runs this in their own terminal.

### 3. Point the tooling at the mirror
`scripts/mirror_env.sh` holds the connection (IP fixed, port overridable):
```bash
MIRROR_PORT=<port> ./scripts/setup_mirror.sh <genome>
```
If a new mirror has a different IP, set `MIRROR_IP` too. Verify:
```bash
ssh -i ~/.ssh/ecr_navigator -p <port> root@<ip> 'echo OK'
```

### 4. Install the `.bashrc` self-heal (makes access survive restarts)
`scripts/setup_mirror.sh` adds a silent, idempotent block to `/root/.bashrc` that
re-appends the key whenever a shell starts. Because `.bashrc` can't authenticate
the *first* SSH login, the platform **web terminal** (which logs in without SSH and
sources `.bashrc`) is what triggers it.

**Per-restart procedure from then on: open the web terminal once → key restored →
automated SSH works.** (Fully automatic if the platform runs a login shell at boot,
or if it ever exposes a dashboard SSH-key / startup-script field — prefer that.)

Then **save/snapshot the mirror** so `.bashrc`, `/root/bin`, and downloaded model
data persist. Snapshot rarely — see the snapshot policy in `server_mirrors.md`.

### 5. Provision the model runtime (the only per-model part)
`setup_mirror.sh` also handles the ChromBERT-generic bits (bedtools, model
download). For a NEW model, add a sibling provisioning path — its env, deps, and
data — following the same "keep scripts in the repo, run on the mirror, keep the
navigator env light, hand off via the `.npz` embedding-artifact contract" pattern
(`docs/embedding_artifact.md`).

## Quick reference

| Fact | Value |
|---|---|
| Key | `~/.ssh/ecr_navigator` (shared across mirrors) |
| IP | fixed per mirror; PORT may change → `MIRROR_PORT=<port>` |
| authorized_keys | wiped every boot → self-healed by `.bashrc` |
| Re-connect after restart | open web terminal once |
| `.bashrc` backup | `~/.bashrc.ecr_bak` |
| Per-model detail | `server_mirrors.md` (one section per model) |
