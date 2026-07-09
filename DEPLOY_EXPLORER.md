# Deploying SAAP Database on Northeastern Explorer (lab-wide)

Goal: one shared instance the whole lab reaches, where **anyone can browse/export**
but **importing, deleting, and editing DOIs require a password**.

> ⚠️ **Confirm with Research Computing (RC) first.** HPC clusters generally do **not**
> allow long-running network services on **login nodes**, and open ports are usually
> firewalled. The setup below works today via an SSH tunnel (per user). For a single
> stable URL the whole lab opens without tunneling, you'll need RC to provide a
> reverse proxy / Open OnDemand app or a dedicated always-on host — open a ticket for that.
> Everything here also assumes lab members are on the **Northeastern network or VPN**.

---

## 1. One-time setup

SSH in and put the code on a **shared filesystem** your lab can all reach (e.g. a
`/work/<yourlab>/…` project space), so both the app and its database are central:

```bash
ssh <yourusername>@login.explorer.northeastern.edu

# choose a shared project location
cd /work/<yourlab>            # adjust to your actual project path
git clone <your-repo-url> saap        # or copy the SAAP_Database folder here
cd saap

module load python/3.12       # or: module load anaconda3 ; whichever RC provides
python3 -m venv backend/.venv
source backend/.venv/bin/activate
pip install -r backend/requirements.txt
```

## 2. Choose the shared database location and password

The database is a single SQLite file. Point every launch at the **same shared path**
so all lab members read/write one database:

```bash
export SAAP_DB_PATH=/work/<yourlab>/saap/data/saap.db
export SAAP_WRITE_PASSWORD='pick-a-strong-shared-password'
```

The parent directory is created automatically. Anyone with write access to that
directory (and the password) can modify data; grant filesystem permissions to your
lab group accordingly (e.g. `chmod -R g+rw /work/<yourlab>/saap`).

## 3. Run the server

`serve.sh` binds to the network and enables the password:

```bash
SAAP_DB_PATH=/work/<yourlab>/saap/data/saap.db \
SAAP_WRITE_PASSWORD='pick-a-strong-shared-password' \
PORT=8000 ./serve.sh
```

Keep it running across your SSH session with **tmux** (so it survives disconnects):

```bash
tmux new -s saap
# … run serve.sh inside …
# detach with Ctrl-b then d ; reattach later with:  tmux attach -t saap
```

> If RC prohibits services on login nodes, run it inside a **Slurm allocation**
> instead (a long-walltime interactive or batch job on a compute node) and note the
> node name (e.g. `d1001`) for the tunnel below. Ask RC which pattern they allow.

## 4. How the lab reaches it

**Option A — SSH tunnel (works today, each user does this once):**

From a lab member's own computer:

```bash
# if serve.sh runs on the login node:
ssh -L 8000:localhost:8000 <username>@login.explorer.northeastern.edu

# if it runs on a compute node named e.g. d1001:
ssh -L 8000:d1001:8000 <username>@login.explorer.northeastern.edu
```

Then open **http://localhost:8000** in their browser. Browsing/exporting works with
no password; clicking **“Read-only · unlock to edit”** in the header and entering the
shared password enables import/delete.

**Option B — one shared URL (no tunneling):** requires RC to expose the port through a
reverse proxy or Open OnDemand. Open a ticket describing a small internal FastAPI web
app on port 8000 that should be reachable on the campus network. Then everyone just
opens the URL RC gives you.

---

## Backups & maintenance

- **Back up** the single file `$SAAP_DB_PATH` (plus `-wal`/`-shm` siblings if present).
  A cron `cp` to a dated filename is enough:
  `cp /work/<yourlab>/saap/data/saap.db backups/saap-$(date +%F).db`.
- **Update the app:** `git pull` (or recopy), then restart `serve.sh`. Schema changes
  that only *add* columns migrate automatically on startup without data loss.
- **Reset:** stop the server and delete the database file, or use **Datasets → Clear all
  data** in the UI (dataset DOIs are preserved).

## Security notes

- The write password is sent in a request header. On a tunneled/VPN connection this is
  fine; if you ever expose a public URL, insist on **HTTPS** (RC's proxy should provide
  it) so the password isn't sent in clear text.
- `serve.sh` warns if `SAAP_WRITE_PASSWORD` is unset — never run the shared instance
  without it, or anyone who can reach the port could wipe the database.
- Read access is open by design (your choice). If you later want to gate reading too,
  ask and I'll add a second, read password.
