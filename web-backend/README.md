# AlphaGEM Web Backend

## Scope

This is the first backend skeleton for exposing AlphaGEM as a web service.

Current capabilities:

- FastAPI app
- static frontend served from `/`
- local registration and login
- seeded admin account
- SQLite metadata store
- local single-thread job queue
- job submission by file upload
- job status polling
- log retrieval
- artifact download

## Run

From the repository root:

```bash
cd web-backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

Open:

- `http://localhost:8000/`
- `http://localhost:8000/healthz`
- `http://localhost:8000/docs`

## Authentication

This backend now supports:

- local registration and login using bearer tokens
- `Cf-Access-Authenticated-User-Email` for future Cloudflare Access integration
- optional `X-Dev-User-Email` for local development fallback only if you explicitly use it

Seeded admin account by default:

- email: `admin@alphagem.local`
- password: `AlphaGEM_admin_123`

Override with:

- `ALPHAGEM_ADMIN_EMAIL`
- `ALPHAGEM_ADMIN_PASSWORD`
- `ALPHAGEM_AUTH_SECRET`

By default, anonymous access is disabled. If you set `ALPHAGEM_DEFAULT_USER_EMAIL`, the backend will allow a fallback local identity; do not enable that in user-facing deployments.

## Example Job Submission

`plmsearch` mode:

```bash
TOKEN="<login returned access_token>"
curl -X POST http://localhost:8000/api/jobs \
  -H "Authorization: Bearer ${TOKEN}" \
  -F "public_name=test-job" \
  -F "mode=plmsearch" \
  -F "refname=yeast" \
  -F "grothmedium=min" \
  -F "fasta=@/absolute/path/to/input.fasta"
```

## Important Limitations

- only one worker thread processes jobs
- cancellation is not implemented yet
- queue position is approximate
- `structure_alignment` requires both `maplist` and a supported archive for structures
- Cloudflare Access JWT validation is not implemented yet; the backend currently trusts forwarded Cloudflare identity headers

## Storage

The backend creates:

- `storage/jobs/`
- `var/alphagem.db`

Job logs and uploaded files are stored under `storage/jobs/<internal_name>/`.
