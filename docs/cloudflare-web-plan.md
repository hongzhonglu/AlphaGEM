# AlphaGEM Web + Cloudflare Tunnel Plan

## Goal

Expose a locally deployed AlphaGEM instance through Cloudflare Tunnel so users can:

- log in
- submit jobs
- track job status
- download results

The compute must continue to run on an Ubuntu host that already satisfies AlphaGEM's local requirements. Cloudflare is only the secure ingress layer.

## Current State

The repository is currently a CLI-style pipeline, not a web service.

- Main workflow entrypoint: `AlphaGEM.py`
- Outputs are written under `working/<name>`
- No API layer
- No authentication layer
- No job queue
- No persistent job metadata store
- No multi-user isolation

Relevant code:

- `AlphaGEM.py` parses CLI arguments and executes the full pipeline
- `README.md` documents Ubuntu-only deployment and local dependency requirements

## Feasibility

### Feasible

- Publish a local web app with Cloudflare Tunnel
- Protect the app with Cloudflare Access login
- Build a lightweight frontend for job submission and status tracking
- Use the current AlphaGEM CLI as the execution backend

### Not Feasible

- Running AlphaGEM inside Cloudflare Workers or on Tunnel infrastructure itself
- Allowing unrestricted concurrent jobs on one machine without resource controls
- Treating the current CLI entrypoint as a production web backend without a wrapper

## Recommended Architecture

### MVP Architecture

```text
Browser
  -> Cloudflare Access
  -> Cloudflare Tunnel
  -> FastAPI backend on localhost:8000
  -> Job runner on localhost
  -> AlphaGEM.py subprocess
  -> job-specific output directory
```

### Components

1. Frontend

- Next.js or React SPA
- Login handled by Cloudflare Access before traffic reaches the app
- Calls backend REST API

2. Backend

- FastAPI
- Validates Cloudflare Access identity
- Accepts uploads
- Stores job metadata
- Starts and monitors local jobs

3. Job execution

- Start AlphaGEM via `subprocess.Popen(...)`
- One job directory per task
- Persist logs and exit status

4. Storage

- SQLite for MVP
- Local disk for input, logs, outputs

5. Access control

- Cloudflare Access as the first auth gate
- Backend maps authenticated email to a local user record

## Why This Architecture Fits AlphaGEM

AlphaGEM has the following properties:

- long-running jobs
- local file inputs and outputs
- GPU or heavy CPU dependency
- Ubuntu-only runtime
- external toolchain dependency

This means:

- job submission must be asynchronous
- HTTP requests must not wait for job completion
- job state must be persisted
- each job must run in an isolated directory

## Proposed Repository Additions

```text
docs/
  cloudflare-web-plan.md

web-backend/
  app/
    main.py
    api/
    core/
    models/
    services/
    workers/
  requirements.txt

web-frontend/
  src/
  package.json

storage/
  jobs/
  uploads/

var/
  alphagem.db
```

## Backend Design

### Core Responsibilities

- verify Cloudflare Access identity
- create and list jobs
- persist job metadata
- store uploads
- launch jobs
- collect logs
- expose downloadable outputs

### Suggested FastAPI Endpoints

#### Authentication / identity

- `GET /api/me`
  - returns authenticated user info

#### Jobs

- `POST /api/jobs`
  - create a new job
  - accepts form fields and file uploads

- `GET /api/jobs`
  - list current user's jobs

- `GET /api/jobs/{job_id}`
  - get one job's metadata and status

- `GET /api/jobs/{job_id}/logs`
  - fetch recent logs

- `POST /api/jobs/{job_id}/cancel`
  - request cancellation

- `GET /api/jobs/{job_id}/artifacts`
  - list files available for download

- `GET /api/jobs/{job_id}/artifacts/{artifact_id}`
  - download one output file

#### Health

- `GET /healthz`
  - simple liveness check for Tunnel origin

### Job State Model

Use a finite state machine:

- `queued`
- `running`
- `succeeded`
- `failed`
- `canceled`

### SQLite Tables

#### `users`

- `id`
- `email`
- `display_name`
- `role`
- `created_at`
- `last_seen_at`

#### `jobs`

- `id`
- `user_id`
- `public_name`
- `mode`
- `refname`
- `growth_medium`
- `status`
- `created_at`
- `started_at`
- `finished_at`
- `working_dir`
- `exit_code`
- `error_message`
- `log_path`

#### `job_inputs`

- `id`
- `job_id`
- `kind`
- `stored_path`
- `original_name`

#### `job_artifacts`

- `id`
- `job_id`
- `kind`
- `stored_path`
- `display_name`
- `size_bytes`

## Job Directory Layout

Do not continue relying on shared `working/<name>` semantics for the web version.

Recommended layout:

```text
storage/jobs/<job_id>/
  input/
    genome.fasta
    maplist.xlsx
    structure.zip
  output/
    final/
    intermediate/
  logs/
    stdout.log
    stderr.log
    state.json
  meta/
    request.json
```

Internally, the backend can still prepare a working directory compatible with AlphaGEM.

## Wrapping the Existing AlphaGEM CLI

### Short-Term Approach

Keep `AlphaGEM.py` intact and execute it as a subprocess.

Example shape:

```python
cmd = [
    "python",
    "AlphaGEM.py",
    "--name", internal_job_name,
    "--mode", mode,
    "--refname", refname,
    "--fasta", fasta_path,
]
```

Add optional arguments only if present and validated.

### Medium-Term Refactor

Refactor `AlphaGEM.py` into:

- argument parser
- pure execution function such as `run_alphagem(config)`
- CLI wrapper calling the same function

This enables:

- cleaner worker integration
- easier testing
- better error handling
- future batch scheduling

## Frontend Design

### MVP Pages

1. Login landing page

- explains what AlphaGEM does
- sign-in is handled by Cloudflare Access

2. Dashboard

- shows recent jobs
- shows current user
- button to create a new job

3. New job page

- form for required parameters
- upload FASTA
- optional advanced settings
- optional structure inputs

4. Job detail page

- job state
- timestamps
- recent logs
- download buttons
- error message if failed

### Suggested Form Fields

Required:

- `public_name`
- `mode`
- `refname`
- `fasta`

Optional:

- `maplist`
- `structure archive`
- `TMscore`
- `coverage`
- `pLDDT`
- `grothmedium`

### UX Constraints

- never block on job completion during submission
- validate file size and file type early
- show that jobs may take a long time
- clearly separate required and advanced parameters

## Authentication Strategy

### Recommended

Use Cloudflare Access for sign-in, then validate identity in the backend.

Backend trust chain:

1. Cloudflare Access authenticates user
2. Access forwards identity material to origin
3. FastAPI verifies the Access JWT
4. FastAPI resolves or creates a local user record

### Why Not Build Password Login First

- more code
- more security surface
- password reset and email verification overhead
- unnecessary if the app is already behind Cloudflare Access

## Cloudflare Tunnel Setup

### Origin Assumption

Backend is listening on:

```text
http://localhost:8000
```

### Example `cloudflared` config

```yaml
tunnel: <tunnel-id>
credentials-file: /root/.cloudflared/<tunnel-id>.json

ingress:
  - hostname: alphagem.example.com
    service: http://localhost:8000
  - service: http_status:404
```

### Cloudflare Steps

1. Install `cloudflared` on the AlphaGEM host
2. Authenticate `cloudflared` with your Cloudflare account
3. Create a named tunnel
4. Add DNS route for `alphagem.example.com`
5. Start `cloudflared` as a service
6. In Zero Trust, create a self-hosted application
7. Add Access policies for allowed users or domains

## Security Controls

### Minimum Controls

- Cloudflare Access in front of the app
- JWT validation in backend
- user-to-job ownership checks
- file type validation
- file size limits
- sanitized job names
- no shell interpolation
- download authorization checks

### Strongly Recommended

- admin-only access at first
- single concurrent job until resource behavior is understood
- per-user quota
- upload retention rules
- periodic cleanup for finished jobs

## Concurrency and Scheduling

### MVP

- allow only one running job at a time
- queue additional jobs
- show queue position in UI if desired

### Later

- use Redis + Celery or RQ
- support worker isolation
- support per-machine capacity limits

## Operational Risks

### Resource contention

AlphaGEM may consume substantial CPU, RAM, GPU, disk, and temporary storage. Multiple concurrent jobs can make the host unstable.

### Disk growth

The repo already contains large data assets and many output directories. Web usage will increase storage pressure quickly.

### Long runtimes

Users need explicit expectation setting. Jobs may take minutes to hours depending on input and mode.

### Input complexity

`structure_alignment` requires more complex inputs than `plmsearch`. The UI must reflect that.

### User isolation

The current `working/<name>` layout is not safe enough as a long-term multi-user contract.

## Recommended Delivery Phases

### Phase 0: Preparation

- identify host machine that will run the service
- confirm AlphaGEM runs end-to-end there
- document expected runtime and peak resource use

### Phase 1: API Wrapper

- create FastAPI app
- add SQLite models
- add job submission endpoint
- add background runner
- persist logs and job states

### Phase 2: Web UI

- build dashboard
- build new job form
- build job detail page
- build artifact download page

### Phase 3: Cloudflare

- set up Tunnel
- set up Access
- validate JWT in backend
- test remote access

### Phase 4: Hardening

- add quotas
- add cleanup tasks
- improve cancellation
- improve error surfacing
- add monitoring

## First Implementation Slice

The best first slice is:

1. Create backend skeleton
2. Add `POST /api/jobs`
3. Add `GET /api/jobs/{job_id}`
4. Run AlphaGEM as a subprocess in single-job mode
5. Store logs and outputs
6. Expose basic HTML or JSON only
7. Add frontend after the backend contract stabilizes

This is better than starting with frontend because the difficult part is execution control, not page layout.

## Success Criteria For MVP

- a logged-in user can submit one valid FASTA job
- the backend stores metadata and files safely
- the job runs without blocking the HTTP request
- the user can refresh and see status changes
- the user can download at least one final result artifact
- the app is reachable only through Cloudflare Access

## Recommendation

Start with:

- FastAPI backend
- SQLite
- single-host single-worker queue
- Cloudflare Access
- minimal React or Next.js frontend

Do not start with:

- custom password auth
- multi-worker scaling
- object storage migration
- complicated frontend polish

## Useful Official References

- Cloudflare Tunnel configuration file:
  - https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/local-management/configuration-file/
- Cloudflare Access for self-hosted applications:
  - https://developers.cloudflare.com/cloudflare-one/access-controls/applications/http-apps/self-hosted-public-app/
- Validating Cloudflare Access JWTs at the origin:
  - https://developers.cloudflare.com/cloudflare-one/access-controls/applications/http-apps/authorization-cookie/validating-json/
