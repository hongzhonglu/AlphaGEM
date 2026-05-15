from __future__ import annotations

import queue
import threading
from pathlib import Path

from app.core.config import settings
from app.core.db import utc_now
from app.services.alphagem_runner import collect_artifacts, copy_result_snapshot, resolve_output_dir, run_subprocess
from app.services.repository import add_artifact, get_job, update_job


class JobRunner:
    def __init__(self) -> None:
        self._queue: queue.Queue[str] = queue.Queue()
        self._thread: threading.Thread | None = None
        self._stop_event = threading.Event()

    def start(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._thread = threading.Thread(target=self._run_loop, name="alphagem-job-runner", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        self._queue.put("__stop__")
        if self._thread:
            self._thread.join(timeout=2)

    def enqueue(self, job_id: str) -> None:
        self._queue.put(job_id)

    def _run_loop(self) -> None:
        while not self._stop_event.is_set():
            job_id = self._queue.get()
            if job_id == "__stop__":
                break
            self._run_job(job_id)
            self._queue.task_done()

    def _run_job(self, job_id: str) -> None:
        job = get_job(job_id)
        if not job:
            return

        log_path = Path(job["log_path"])
        update_job(job_id, status="running", started_at=utc_now(), queue_position=None, error_message=None)
        try:
            exit_code = run_subprocess(job, log_path)
            output_dir = resolve_output_dir(job)
            snapshot_dir = Path(job["working_dir"]) / "output"
            copy_result_snapshot(output_dir, snapshot_dir)

            if exit_code == 0:
                update_job(
                    job_id,
                    status="succeeded",
                    finished_at=utc_now(),
                    exit_code=exit_code,
                    output_dir=str(output_dir),
                )
            else:
                update_job(
                    job_id,
                    status="failed",
                    finished_at=utc_now(),
                    exit_code=exit_code,
                    output_dir=str(output_dir),
                    error_message=f"AlphaGEM exited with code {exit_code}",
                )

            for artifact_path in collect_artifacts(snapshot_dir):
                add_artifact(
                    job_id=job_id,
                    kind="output",
                    stored_path=artifact_path,
                    display_name=artifact_path.name,
                    size_bytes=artifact_path.stat().st_size,
                )
        except Exception as exc:  # pragma: no cover - safety path
            update_job(
                job_id,
                status="failed",
                finished_at=utc_now(),
                error_message=str(exc),
            )
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write(f"\n[runner-error] {exc}\n")


job_runner = JobRunner()
