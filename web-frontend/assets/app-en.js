const state = {
  autoRefresh: true,
  timerId: null,
  jobs: [],
  authMode: "login",
  token: localStorage.getItem("alphagem_token_en") || "",
  user: readStoredUser(),
};

const elements = {
  authForm: document.getElementById("auth-form"),
  authEmail: document.getElementById("auth-email"),
  authPassword: document.getElementById("auth-password"),
  authDisplayName: document.getElementById("auth-display-name"),
  authSubmit: document.getElementById("auth-submit"),
  logout: document.getElementById("logout"),
  tabLogin: document.getElementById("tab-login"),
  tabRegister: document.getElementById("tab-register"),
  meOutput: document.getElementById("me-output"),
  authHint: document.getElementById("auth-hint"),
  refreshMe: document.getElementById("refresh-me"),
  refreshJobs: document.getElementById("refresh-jobs"),
  toggleAutoRefresh: document.getElementById("toggle-auto-refresh"),
  jobForm: document.getElementById("job-form"),
  submitStatus: document.getElementById("submit-status"),
  jobsList: document.getElementById("jobs-list"),
  jobsEmpty: document.getElementById("jobs-empty"),
  template: document.getElementById("job-card-template"),
  modeSelect: document.getElementById("mode-select"),
};

function readStoredUser() {
  try {
    return JSON.parse(localStorage.getItem("alphagem_user_en") || "null");
  } catch {
    return null;
  }
}

function writeStoredUser(user, token) {
  if (user) {
    localStorage.setItem("alphagem_user_en", JSON.stringify(user));
  } else {
    localStorage.removeItem("alphagem_user_en");
  }
  if (token) {
    localStorage.setItem("alphagem_token_en", token);
  } else {
    localStorage.removeItem("alphagem_token_en");
  }
}

function buildHeaders() {
  const headers = {};
  if (state.token) {
    headers.Authorization = `Bearer ${state.token}`;
  }
  return headers;
}

function setAuthMode(mode) {
  state.authMode = mode;
  elements.tabLogin.classList.toggle("active", mode === "login");
  elements.tabRegister.classList.toggle("active", mode === "register");
  elements.authDisplayName.style.display = mode === "register" ? "block" : "none";
  elements.authSubmit.textContent = mode === "register" ? "Create account" : "Login";
  elements.authHint.textContent =
    mode === "register"
      ? "Create an account to start using AlphaGEM."
      : "Log in to view jobs, logs, and result files.";
}

function setStructureFieldsVisibility() {
  const structureOnly = document.querySelectorAll("[data-structure-only='true']");
  const visible = elements.modeSelect.value === "structure_alignment";
  for (const block of structureOnly) {
    block.style.display = visible ? "block" : "none";
    const input = block.querySelector("input");
    if (input) {
      input.required = visible;
    }
  }
}

async function fetchJSON(url, options = {}) {
  const response = await fetch(url, options);
  if (!response.ok) {
    const contentType = response.headers.get("content-type") || "";
    if (contentType.includes("application/json")) {
      const data = await response.json();
      throw new Error(data.detail || JSON.stringify(data));
    }
    const body = await response.text();
    throw new Error(body || `Request failed: ${response.status}`);
  }
  return response.json();
}

function describeUser(user) {
  return `Logged in as: ${user.display_name || user.email}`;
}

async function refreshIdentity() {
  if (!state.token) {
    state.user = null;
    elements.meOutput.textContent = "Please log in or create an account.";
    return;
  }

  try {
    const user = await fetchJSON("/api/me", { headers: buildHeaders() });
    state.user = user;
    elements.meOutput.textContent = describeUser(user);
    writeStoredUser(user, state.token);
  } catch (error) {
    state.token = "";
    state.user = null;
    writeStoredUser(null, "");
    elements.meOutput.textContent = `Authentication failed: ${error.message}`;
  }
}

function renderJobs() {
  elements.jobsList.innerHTML = "";
  elements.jobsEmpty.style.display = state.jobs.length ? "none" : "block";

  if (!state.jobs.length) {
    elements.jobsEmpty.textContent = "No jobs yet. Log in and upload a FASTA file to create a new AlphaGEM run.";
  }

  for (const job of state.jobs) {
    const node = elements.template.content.firstElementChild.cloneNode(true);
    node.querySelector(".job-name").textContent = job.public_name;
    node.querySelector(".job-meta").textContent = `${job.mode} · ${job.refname} · ${job.growthMedium}`;
    node.querySelector(".job-id").textContent = job.id.slice(0, 8);
    node.querySelector(".job-created").textContent = new Date(job.created_at).toLocaleString();
    node.querySelector(".job-message").textContent = translateStatus(job.status, job.user_message);

    const status = node.querySelector(".status-pill");
    status.textContent = translateStatusWord(job.status);
    status.dataset.status = job.status;

    const logEl = node.querySelector(".job-log");
    logEl.textContent = "Open this panel to load logs...";

    const details = node.querySelector(".job-details");
    details.addEventListener("toggle", async () => {
      if (!details.open) return;
      try {
        const logs = await fetchJSON(`/api/jobs/${job.id}/logs`, { headers: buildHeaders() });
        logEl.textContent = logs.text || "No logs yet";
      } catch (error) {
        logEl.textContent = `Failed to load logs: ${error.message}`;
      }
    });

    const artifactList = node.querySelector(".artifact-list");
    artifactList.innerHTML = "";
    if (job.artifacts?.length) {
      for (const artifact of job.artifacts) {
        const link = document.createElement("a");
        link.className = "artifact-link";
        link.href = `/api/jobs/${job.id}/artifacts/${artifact.id}`;
        link.textContent = `${artifact.display_name} (${formatBytes(artifact.size_bytes)})`;
        artifactList.appendChild(link);
      }
    } else {
      artifactList.textContent = "No downloadable result files yet";
    }

    elements.jobsList.appendChild(node);
  }
}

function translateStatusWord(status) {
  const map = {
    queued: "Queued",
    running: "Running",
    succeeded: "Done",
    failed: "Failed",
  };
  return map[status] || status;
}

function translateStatus(status, fallback) {
  const map = {
    queued: "Your job has been received and is waiting to start.",
    running: "Your job is currently running.",
    succeeded: "Your job finished successfully. Files are ready to download.",
    failed: "This job did not finish successfully.",
  };
  return map[status] || fallback || "Processing";
}

async function refreshJobs() {
  if (!state.token) {
    state.jobs = [];
    elements.jobsList.innerHTML = "";
    elements.jobsEmpty.style.display = "block";
    elements.jobsEmpty.textContent = "Please log in to view your jobs.";
    return;
  }

  try {
    const payload = await fetchJSON("/api/jobs", { headers: buildHeaders() });
    state.jobs = payload.jobs || [];
    renderJobs();
  } catch (error) {
    state.jobs = [];
    elements.jobsList.innerHTML = "";
    elements.jobsEmpty.style.display = "block";
    elements.jobsEmpty.textContent = `Failed to load jobs: ${error.message}`;
  }
}

async function submitJob(event) {
  event.preventDefault();
  if (!state.token) {
    elements.submitStatus.textContent = "Please log in first.";
    return;
  }

  const formData = new FormData(elements.jobForm);
  elements.submitStatus.textContent = "Submitting...";
  try {
    const job = await fetchJSON("/api/jobs", {
      method: "POST",
      headers: buildHeaders(),
      body: formData,
    });
    elements.submitStatus.textContent = `Submitted: ${job.public_name}`;
    elements.jobForm.reset();
    setStructureFieldsVisibility();
    await refreshJobs();
  } catch (error) {
    elements.submitStatus.textContent = `Submission failed: ${error.message}`;
  }
}

async function handleAuthSubmit(event) {
  event.preventDefault();
  const payload = {
    email: elements.authEmail.value.trim(),
    password: elements.authPassword.value,
    display_name: elements.authDisplayName.value.trim() || null,
  };
  const endpoint = state.authMode === "register" ? "/api/auth/register" : "/api/auth/login";

  try {
    const data = await fetchJSON(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    state.token = data.access_token;
    state.user = data.user;
    writeStoredUser(data.user, data.access_token);
    elements.authPassword.value = "";
    elements.meOutput.textContent = describeUser(data.user);
    await refreshJobs();
  } catch (error) {
    elements.meOutput.textContent = `Authentication failed: ${error.message}`;
  }
}

function logout() {
  state.token = "";
  state.user = null;
  writeStoredUser(null, "");
  elements.meOutput.textContent = "You have been logged out.";
  refreshJobs();
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

function toggleAutoRefresh() {
  state.autoRefresh = !state.autoRefresh;
  elements.toggleAutoRefresh.textContent = `Auto refresh: ${state.autoRefresh ? "On" : "Off"}`;
}

function startPolling() {
  if (state.timerId) {
    window.clearInterval(state.timerId);
  }
  state.timerId = window.setInterval(() => {
    if (state.autoRefresh) {
      refreshJobs();
    }
  }, 10000);
}

elements.authForm.addEventListener("submit", handleAuthSubmit);
elements.logout.addEventListener("click", logout);
elements.tabLogin.addEventListener("click", () => setAuthMode("login"));
elements.tabRegister.addEventListener("click", () => setAuthMode("register"));
elements.refreshMe.addEventListener("click", refreshIdentity);
elements.refreshJobs.addEventListener("click", refreshJobs);
elements.toggleAutoRefresh.addEventListener("click", toggleAutoRefresh);
elements.jobForm.addEventListener("submit", submitJob);
elements.modeSelect.addEventListener("change", setStructureFieldsVisibility);

setAuthMode("login");
setStructureFieldsVisibility();
if (state.user && state.token) {
  elements.meOutput.textContent = describeUser(state.user);
}
refreshIdentity();
refreshJobs();
startPolling();
