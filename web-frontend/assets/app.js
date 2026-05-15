const state = {
  autoRefresh: true,
  timerId: null,
  jobs: [],
  authMode: "login",
  token: localStorage.getItem("alphagem_token") || "",
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
    return JSON.parse(localStorage.getItem("alphagem_user") || "null");
  } catch {
    return null;
  }
}

function writeStoredUser(user, token) {
  if (user) {
    localStorage.setItem("alphagem_user", JSON.stringify(user));
  } else {
    localStorage.removeItem("alphagem_user");
  }
  if (token) {
    localStorage.setItem("alphagem_token", token);
  } else {
    localStorage.removeItem("alphagem_token");
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
  elements.authSubmit.textContent = mode === "register" ? "创建账户" : "登录";
  elements.authHint.textContent =
    mode === "register"
      ? "创建账户后即可提交任务。"
      : "登录后可查看任务、日志和结果文件。";
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
  return `已登录：${user.display_name || user.email}`;
}

async function refreshIdentity() {
  if (!state.token) {
    state.user = null;
    elements.meOutput.textContent = "请先登录或注册。";
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
    elements.meOutput.textContent = `身份验证失败：${error.message}`;
  }
}

function renderJobs() {
  elements.jobsList.innerHTML = "";
  elements.jobsEmpty.style.display = state.jobs.length ? "none" : "block";

  if (!state.jobs.length) {
    elements.jobsEmpty.textContent = "还没有任务。登录后上传 FASTA 文件即可创建新的 AlphaGEM 任务。";
  }

  for (const job of state.jobs) {
    const node = elements.template.content.firstElementChild.cloneNode(true);
    node.querySelector(".job-name").textContent = job.public_name;
    node.querySelector(".job-meta").textContent = `${job.mode} · ${job.refname} · ${job.growthMedium}`;
    node.querySelector(".job-id").textContent = job.id.slice(0, 8);
    node.querySelector(".job-created").textContent = new Date(job.created_at).toLocaleString();
    node.querySelector(".job-message").textContent = job.user_message || "处理中";

    const status = node.querySelector(".status-pill");
    status.textContent = job.status;
    status.dataset.status = job.status;

    const logEl = node.querySelector(".job-log");
    logEl.textContent = "点击展开后自动加载日志...";

    const details = node.querySelector(".job-details");
    details.addEventListener("toggle", async () => {
      if (!details.open) return;
      try {
        const logs = await fetchJSON(`/api/jobs/${job.id}/logs`, { headers: buildHeaders() });
        logEl.textContent = logs.text || "暂无日志";
      } catch (error) {
        logEl.textContent = `日志加载失败：${error.message}`;
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
      artifactList.textContent = "暂无可下载结果文件";
    }

    elements.jobsList.appendChild(node);
  }
}

async function refreshJobs() {
  if (!state.token) {
    state.jobs = [];
    elements.jobsList.innerHTML = "";
    elements.jobsEmpty.style.display = "block";
    elements.jobsEmpty.textContent = "请先登录后查看任务。";
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
    elements.jobsEmpty.textContent = `任务加载失败：${error.message}`;
  }
}

async function submitJob(event) {
  event.preventDefault();
  if (!state.token) {
    elements.submitStatus.textContent = "请先登录。";
    return;
  }

  const formData = new FormData(elements.jobForm);
  elements.submitStatus.textContent = "提交中...";
  try {
    const job = await fetchJSON("/api/jobs", {
      method: "POST",
      headers: buildHeaders(),
      body: formData,
    });
    elements.submitStatus.textContent = `任务已提交：${job.public_name}`;
    elements.jobForm.reset();
    setStructureFieldsVisibility();
    await refreshJobs();
  } catch (error) {
    elements.submitStatus.textContent = `提交失败：${error.message}`;
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
    elements.meOutput.textContent = `认证失败：${error.message}`;
  }
}

function logout() {
  state.token = "";
  state.user = null;
  writeStoredUser(null, "");
  elements.meOutput.textContent = "你已退出登录。";
  refreshJobs();
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

function toggleAutoRefresh() {
  state.autoRefresh = !state.autoRefresh;
  elements.toggleAutoRefresh.textContent = `自动刷新: ${state.autoRefresh ? "开" : "关"}`;
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
