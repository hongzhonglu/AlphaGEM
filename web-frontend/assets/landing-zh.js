const phrases = [
  "结构感知的模型重建",
  "PLM 引导的暗代谢发现",
  "从序列走向 GEM",
  "精确、可追踪、可网页化的建模流程",
];

const phraseNode = document.getElementById("rotating-phrase");
let phraseIndex = 0;

if (phraseNode) {
  window.setInterval(() => {
    phraseNode.classList.add("is-fading");
    window.setTimeout(() => {
      phraseIndex = (phraseIndex + 1) % phrases.length;
      phraseNode.textContent = phrases[phraseIndex];
      phraseNode.classList.remove("is-fading");
    }, 220);
  }, 2600);
}
