const phrases = [
  "Structure-aware reconstruction",
  "PLM-guided dark metabolism discovery",
  "From sequence to GEM",
  "Precise, trackable, web-based modeling",
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
