Reveal.initialize({
  hash: true,
  history: true,
  margin: 0.06,
  controls: true,
  progress: true,
  slideNumber: "c/t",
  transition: "fade",
  backgroundTransition: "fade",
  plugins: [RevealMath.MathJax3],
  mathjax3: {
    tex: {
      inlineMath: [["\\(", "\\)"]],
      displayMath: [["\\[", "\\]"]]
    }
  }
});
