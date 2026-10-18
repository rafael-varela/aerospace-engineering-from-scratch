// Small enhancements for the generated pages. Everything works without it.
(function () {
  // Copy buttons on code blocks.
  document.querySelectorAll("pre.code").forEach(function (pre) {
    if (!navigator.clipboard) return;
    var b = document.createElement("button");
    b.type = "button";
    b.className = "copy";
    b.textContent = "Copy";
    b.addEventListener("click", function () {
      navigator.clipboard.writeText(pre.querySelector("code").innerText).then(function () {
        b.textContent = "Copied";
        setTimeout(function () { b.textContent = "Copy"; }, 1400);
      });
    });
    pre.appendChild(b);
  });

  // On narrow screens the lesson list starts collapsed; keep the current lesson in view.
  var side = document.querySelector(".side-box");
  if (side) {
    if (window.matchMedia("(max-width: 920px)").matches) side.open = false;
    else {
      side.addEventListener("click", function (e) { if (e.target.closest("summary")) e.preventDefault(); });
      var cur = side.querySelector('[aria-current="page"]');
      var pane = side.parentElement;
      if (cur) pane.scrollTop = cur.offsetTop - pane.offsetTop - pane.clientHeight / 2;
    }
  }

  // Highlight the section being read in "On this page".
  var links = Array.prototype.slice.call(document.querySelectorAll(".toc a"));
  var targets = links.map(function (a) { return document.getElementById(a.getAttribute("href").slice(1)); });
  if (!links.length) return;
  var ticking = false;
  function spy() {
    ticking = false;
    var current = -1;
    for (var i = 0; i < targets.length; i++) {
      if (targets[i] && targets[i].getBoundingClientRect().top < window.innerHeight * 0.3) current = i;
    }
    if (window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 4) current = targets.length - 1;
    links.forEach(function (a, i) { a.classList.toggle("active", i === current); });
  }
  window.addEventListener("scroll", function () { if (!ticking) { ticking = true; requestAnimationFrame(spy); } }, { passive: true });
  spy();
})();
