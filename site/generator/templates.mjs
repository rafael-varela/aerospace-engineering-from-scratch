// HTML templates. Every lesson page, every phase page and the curriculum index are rendered by
// these functions; the styles live in site/assets/lesson.css and share the landing page's tokens.
import { esc } from "./markdown.mjs";

const strip = (html) => html.replace(/<[^>]+>/g, "").replace(/\s+/g, " ").trim();

function page({ SITE, title, description, body, math = false }) {
  return `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>${esc(title)}</title>
<meta name="description" content="${esc(description)}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wdth,wght@100,500..900&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&family=JetBrains+Mono:wght@400;500&display=swap">
${math ? '<link rel="stylesheet" href="/assets/katex/katex.min.css">\n' : ""}<link rel="stylesheet" href="/assets/lesson.css">
<script src="/assets/lesson.js" defer></script>
</head>
<body>
<a class="skip" href="#content">Skip to content</a>
<header class="nav">
  <div class="nav-in">
    <a class="brand" href="/">
      <svg viewBox="0 0 32 32" aria-hidden="true"><path d="M3 22 L29 10" stroke="currentColor" stroke-width="2" fill="none"/><path d="M3 22 C 10 12, 22 8, 29 10 C 24 14, 12 20, 3 22 Z" fill="var(--orange)"/></svg>
      <span>${esc(SITE.name)}</span>
    </a>
    <nav class="nav-links" aria-label="Site">
      <a href="/curriculum/">Curriculum</a>
      <a href="${SITE.repo}/blob/main/ROADMAP.md">Roadmap</a>
      <a class="gh" href="${SITE.repo}">GitHub ↗</a>
    </nav>
  </div>
</header>
${body}
<footer class="foot">
  <div class="foot-in">
    <span>${esc(SITE.name)}. Code MIT, text CC BY-SA 4.0. Scope comparable to a degree; no credential is conferred.</span>
    <a href="${SITE.repo}">Source on GitHub</a>
  </div>
</footer>
</body>
</html>
`;
}

const crumbs = (items) =>
  `<nav class="crumbs" aria-label="Breadcrumb"><ol>${items
    .map(([label, url]) => (url ? `<li><a href="${url}">${label}</a></li>` : `<li aria-current="page">${label}</li>`))
    .join("")}</ol></nav>`;

const modeLabel = (mode) => {
  const m = strip(mode).toLowerCase();
  if (!m) return "";
  const kind = m.includes("home") || m.includes("physical") || m.includes("🔬") ? "lab" : m.includes("data") ? "data" : "comp";
  const full = strip(mode);
  const short = full.split(/[:,;(]/)[0].trim();
  return `<span class="mode mode-${kind}"${short !== full ? ` title="${esc(full)}"` : ""}>${esc(short)}</span>`;
};

// "1. Do this. 2. Do that." → an ordered list, when the cell is written that way.
function outcomesList(html) {
  const parts = html.split(/(?:^|\s)\d{1,2}\.\s(?=[A-Z★`*_<])/).map((s) => s.trim()).filter(Boolean);
  if (parts.length < 2) return `<p>${html}</p>`;
  return `<ol>${parts.map((p) => `<li>${p}</li>`).join("")}</ol>`;
}

const PRIMARY = ["Learning outcomes", "Prerequisites", "Engineering question"];

function spec(metaRows) {
  const get = (k) => metaRows.find((r) => r.key === k);
  const primary = PRIMARY.map(get).filter(Boolean);
  const safety = get("Safety");
  const safetyNeeded = safety && !/^none\b/i.test(strip(safety.html));
  const others = metaRows.filter((r) => !PRIMARY.includes(r.key) && !(safetyNeeded && r === safety));
  return `<section class="spec" id="spec" aria-label="Lesson specification">
  <dl class="spec-main">
    ${primary
      .map((r) => `<div class="spec-row spec-${r.key.split(" ")[0].toLowerCase()}"><dt>${esc(r.key)}</dt><dd>${r.key === "Learning outcomes" ? outcomesList(r.html) : `<p>${r.html}</p>`}</dd></div>`)
      .join("\n    ")}
  </dl>
  ${safetyNeeded ? `<aside class="callout callout-safety"><p class="callout-k">Safety</p><p>${safety.html}</p></aside>` : ""}
  <details class="spec-more">
    <summary>Full lesson specification <span class="muted">(${others.length} more items)</span></summary>
    <dl>
      ${others.map((r) => `<div class="spec-row"><dt>${esc(r.key)}</dt><dd><p>${r.html}</p></dd></div>`).join("\n      ")}
    </dl>
  </details>
</section>`;
}

function sidebar(phase, current) {
  return `<aside class="side" aria-label="Lessons in Phase ${phase.num}">
  <details class="side-box" open>
    <summary><span class="eyebrow">Phase ${phase.num}</span><span class="side-title">${esc(phase.title)}</span></summary>
    <ol class="side-list">
      ${phase.lessons
        .map(
          (l) =>
            `<li><a href="${l.url}"${l === current ? ' aria-current="page"' : ""}><span class="n">${l.num}</span><span class="t">${esc(l.title)}</span></a></li>`,
        )
        .join("\n      ")}
    </ol>
    <a class="side-phase" href="${phase.url}">Phase overview, exam and kit →</a>
  </details>
</aside>`;
}

function tocAside(toc, hasProblems) {
  const items = [...toc];
  if (hasProblems) items.push({ id: "problem-set", text: "Problem set" });
  items.push({ id: "get-the-code", text: "Get the code" });
  return `<aside class="toc" aria-label="On this page">
  <p class="eyebrow">On this page</p>
  <ol>${items.map((t) => `<li><a href="#${t.id}">${esc(t.text)}</a></li>`).join("")}</ol>
</aside>`;
}

export function lessonPage(d) {
  const { SITE, lesson, phase, prev, next, links } = d;
  const time = d.metaRows.find((r) => r.key === "Time");
  const description = strip(d.taglineHtml) || `${lesson.id} ${lesson.title}`;
  const repoPath = `curriculum/${phase.dir}/${lesson.dir}`;
  const body = `
<div class="shell" id="content">
${sidebar(phase, lesson)}
<main class="article">
  ${crumbs([["Curriculum", "/curriculum/"], [`Phase ${phase.num} · ${esc(phase.title)}`, phase.url], [lesson.id, null]])}
  <header class="lesson-head">
    <p class="eyebrow"><span class="mono">${esc(lesson.id)}</span> · Lesson ${lesson.num} of ${phase.lessons.length}</p>
    <h1>${esc(lesson.title)}</h1>
    ${d.taglineHtml ? `<div class="tagline">${d.taglineHtml}</div>` : ""}
    <div class="chips">
      ${modeLabel(lesson.mode)}
      ${time ? `<span class="chip"><b>Time</b> ${esc(strip(time.html))}</span>` : ""}
    </div>
    ${lesson.buildIt ? `<p class="buildit"><span class="k">Build it</span> ${esc(strip(lesson.buildIt))}</p>` : ""}
  </header>
  ${d.metaRows.length ? spec(d.metaRows) : ""}
  ${d.introHtml ? `<div class="prose">${d.introHtml}</div>` : ""}
  <div class="prose">
${d.bodyHtml}
  </div>
  ${
    d.problemsHtml
      ? `<section class="prose problems" aria-labelledby="problem-set">
    <h2 id="problem-set"><a class="anchor" href="#problem-set" aria-hidden="true" tabindex="-1">#</a>Problem set</h2>
${d.problemsHtml}
    <aside class="callout">
      <p class="callout-k">Solutions</p>
      <p>Attempt every problem before you look. ${
        links.solutions
          ? `Hints, the method and worked solutions are in <a href="${links.solutions}">solutions/${phase.dir}/${lesson.dir}/</a>, in that order: open the hints first.`
          : "Worked solutions for this lesson are not published yet."
      }</p>
    </aside>
  </section>`
      : ""
  }
  <section class="codebox" aria-labelledby="get-the-code">
    <h2 id="get-the-code">Get the code</h2>
    <p>Clone the repository once, write your code in <code>starter/main.py</code>, then run the tests from the repository root.</p>
    <pre class="code"><code>pytest ${repoPath}</code></pre>
    <p class="muted">To run the same tests against the reference solution instead:</p>
    <pre class="code"><code>AEFS_IMPL=solution pytest ${repoPath}</code></pre>
    <div class="btn-row">
      ${links.starter ? `<a class="btn btn-primary" href="${links.starter}">Starter code ↗</a>` : ""}
      ${links.tests ? `<a class="btn btn-ghost" href="${links.tests}">Tests ↗</a>` : ""}
      <a class="btn btn-ghost" href="${links.source}">Lesson source ↗</a>
    </div>
    <p class="muted small">Found a mistake? Every page is built from the Markdown in the repository; <a href="${links.source}">edit it on GitHub</a>.</p>
  </section>
  <nav class="pager" aria-label="Lessons">
    ${prev ? `<a class="prev" href="${prev.url}"><span class="k">← Previous</span><span class="t">${esc(prev.id)} ${esc(prev.title)}</span></a>` : `<a class="prev" href="${phase.url}"><span class="k">← Phase ${phase.num}</span><span class="t">${esc(phase.title)}</span></a>`}
    ${next ? `<a class="next" href="${next.url}"><span class="k">Next →</span><span class="t">${esc(next.id)} ${esc(next.title)}</span></a>` : `<a class="next" href="${phase.examUrl ?? phase.url}"><span class="k">Phase ${phase.num} done →</span><span class="t">Exam and mini-project</span></a>`}
  </nav>
</main>
${tocAside(d.toc, !!d.problemsHtml)}
</div>`;
  return page({ SITE, title: `${lesson.id} ${lesson.title} · ${SITE.name}`, description, body, math: true });
}

export function phasePage(d) {
  const { SITE, phase } = d;
  const first = phase.lessons[0];
  const body = `
<main class="wide" id="content">
  ${crumbs([["Curriculum", "/curriculum/"], [`Phase ${phase.num}`, null]])}
  <header class="phase-head">
    <p class="eyebrow">Phase ${phase.num} · ${phase.lessons.length} lessons</p>
    <h1>${esc(phase.title)}</h1>
    <p class="lede">${d.introHtml}</p>
    <div class="btn-row">
      ${first ? `<a class="btn btn-primary" href="${first.url}">Start with lesson ${first.num}</a>` : ""}
      ${phase.examUrl ? `<a class="btn btn-ghost" href="${phase.examUrl}">Exam and mini-project ↗</a>` : ""}
    </div>
  </header>
  <section aria-labelledby="lessons">
    <h2 id="lessons" class="section-h">Lessons</h2>
    <ol class="cards">
      ${d.lessonCards
        .map(
          (l) => `<li><a class="card" href="${l.url}">
        <span class="card-top"><span class="mono">${esc(l.id)}</span>${modeLabel(l.mode)}</span>
        <span class="card-t">${esc(l.title)}</span>
        ${l.buildItHtml ? `<span class="card-b"><b>Build it:</b> ${strip(l.buildItHtml)}</span>` : ""}
      </a></li>`,
        )
        .join("\n      ")}
    </ol>
  </section>
  <div class="prose prose-wide">
${d.restHtml}
  </div>
  <p class="muted small"><a href="${d.sourceUrl}">This page on GitHub</a></p>
</main>`;
  return page({ SITE, title: `Phase ${phase.num}: ${phase.title} · ${SITE.name}`, description: strip(d.introHtml).slice(0, 300), body, math: true });
}

export function curriculumPage(d) {
  const { SITE, phases } = d;
  const total = phases.reduce((n, p) => n + p.lessons.length, 0);
  const body = `
<main class="wide" id="content">
  <header class="phase-head">
    <p class="eyebrow">Curriculum</p>
    <h1>Every open lesson</h1>
    <p class="lede">${phases.length} phases and ${total} lessons are written so far. Each lesson goes problem → derive → model → experiment → validate → defend, with tested code and a problem set. The later phases are being written in <a href="${d.roadmapUrl}">roadmap</a> order.</p>
  </header>
  <ol class="phase-cards">
    ${phases
      .map(
        (p) => `<li><a class="card phase-card" href="${p.url}">
      <span class="card-top"><span class="mono">Phase ${p.num}</span><span class="muted mono">${p.lessons.length} lessons</span></span>
      <span class="card-t">${esc(p.title)}</span>
      <span class="card-b">${strip(p.introHtml)}</span>
    </a></li>`,
      )
      .join("\n    ")}
  </ol>
</main>`;
  return page({ SITE, title: `Curriculum · ${SITE.name}`, description: `All ${total} open lessons of ${SITE.name}.`, body, math: true });
}
