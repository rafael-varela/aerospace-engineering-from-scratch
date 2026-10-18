// Builds the site into site/dist:
//   /                         the landing page (site/index.html), with finished phases linked
//   /curriculum/              every phase that exists in curriculum/
//   /curriculum/phase-NN/     one page per phase, from its README.md
//   /curriculum/phase-NN/lesson-MM-slug/
//                             one page per lesson, from its lesson.md and problems.md
// Every lesson page comes from the same template, so a new phase or lesson needs no UI work:
// add the Markdown, push, and the next deploy has its pages.
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";
import { lexer, renderTokens, renderInline } from "./markdown.mjs";
import * as T from "./templates.mjs";

const siteDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const repoDir = path.resolve(siteDir, "..");
const outDir = path.join(siteDir, "dist");
const curriculumDir = path.join(repoDir, "curriculum");

if (!fs.existsSync(curriculumDir)) {
  console.error(`No curriculum/ folder at ${curriculumDir}. The build reads the lessons from the repository root,`);
  console.error("so the host must include files outside site/ (Vercel: Settings → Build → 'Include files outside the root directory').");
  process.exit(1);
}

// ---- site name and repository, from config.js (the same file the landing page reads) ----
const sandbox = { window: {} };
vm.runInNewContext(fs.readFileSync(path.join(siteDir, "config.js"), "utf8"), sandbox);
const SITE = {
  name: sandbox.window.SITE?.name ?? "Aerospace Engineering from Scratch",
  repo: (sandbox.window.SITE?.repo ?? "https://github.com/rafael-varela/aerospace-engineering-from-scratch").replace(/\/$/, ""),
};

const read = (p) => fs.readFileSync(path.join(repoDir, p), "utf8");
const exists = (p) => fs.existsSync(path.join(repoDir, p));
const isDir = (p) => exists(p) && fs.statSync(path.join(repoDir, p)).isDirectory();
const ls = (p) => (isDir(p) ? fs.readdirSync(path.join(repoDir, p)).sort() : []);

// ---- discover phases and lessons ----
const phases = [];
const lessonsById = new Map();
const lessonsByDir = new Map();

for (const phaseDir of ls("curriculum").filter((d) => /^phase-\d\d$/.test(d))) {
  const readme = `curriculum/${phaseDir}/README.md`;
  if (!exists(readme)) continue;
  const num = phaseDir.slice(6);
  const tokens = lexer(read(readme));
  const h1 = tokens.find((t) => t.type === "heading" && t.depth === 1);
  const title = (h1?.text ?? `Phase ${num}`).replace(/^Phase\s+\d+\s*[:·.-]\s*/, "");
  const intro = tokens.find((t) => t.type === "paragraph");

  // The lesson table (header "# | Lesson | Build it | Mode") gives each lesson's one-liners.
  const table = tokens.find((t) => t.type === "table" && /lesson/i.test(t.header[1]?.text ?? ""));
  const rows = new Map();
  for (const r of table?.rows ?? []) {
    const m = /\]\(([^)]+?)\/lesson\.md\)/.exec(r[1].text);
    if (m) rows.set(m[1], { buildIt: r[2]?.text ?? "", mode: r[3]?.text ?? "" });
  }

  const phase = { num, dir: phaseDir, title, readme, tokens, table, intro, url: `/curriculum/${phaseDir}/`, lessons: [] };
  phase.exam = exists(`assessments/${phaseDir}`) ? `assessments/${phaseDir}` : null;
  for (const lessonDir of ls(`curriculum/${phaseDir}`).filter((d) => /^lesson-\d\d/.test(d))) {
    const src = `curriculum/${phaseDir}/${lessonDir}/lesson.md`;
    if (!exists(src)) continue;
    const firstLine = read(src).split("\n", 1)[0];
    const m = /^#\s+(P\d\d\.\d\d)\s+(.*)$/.exec(firstLine);
    const lesson = {
      phase,
      dir: lessonDir,
      src,
      num: lessonDir.slice(7, 9),
      id: m?.[1] ?? `P${num}.${lessonDir.slice(7, 9)}`,
      title: (m?.[2] ?? lessonDir).trim(),
      url: `/curriculum/${phaseDir}/${lessonDir}/`,
      ...(rows.get(lessonDir) ?? { buildIt: "", mode: "" }),
    };
    phase.lessons.push(lesson);
    lessonsById.set(lesson.id, lesson);
    lessonsByDir.set(`curriculum/${phaseDir}/${lessonDir}`, lesson);
  }
  phases.push(phase);
}

// ---- link resolution shared by every page ----
const mediaCopies = new Map();
const site = {
  repoDir,
  lessonById: (id) => lessonsById.get(id),
  urlFor(target, hash) {
    const t = target.replace(/\/$/, "");
    const h = hash ? `#${hash}` : "";
    let m = /^(curriculum\/phase-\d\d\/lesson-[^/]+)(?:\/(lesson|problems)\.md)?$/.exec(t);
    if (m && lessonsByDir.has(m[1])) return lessonsByDir.get(m[1]).url + (h || (m[2] === "problems" ? "#problem-set" : ""));
    m = /^curriculum\/(phase-\d\d)(?:\/README\.md)?$/.exec(t);
    if (m && phases.some((p) => p.dir === m[1])) return `/curriculum/${m[1]}/${h}`;
    if (/^curriculum(\/README\.md)?$/.test(t)) return `/curriculum/${h}`;
    if (/\.(png|jpe?g|gif|svg|webp)$/i.test(t) && exists(t)) {
      mediaCopies.set(t, `/media/${t}`);
      return `/media/${t}`;
    }
    return `${SITE.repo}/${isDir(t) ? "tree" : "blob"}/main/${t}${h}`;
  },
};
const gh = (p) => site.urlFor(p, "");
for (const p of phases) p.examUrl = p.exam ? gh(p.exam) : null;

// ---- page builders ----
function splitLesson(md) {
  // Peel the title, the one-line tagline and the metadata table off the front of lesson.md;
  // everything from the first ## heading on is the lesson body.
  const tokens = lexer(md);
  const firstH2 = tokens.findIndex((t) => t.type === "heading" && t.depth === 2);
  const head = firstH2 < 0 ? tokens : tokens.slice(0, firstH2);
  const body = firstH2 < 0 ? [] : tokens.slice(firstH2);
  const tagline = head.find((t) => t.type === "blockquote");
  const meta = head.find((t) => t.type === "table");
  const rest = head.filter((t) => t !== tagline && t !== meta && !(t.type === "heading" && t.depth === 1));
  body.links = tokens.links;
  rest.links = tokens.links;
  return { tagline, meta, rest, body };
}

function buildLesson(lesson, prev, next) {
  const { phase, src } = lesson;
  const { tagline, meta, rest, body } = splitLesson(read(src));
  const ids = new Set(["problem-set", "get-the-code", "spec"]);
  const toc = [];
  const opts = { src, site, ids, toc };

  const metaRows = [];
  for (const r of meta?.rows ?? []) {
    const key = r[0].text.replace(/\*\*/g, "").trim();
    metaRows.push({ key, text: r[1].text, html: renderInline(r[1].text, { src, site }) });
  }

  const intro = rest.length ? renderTokens(rest, opts) : "";
  const bodyHtml = renderTokens(body, opts);
  const lessonToc = toc.splice(0);

  let problemsHtml = "";
  const problemsSrc = src.replace(/lesson\.md$/, "problems.md");
  if (exists(problemsSrc)) {
    const ptokens = lexer(read(problemsSrc)).filter((t) => !(t.type === "heading" && t.depth === 1));
    problemsHtml = renderTokens(ptokens, { src: problemsSrc, site, ids, toc: [], shift: 1 });
  }

  const solutionsDir = `solutions/${phase.dir}/${lesson.dir}`;
  const html = T.lessonPage({
    SITE,
    lesson,
    phase,
    prev,
    next,
    taglineHtml: tagline ? renderTokens(tagline.tokens, { src, site, ids: new Set() }) : "",
    metaRows,
    introHtml: intro,
    bodyHtml,
    toc: lessonToc,
    problemsHtml,
    links: {
      starter: exists(`curriculum/${phase.dir}/${lesson.dir}/starter`) ? gh(`curriculum/${phase.dir}/${lesson.dir}/starter`) : null,
      tests: exists(`curriculum/${phase.dir}/${lesson.dir}/tests`) ? gh(`curriculum/${phase.dir}/${lesson.dir}/tests`) : null,
      source: gh(src),
      problemsSource: exists(problemsSrc) ? gh(problemsSrc) : null,
      solutions: exists(solutionsDir) ? gh(solutionsDir) : null,
    },
  });
  write(lesson.url, html);
}

function buildPhase(phase) {
  // README minus its h1, its intro paragraph and the lesson table (shown as cards instead).
  const tokens = phase.tokens.filter((t) => t !== phase.table && t !== phase.intro && !(t.type === "heading" && t.depth === 1));
  tokens.links = phase.tokens.links;
  const ids = new Set(["lessons"]);
  const restHtml = renderTokens(tokens, { src: phase.readme, site, ids, toc: [] });
  const html = T.phasePage({
    SITE,
    phase,
    introHtml: phase.intro ? renderInline(phase.intro.text, { src: phase.readme, site }) : "",
    lessonCards: phase.lessons.map((l) => ({ ...l, buildItHtml: renderInline(l.buildIt, { src: phase.readme, site }) })),
    restHtml,
    sourceUrl: gh(phase.readme),
  });
  write(phase.url, html);
}

function buildIndex() {
  const html = T.curriculumPage({
    SITE,
    phases: phases.map((p) => ({ ...p, introHtml: p.intro ? renderInline(firstSentence(p.intro.text), { src: p.readme, site }) : "" })),
    roadmapUrl: gh("ROADMAP.md"),
  });
  write("/curriculum/", html);
}

const firstSentence = (s) => (/^(.*?[.!?])(\s|$)/s.exec(s)?.[1] ?? s);

// The landing page is hand-written; the build only links every phase row that has pages
// and fills the status line, so finishing a phase needs no edit to index.html.
function buildLanding() {
  let html = fs.readFileSync(path.join(siteDir, "index.html"), "utf8");
  const byNum = new Map(phases.map((p) => [p.num, p]));
  html = html.replace(
    /<(a|div) class="phase(?: ready)?"[^>]*>\s*<span class="n">(\d\d)<\/span>\s*(<span class="t">.*?<\/span>)\s*<span class="l">[^<]*<\/span>\s*<\/\1>/g,
    (whole, _tag, num, title) => {
      const p = byNum.get(num);
      if (!p) return whole;
      return `<a class="phase ready" href="${p.url}"><span class="n">${num}</span>${title}<span class="l">${p.lessons.length} · open →</span></a>`;
    },
  );
  const total = phases.reduce((n, p) => n + p.lessons.length, 0);
  const range = phases.length ? `${phases[0].num}–${phases.at(-1).num}` : "";
  html = html.replace(
    /(<span data-build="status">)[\s\S]*?(<\/span>)/,
    `$1Status: Phases ${range} are open, ${total} lessons. Later phases are being written in roadmap order.$2`,
  );
  write("/", html);
}

// ---- output ----
const warnings = [];
function write(url, html) {
  const bad = html.match(/class="katex-error"/g);
  if (bad) warnings.push(`${url}: ${bad.length} formula(s) KaTeX could not parse; shown as raw TeX`);
  const file = path.join(outDir, url, "index.html");
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, html);
}

function copy(from, to) {
  fs.mkdirSync(path.dirname(to), { recursive: true });
  fs.cpSync(from, to, { recursive: true });
}

fs.rmSync(outDir, { recursive: true, force: true });
fs.mkdirSync(outDir, { recursive: true });

for (const phase of phases) {
  buildPhase(phase);
  phase.lessons.forEach((l, i) => buildLesson(l, phase.lessons[i - 1] ?? null, phase.lessons[i + 1] ?? null));
}
buildIndex();
buildLanding();

copy(path.join(siteDir, "config.js"), path.join(outDir, "config.js"));
copy(path.join(siteDir, "assets"), path.join(outDir, "assets"));
const katexDist = path.join(siteDir, "node_modules", "katex", "dist");
copy(path.join(katexDist, "katex.min.css"), path.join(outDir, "assets", "katex", "katex.min.css"));
copy(path.join(katexDist, "fonts"), path.join(outDir, "assets", "katex", "fonts"));
for (const [srcPath, url] of mediaCopies) copy(path.join(repoDir, srcPath), path.join(outDir, url));

const count = phases.reduce((n, p) => n + p.lessons.length, 0);
for (const w of warnings) console.warn(`warning: ${w}`);
console.log(`Built ${phases.length} phases and ${count} lesson pages into ${path.relative(repoDir, outDir)}/`);
