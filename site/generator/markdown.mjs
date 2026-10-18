// Markdown → HTML for curriculum files.
// One Marked instance, with a per-file context (set by renderFile) that the renderer uses to
// resolve relative links, give headings stable ids, and link lesson ids like `P10.03`.
import fs from "node:fs";
import path from "node:path";
import { Marked } from "marked";
import katex from "katex";
import hljs from "highlight.js/lib/core";
import python from "highlight.js/lib/languages/python";
import bash from "highlight.js/lib/languages/bash";
import yaml from "highlight.js/lib/languages/yaml";
import json from "highlight.js/lib/languages/json";
import ini from "highlight.js/lib/languages/ini";

hljs.registerLanguage("python", python);
hljs.registerLanguage("bash", bash);
hljs.registerLanguage("yaml", yaml);
hljs.registerLanguage("json", json);
hljs.registerLanguage("ini", ini);
hljs.registerAliases(["py"], { languageName: "python" });
hljs.registerAliases(["sh", "shell", "console", "zsh"], { languageName: "bash" });
hljs.registerAliases(["yml"], { languageName: "yaml" });
hljs.registerAliases(["toml"], { languageName: "ini" });

export const esc = (s) =>
  String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");

// ---- per-file context ----
let ctx = null;

export function slugify(text) {
  return text
    .toLowerCase()
    .replace(/<[^>]+>/g, "")
    .replace(/&[a-z]+;|&#\d+;/g, "")
    .replace(/[^\p{L}\p{N}\s-]/gu, "")
    .trim()
    .replace(/\s+/g, "-");
}

function uniqueId(text) {
  let base = slugify(text) || "section";
  let id = base;
  for (let i = 1; ctx.ids.has(id); i++) id = `${base}-${i}`;
  ctx.ids.add(id);
  return id;
}

// Resolve a repo-relative link from the file being rendered to a site URL or a GitHub URL.
function resolveHref(href) {
  if (!href || /^[a-z][a-z0-9+.-]*:/i.test(href) || href.startsWith("#") || href.startsWith("//")) return href;
  const [p, hash = ""] = href.split("#");
  const target = path.posix.normalize(path.posix.join(path.posix.dirname(ctx.src), decodeURI(p)));
  if (target.startsWith("..")) return href;
  return ctx.site.urlFor(target, hash) ?? href;
}

function mathHtml(tex, displayMode) {
  try {
    return katex.renderToString(tex, { displayMode, throwOnError: false, output: "htmlAndMathml", strict: "ignore" });
  } catch {
    return `<code>${esc(tex)}</code>`;
  }
}

const mathBlock = {
  name: "mathBlock",
  level: "block",
  start(src) {
    const m = /(^|\n)[ \t]*\$\$/.exec(src);
    return m ? m.index + m[1].length : undefined;
  },
  // "$$ … $$" on its own lines, optionally followed by an equation number: "$$ … $$   (1)".
  tokenizer(src) {
    const m = /^[ \t]*\$\$([\s\S]+?)\$\$[ \t]*(?:\((\d+[a-z]?)\))?[ \t]*(?:\n|$)/.exec(src);
    if (m) return { type: "mathBlock", raw: m[0], text: m[1].trim(), tag: m[2] };
  },
  renderer: (t) => `<div class="math-display">${mathHtml(t.tag ? `${t.text} \\tag{${t.tag}}` : t.text, true)}</div>\n`,
};

// $…$ inline. The opening $ must not be followed by a space and the closing one must not be
// preceded by a space or followed by a digit, so "$5 and $10" and shell prompts stay text.
const mathInline = {
  name: "mathInline",
  level: "inline",
  start(src) {
    const i = src.indexOf("$");
    return i < 0 ? undefined : i;
  },
  tokenizer(src) {
    let m = /^\$\$([^$]+?)\$\$/.exec(src);
    if (m) return { type: "mathInline", raw: m[0], text: m[1].trim(), display: true };
    m = /^\$(?!\s)((?:\\\$|[^$\n])+?)(?<![\s\\])\$(?![\d$])/.exec(src);
    if (m) return { type: "mathInline", raw: m[0], text: m[1], display: false };
  },
  renderer: (t) => mathHtml(t.text, t.display),
};

const LESSON_ID = /^P(\d\d)\.(\d\d)$/;
const PATHLIKE = /^[\w.-]+(\/[\w.<>-]*)+\/?$|^[\w-]+\.(py|md|yaml|yml|csv|json|txt|toml)$/;

const marked = new Marked({ gfm: true });
marked.use({
  extensions: [mathBlock, mathInline],
  renderer: {
    heading({ tokens, depth }) {
      const level = Math.min(6, depth + ctx.shift);
      const html = this.parser.parseInline(tokens);
      const plain = html.replace(/<[^>]+>/g, "");
      const id = uniqueId(plain);
      if (level === 2) ctx.toc.push({ id, text: plain });
      return `<h${level} id="${id}"><a class="anchor" href="#${id}" aria-hidden="true" tabindex="-1">#</a>${html}</h${level}>\n`;
    },
    link({ href, title, tokens }) {
      const url = resolveHref(href);
      const ext = /^https?:/.test(url);
      const t = title ? ` title="${esc(title)}"` : "";
      return `<a href="${esc(url)}"${t}${ext ? ' rel="noopener"' : ""}>${this.parser.parseInline(tokens)}</a>`;
    },
    image({ href, title, text }) {
      const url = resolveHref(href);
      const t = title ? ` title="${esc(title)}"` : "";
      return `<img src="${esc(url)}" alt="${esc(text)}"${t} loading="lazy">`;
    },
    codespan({ text }) {
      // marked has already escaped the text; undo it to test it against ids and paths.
      const raw = text.replace(/&amp;/g, "&").replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&quot;/g, '"').replace(/&#39;/g, "'");
      const id = LESSON_ID.exec(raw);
      if (id) {
        const lesson = ctx.site.lessonById(raw);
        if (lesson) return `<a class="lref" href="${lesson.url}" title="${esc(lesson.title)}"><code>${text}</code></a>`;
      }
      if (raw.startsWith("@")) return `<code class="cite">${text}</code>`;
      if (PATHLIKE.test(raw) && !raw.includes("<")) {
        const fromFile = path.posix.normalize(path.posix.join(path.posix.dirname(ctx.src), raw));
        for (const cand of [fromFile, path.posix.normalize(raw)]) {
          if (!cand.startsWith("..") && fs.existsSync(path.join(ctx.site.repoDir, cand))) {
            return `<a class="pref" href="${esc(ctx.site.urlFor(cand, ""))}"><code>${text}</code></a>`;
          }
        }
      }
      return `<code>${text}</code>`;
    },
    code({ text, lang }) {
      const l = (lang || "").trim().split(/\s+/)[0].toLowerCase();
      if (l && hljs.getLanguage(l)) {
        const html = hljs.highlight(text, { language: l, ignoreIllegals: true }).value;
        return `<pre class="code" data-lang="${esc(l)}"><code class="hljs">${html}</code></pre>\n`;
      }
      return `<pre class="code"><code>${esc(text)}</code></pre>\n`;
    },
  },
});

export const lexer = (md) => marked.lexer(md);

// Render tokens (or a Markdown string) as if they came from repo file `src`.
// `shift` demotes headings; `ids` and `toc` are shared across the parts of one page.
export function renderTokens(tokens, { src, site, shift = 0, ids, toc }) {
  const prev = ctx;
  ctx = { src, site, shift, ids: ids ?? new Set(), toc: toc ?? [] };
  try {
    return wrapTables(marked.parser(tokens));
  } finally {
    ctx = prev;
  }
}

const wrapTables = (html) =>
  html.replace(/<table>/g, '<div class="table-wrap"><table>').replace(/<\/table>/g, "</table></div>");

// Render one line of inline Markdown (a table cell, a tagline) from repo file `src`.
export function renderInline(text, { src, site }) {
  const prev = ctx;
  ctx = { src, site, shift: 0, ids: new Set(), toc: [] };
  try {
    return marked.parseInline(text);
  } finally {
    ctx = prev;
  }
}
