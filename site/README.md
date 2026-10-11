# Landing page

A single static page, no build step. Open `index.html` in a browser or serve the folder with any static host (GitHub Pages pointed at `/site`, Netlify, Cloudflare Pages).

- `index.html` holds the content, styles, and the interactive airfoil.
- `config.js` holds the project name and repository URL. Change the name there once; the page fills it in everywhere. The fallback text in `index.html` (search for `data-site-name`) is only for readers without JavaScript, so update it too when the name is final.

Content is drawn from `ROADMAP.md` and `README.md` at the repository root. Keep the page's claims in step with them: scope comparable to a degree, competence demonstrated rather than conferred, no credential, and no total time estimate until pilots have measured one.
