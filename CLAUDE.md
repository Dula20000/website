# CLAUDE.md

Instructions for Claude when working in this repo.

## What this is

`dulran.com`, a personal site owned by Dulran Kariyawasam. Plain static HTML
served by **GitHub Pages** from the **root of `main`**. The custom domain comes
from the `CNAME` file.

There is no build step, no framework, no package manager, no bundler. A page is
one `.html` file at the repo root, and pushing to `main` deploys it. A new page
at `/foo.html` is reachable at `https://dulran.com/foo.html` within a minute or
two of the push.

Pages are mostly self-contained: styles live in a `<style>` block in the file
rather than in a shared stylesheet. Follow that pattern unless told otherwise.
Shared assets that do exist: `styles/`, `scripts/`, `assets/`, `logofiles/`.

## Non-negotiable: the tracker goes on every page

Every HTML page must carry this immediately before `</body>`, comments and all:

```html
<!--this is the tracking code snippet from cloudfair worker and worker Kv. It uses FingerprintJS to log everything needed to fingerprint users. have this in all html files going forward.-->
<script src="https://dulran.com/tracker.js"></script>
<!--this is the tracking code snippet from cloudfair worker and worker Kv. It uses FingerprintJS to log everything needed to fingerprint users. have this in all html files going forward.-->
```

It reports to a Cloudflare Worker backed by Worker KV and uses FingerprintJS.
Add it to any new page without being asked. When editing an existing page, check
it is present and add it if it is missing.

## Media does not go in git

The repo is already heavy (`.git` is over 200 MB, and individual tracked videos
run to 25 MB). More importantly, **anything committed stays recoverable from the
history forever**, even after deletion.

So:

- Large or private media belongs in **Cloudflare R2**, served from
  **`media.dulran.com`**, not committed here.
- Check `.gitignore` before adding any video. Several patterns already keep
  clips out of the repo deliberately.
- If the owner says a file should be deletable later with no trace, that rules
  out git entirely. Use R2, give the object an unguessable name, and remember
  that the CDN cache needs purging after the object is deleted.
- Pages that depend on external media should degrade gracefully when the file is
  gone rather than showing a broken player. `pronunciation.html` is the worked
  example: it falls back to a "redacted" card.

## Crawlers are blocked site-wide

`robots.txt` disallows every user agent, including GPTBot, ClaudeBot, CCBot,
Google-Extended and PerplexityBot. Do not add sitemaps, SEO meta, or anything
premised on the site being indexed. For pages that should be extra quiet, also
set `<meta name="robots" content="noindex, nofollow, noarchive, nosnippet">`.

## Writing style

In all prose, on pages and in commits alike, **do not use em dashes** or the
other usual tells of machine-written text. This is not a request for simpler
words. It is a request for deliberate ones.

## Typography already in use

Loaded from Google Fonts by `@import` inside the page's own `<style>`:
Space Grotesk, JetBrains Mono, Inter, Fraunces, Cormorant Garamond,
DM Serif Display, IBM Plex Mono. Prefer one of these over introducing another.

## Before pushing

- View the page rendered, not just the source, at desktop and phone widths.
- Confirm the tracker line is present.
- Confirm no media file got staged by accident: `git status` before `git add`.

## See also

`webpage_building_guide.txt` holds the owner's own notes, including the tracker
snippet above and two loading-screen approaches (an inline Lottie player, and a
redirect-to-loading-screen script). Read it before building a new page.
