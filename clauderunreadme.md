# Clawd Run

A pixelated homage to Chrome's offline dino game, starring **Clawd** — the orange
Claude mascot — galloping on four legs across a Claude-creme desert.

Open `index.html` in any browser. No build step, no dependencies, works offline.

## How to play

- **SPACE / ↑ / W / tap** — jump (hold for a higher hop)
- **↓** — fast-fall while airborne
- **P** — pause

## The twist

- Obstacles are **CAPTCHAs**: checkboxes, distorted-text boxes, and image-select grids. Hit one and you fail the CAPTCHA.
- Coins are **diamond-shaped tokens**. The longer you survive, the more each coin is worth (+1 per 15 seconds) — the coins change colour as their value tier rises: gold → emerald → sapphire → amethyst.
- Tokens accumulate in your **token bank** across runs (saved in localStorage, along with your high score).

All art is generated at runtime from pixel maps in `game.js` — no image files.
