# Visual Transfer

Move files between devices by pointing a camera at a screen. No accounts, no cables,
no file ever stored on a server.

Two transports share one interface — you scan a code either way:

| | **Link** | **Beam** |
|---|---|---|
| What the code holds | a one-time pairing address | the file data itself |
| Carries bytes over | direct WebRTC peer-to-peer | the camera, optically |
| Speed | network speed (MB/s) | ~1–10 KB/s |
| Practical size | any | up to a few MB |
| Needs a network | yes (same LAN is enough) | **no** — fully air-gapped |
| Receiver needs this app | no, any camera app opens it | yes, for the scanner |

Both paths move raw bytes and preserve the original filename and media type, so any
format works — photos, video, PDFs, archives, code, 3D models, databases. Nothing is
re-encoded, and every transfer is CRC-32 verified before it is offered for saving.

## Run it locally

```bash
npm install && npm start
```

That serves the app on `http://localhost:8080` **and** on `https://localhost:8443`,
printing your LAN address for the phone:

```
https://10.0.0.70:8443   ← open this on your phone
```

Use the HTTPS address on any other device. Browsers only grant camera access in a
secure context, so a phone loading `http://192.168.x.x` cannot open its camera at
all. A self-signed certificate covering localhost and your current LAN IPs is
generated on first run (needs `openssl`); the phone will warn once — tap
**Advanced → Proceed**.

## Deploy it

It is a static frontend plus a tiny WebSocket relay, so any Node host works. Bind
to the platform's port and let the platform terminate TLS:

```bash
PORT=$PORT node server.js
```

`NO_HTTPS=1` skips local certificate generation (automatic when `RENDER` or
`VERCEL` is set). `/healthz` answers `ok` for health checks. Nothing is persisted,
so it scales horizontally only if both peers land on the same instance — use a
sticky-session or single-instance setup.

### On static hosting (GitHub Pages)

The `public/` directory is a complete static app and can be served by GitHub Pages
directly — **Beam works fully**, since the optical path involves no server at all.

Link mode will not work there. It needs the WebSocket relay to introduce the two
devices, and Pages cannot run Node. The app detects this and says so instead of
failing vaguely. To get Link working, run `server.js` anywhere that executes Node
and use that address.

## How Beam works

The interesting part is that the optical path has **no back-channel**. The sender
cannot be told what arrived, so it never asks.

The file is split into `K` blocks, and each frame carries an XOR of a
pseudo-random subset of them — a [Luby transform](https://en.wikipedia.org/wiki/Luby_transform_code)
fountain code. Which blocks a frame combines is derived from its 32-bit seed, so a
frame is entirely self-describing and the receiver needs no ordering information.

That changes the failure mode completely. A camera watching a screen drops frames
unpredictably — motion blur, autofocus hunting, a hand moving — and a numbered
scheme has to detect each gap and request a resend. A fountain stream has no gaps
to detect: the receiver collects whatever it catches, in any order, and finishes
once it has slightly more than `K` useful frames. Measured overhead is 1–2% even
when 60% of frames are lost, so the receiver is just a camera you point, and the
sender just loops forever.

The first `K` seeds are deliberately systematic (one block each) so a fresh
receiver gets usable data immediately rather than waiting for the peeling decoder
to find a foothold.

Frame density presets trade capacity against how steady a hand the code needs:

| Preset | Payload | QR version |
|---|---|---|
| Gentle | 400 B | 13 |
| Balanced | 800 B | 20 |
| Turbo | 1600 B | 29 |

Sending multiple files rotates through the queue and repeats indefinitely, so
anything missed on one pass lands on the next. The receiver decodes every file it
sees in parallel and ignores ones it has already completed — point the camera and
wait, no coordination.

## How Link works

The QR holds `https://<host>/#j=CODE`. Scanning it with **any** camera app opens
the page, which auto-joins that pairing code, and the two devices negotiate a
WebRTC data channel. The relay forwards only SDP and ICE — file contents never
touch the server. Files stream off disk in chunks sized to the negotiated SCTP
message limit, so a multi-gigabyte file does not have to fit in memory on the
sending side.

Once connected the channel is symmetric: the receiving device can send files back
without re-pairing.

## Tests

```bash
npm test
```

Covers the parts where a bug means silent corruption: CRC-32 against a known
vector, XOR fast paths including unaligned views, packet framing with Unicode
filenames, encoder/decoder seed agreement, fountain round trips at 0–85% frame
loss, late joins, duplicate frames, and a full optical pipeline simulation. The
relay suite boots the real server and drives it with WebSocket clients.

Verified in-browser beyond that: rendered QR frames read back through jsQR
byte-exact, the scan worker's message protocol, and real WebRTC transfers up to
24 MB between two tabs with checksums compared on both ends.

## Layout

```
server.js              static host + WebSocket pairing relay
public/
  index.html           all views
  app.css
  js/
    fountain.js        LT encoder/decoder — the core of Beam
    packet.js          binary frame format
    optical.js         beam sender (render loop) + receiver (camera)
    scanworker.js      jsQR off the main thread
    qrdraw.js          canvas QR rendering
    rtc.js             signalling + peer-to-peer file protocol
    util.js            CRC-32, XOR, formatting
    main.js            views, staging, wiring
  vendor/              bundled qrcode + jsQR (offline-capable)
test/
```

## Limits worth knowing

- Received files are assembled in memory before being offered for saving, so the
  receiving device's RAM caps a single file. The sending side streams and is not
  affected.
- Beam is genuinely slow. At Balanced density and 10 fps the ceiling is about
  8 KB/s. It is for a config file, a key, a photo across an air gap — not a video.
- Link needs a peer-to-peer path. Two devices on the same Wi-Fi are fine; some
  corporate or mobile-carrier networks block it, and there is no TURN relay
  configured to fall back on.
- Browsers may block the automatic download of a received file. Each file also
  appears as a card with an explicit Save button for exactly that reason.
