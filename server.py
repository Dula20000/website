#!/usr/bin/env python3
"""Local quote server for the gold hedge calculator.

Fetches server-side; browsers can't hit these directly (no CORS headers):
  spot     : api.gold-api.com          -> XAU/USD spot
  futures  : TradingView COMEX:1OZV2026 -> Oct 2026 1-Ounce Gold, 10-min delayed
             (falls back to Yahoo 1OZV26.CMX, which 429s under heavy polling)
             also returns COMEX:GCV2026 as a sanity cross-check

Everything is cached, and a failed refresh serves the last good price marked
stale rather than dropping the value.
"""

import http.cookiejar
import json
import os
import time
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

PORT = int(os.environ.get("PORT", 8787))
HERE = os.path.dirname(os.path.abspath(__file__))
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")

FUT_TTL = 25       # seconds between live 1OZV26 fetches
SPOT_TTL = 20

_opener = urllib.request.build_opener(
    urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
_last = {}   # key -> {"value": dict, "at": epoch}


def _get_json(url, timeout=8):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with _opener.open(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def _iso(epoch):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(epoch))


def yahoo_quote(symbol):
    """(price, epoch) for a Yahoo symbol; tries both query hosts on failure."""
    last_err = None
    for host in ("query1", "query2"):
        try:
            url = (f"https://{host}.finance.yahoo.com/v8/finance/chart/"
                   f"{symbol}?interval=1m&range=1d")
            meta = _get_json(url)["chart"]["result"][0]["meta"]
            price = meta.get("regularMarketPrice")
            if price is None:
                raise ValueError("no regularMarketPrice")
            return float(price), int(meta.get("regularMarketTime") or time.time())
        except Exception as e:
            last_err = e
            time.sleep(0.4)
    raise last_err


def spot_price():
    d = _get_json("https://api.gold-api.com/price/XAU")
    return {"price": float(d["price"]), "source": "gold-api.com",
            "asOf": d.get("updatedAt") or _iso(time.time())}


def tv_futures():
    """Primary: TradingView returns 1OZV26 and GCV26 in one call, 10-min delayed."""
    body = json.dumps({
        "symbols": {"tickers": ["COMEX:1OZV2026", "COMEX:GCV2026"]},
        "columns": ["close", "change", "update_mode"],
    }).encode("utf-8")
    req = urllib.request.Request(
        "https://scanner.tradingview.com/futures/scan", data=body,
        headers={"User-Agent": UA, "Content-Type": "application/json"})
    with _opener.open(req, timeout=8) as r:
        rows = json.loads(r.read().decode("utf-8"))["data"]

    by_sym = {row["s"]: row["d"] for row in rows}
    one = by_sym.get("COMEX:1OZV2026")
    if not one or one[0] is None:
        raise ValueError("no 1OZV2026 price")

    mode = one[2] or ""
    delay = 600 if "600" in mode else (900 if "900" in mode else 0)
    cross = by_sym.get("COMEX:GCV2026")
    return {
        "price": float(one[0]),
        "changePct": one[1],
        "symbol": "1OZV26",
        "source": "tradingview (COMEX)",
        "delaySec": delay,
        "asOf": _iso(time.time() - delay),
        "cross": float(cross[0]) if cross and cross[0] is not None else None,
    }


def yahoo_futures():
    """Fallback if TradingView fails. Yahoo 429s under heavy polling."""
    px, ts = yahoo_quote("1OZV26.CMX")
    out = {"price": px, "symbol": "1OZV26", "source": "yahoo (fallback)",
           "delaySec": 0, "asOf": _iso(ts), "cross": None}
    try:
        out["cross"] = yahoo_quote("GCV26.CMX")[0]
    except Exception:
        pass
    return out


def futures_quote():
    try:
        return tv_futures()
    except Exception as tv_err:
        try:
            return yahoo_futures()
        except Exception as y_err:
            raise RuntimeError(f"tradingview: {tv_err}; yahoo: {y_err}")


def cached(key, ttl, fetch, errors):
    """Fetch through cache. On failure, serve the last good value as stale."""
    now = time.time()
    prev = _last.get(key)
    if prev and now - prev["at"] < ttl:
        return dict(prev["value"], stale=False, ageSec=int(now - prev["at"]))
    try:
        value = fetch()
        _last[key] = {"value": value, "at": now}
        return dict(value, stale=False, ageSec=0)
    except Exception as e:
        errors.append(f"{key}: {e}")
        if prev:
            return dict(prev["value"], stale=True, ageSec=int(now - prev["at"]))
        return None


def collect():
    errors = []
    out = {
        "spot": cached("spot", SPOT_TTL, spot_price, errors),
        "futures": cached("futures", FUT_TTL, futures_quote, errors),
        "fetchedAt": _iso(time.time()),
    }
    out["errors"] = errors
    return out


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype):
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        path = self.path.split("?")[0]
        if path == "/api/quotes":
            self._send(200, json.dumps(collect()), "application/json")
            return
        if path in ("/", "/index.html"):
            with open(os.path.join(HERE, "index.html"), "rb") as f:
                self._send(200, f.read(), "text/html; charset=utf-8")
            return
        self._send(404, "not found", "text/plain")

    def log_message(self, fmt, *args):
        return  # quiet: the page polls every 30s


if __name__ == "__main__":
    print(f"gold hedge calculator -> http://localhost:{PORT}")
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
