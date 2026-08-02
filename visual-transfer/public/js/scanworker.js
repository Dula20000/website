/**
 * QR decode worker. Keeping jsQR off the main thread is what keeps the camera
 * preview smooth — a 1280px frame decode costs tens of milliseconds, which would
 * otherwise show up as visible stutter while the user is aiming.
 */
/* eslint-env worker */
importScripts('../vendor/jsqr.js');

self.onmessage = (event) => {
  const { buffer, width, height, id } = event.data;
  const pixels = new Uint8ClampedArray(buffer);

  let result = null;
  try {
    // 'dontInvert' — we always render standard dark-on-light codes, and skipping
    // the inverted pass roughly halves decode time.
    result = self.jsQR(pixels, width, height, { inversionAttempts: 'dontInvert' });
  } catch {
    result = null;
  }

  if (result && result.binaryData && result.binaryData.length) {
    const bytes = new Uint8Array(result.binaryData);
    self.postMessage(
      { id, found: true, bytes, text: result.data || '', buffer },
      [bytes.buffer, buffer],
    );
  } else {
    // Hand the pixel buffer back so the main thread can reuse the allocation.
    self.postMessage({ id, found: false, buffer }, [buffer]);
  }
};
