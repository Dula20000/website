/** Canvas QR rendering, shared by the beam stream and the pairing code. */

/** Modules of white margin. Scanners need a quiet zone to find the symbol. */
const DEFAULT_QUIET = 3;

/**
 * Draw a QR code filling `canvas`, snapped to whole pixels per module.
 *
 * Integer module size matters more than it sounds: a fractional cell size leaves
 * grey half-pixels on every edge, and against a camera's own sampling grid that
 * is the difference between an instant lock and a code that never reads.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {Uint8Array|string} payload bytes (binary mode) or text
 * @param {{ecc?:'L'|'M'|'Q'|'H', quiet?:number, dark?:string, light?:string}} opts
 * @returns {{version:number, modules:number, cell:number}}
 */
export function renderQR(canvas, payload, opts = {}) {
  const { ecc = 'L', quiet = DEFAULT_QUIET, dark = '#000000', light = '#ffffff' } = opts;

  const segments = typeof payload === 'string' ? payload : [{ data: payload, mode: 'byte' }];
  const qr = QRCodeLib.create(segments, { errorCorrectionLevel: ecc });

  const size = qr.modules.size;
  const data = qr.modules.data;
  const ctx = canvas.getContext('2d', { alpha: false });

  const span = Math.min(canvas.width, canvas.height);
  const cell = Math.max(1, Math.floor(span / (size + quiet * 2)));
  const drawn = cell * (size + quiet * 2);
  const padX = Math.floor((canvas.width - drawn) / 2);
  const padY = Math.floor((canvas.height - drawn) / 2);

  ctx.fillStyle = light;
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = dark;
  for (let y = 0; y < size; y++) {
    const row = y * size;
    for (let x = 0; x < size; x++) {
      if (data[row + x]) {
        ctx.fillRect(padX + (x + quiet) * cell, padY + (y + quiet) * cell, cell, cell);
      }
    }
  }
  return { version: qr.version, modules: size, cell };
}
