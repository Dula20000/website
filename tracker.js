//this is the tracking code snippet from cloudfair worker and worker Kv. It uses FingerprintJS to log everything needed to fingerprint users. have this in all html files going forward.

const fpPromise = import('https://openfpcdn.io/fingerprintjs/v4')
  .then(FingerprintJS => FingerprintJS.load())

fpPromise
  .then(fp => fp.get())
  .then(result => {
    fetch('https://site-tracker.junkmails20d.workers.dev', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        fingerprintId: result.visitorId,
        page: window.location.pathname
      })
    })
  })

//this is the tracking code snippet from cloudfair worker and worker Kv. It uses FingerprintJS to log everything needed to fingerprint users. have this in all html files going forward.
