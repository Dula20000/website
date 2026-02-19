const lastVisit = localStorage.getItem("lastVisit");
const now = Date.now();
const oneDay = 24 * 60 * 60 * 1000; // adjust this threshold

if (!lastVisit || (now - lastVisit) > oneDay) {
  localStorage.setItem("lastVisit", now);
  const returnTo = window.location.pathname;
  window.location.replace("/loadingscreen.html?from=" + encodeURIComponent(returnTo));
}
