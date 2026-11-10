// Light or dark. The choice (system, light or dark) is kept in localStorage;
// system follows the OS setting. A plain script in <head>, so the page never
// paints in the wrong theme first. Plots redraw on the "themechange" event.
(() => {
  const KEY = "sigmatau-theme";
  const NEXT = { system: "light", light: "dark", dark: "system" };
  const root = document.documentElement;
  const osDark = matchMedia("(prefers-color-scheme: dark)");
  let choice = "system";
  try { if (NEXT[localStorage.getItem(KEY)]) choice = localStorage.getItem(KEY); } catch {}

  function label() {
    const btn = document.getElementById("theme-btn");
    if (!btn) return;
    btn.dataset.choice = choice;
    btn.title = `Theme: ${choice === "system" ? `system (${root.dataset.theme})` : choice}; click for ${NEXT[choice]}`;
    btn.setAttribute("aria-label", btn.title);
  }

  function apply() {
    const theme = choice === "system" ? (osDark.matches ? "dark" : "light") : choice;
    const changed = root.dataset.theme !== theme;
    root.dataset.theme = theme;
    label();
    if (changed) window.dispatchEvent(new Event("themechange"));
  }

  osDark.addEventListener("change", apply);
  document.addEventListener("DOMContentLoaded", () => {
    label();
    document.getElementById("theme-btn").addEventListener("click", () => {
      choice = NEXT[choice];
      try { localStorage.setItem(KEY, choice); } catch {}
      apply();
    });
  });
  apply();
})();
