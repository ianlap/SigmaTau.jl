// Session files: the clocks (with their file contents), import settings and
// analysis settings in one JSON file, so an analysis can be reopened later or
// handed to someone else. Results are not stored; opening a session reruns.
import { state, newClock } from "./state.js";

const FORMAT = "sigmatau-session";
const VERSION = 1;

// Settings saved as they are.
const SETTINGS = ["devs", "grid", "ci", "confidence", "band", "noiseId", "bias", "axis", "spec", "dyn"];
// Per-clock fields that are results or bookkeeping, not inputs.
const TRANSIENT = new Set(["id", "slot", "record", "result", "dynamic", "error", "loading", "source"]);

export function sessionText() {
  const sources = [];
  const sourceIndex = (s) => {
    let k = sources.findIndex((o) => o.name === s.name && o.content === s.content);
    if (k < 0) k = sources.push({ name: s.name, content: s.content, sha256: s.sha256 ?? "" }) - 1;
    return k;
  };
  const clocks = state.clocks.map((c) => {
    const out = {};
    for (const [k, v] of Object.entries(c)) if (!TRANSIENT.has(k)) out[k] = v;
    out.source = c.source.type === "file" ? { type: "file", file: sourceIndex(c.source) } : { ...c.source };
    out.reference = c.id === state.reference;
    out.active = c.id === state.activeClock;
    return out;
  });
  const settings = Object.fromEntries(SETTINGS.map((k) => [k, state[k]]));
  return JSON.stringify({ format: FORMAT, version: VERSION, saved: new Date().toISOString(), settings, sources, clocks }, null, 1);
}

export const isSession = (text) => /^\s*\{\s*"format"\s*:\s*"sigmatau-session"/.test(text);

// Read and check a session file; throws with a message for the status bar.
export function parseSession(text) {
  let doc;
  try { doc = JSON.parse(text); } catch { throw new Error("this is not a session file (it is not JSON)"); }
  if (doc?.format !== FORMAT) throw new Error("this is not a SigmaTau session file");
  if (doc.version > VERSION) throw new Error(`the session is version ${doc.version}; this explorer reads ${VERSION}`);
  for (const c of doc.clocks ?? []) {
    const f = c.source?.type === "file" ? doc.sources?.[c.source.file] : null;
    if (!c.source || (c.source.type === "file" && typeof f?.content !== "string")) throw new Error("the session is missing a clock's file");
  }
  return doc;
}

// Replace the state's clocks and settings with a session's. The caller then
// reads each clock's record and reruns.
export function applySession(s, doc) {
  for (const k of SETTINGS) if (k in (doc.settings ?? {})) s[k] = doc.settings[k];
  s.clocks = [];
  s.reference = null;
  s.activeClock = null;
  for (const { reference, active, source, ...fields } of doc.clocks ?? []) {
    const src = source.type === "file" ? { type: "file", ...doc.sources[source.file] } : source;
    const c = newClock({ ...fields, source: src });
    s.clocks.push(c);
    if (reference) s.reference = c.id;
    if (active) s.activeClock = c.id;
  }
  s.tableDev = null;
  s.compareTau = null;
  s.dataMode = s.clocks.some((c) => c.source.type === "file") ? "file" : "synth";
  return s.clocks;
}
