// Host integration in one place, so the pages work the same in a browser
// and in the desktop shell (app/desktop):
//   openFile()  -> { name, text, sha256 } or null   (a file input in both)
//   saveFile(name, blob)                      (browser download, or the shell's save dialog)
//   webLinks(onError)                         (links to the web, from the shell)
// The core itself is reached through core.js.

const tauri = window.__TAURI__?.core;

// A link to the web (target="_blank") opens a tab in a browser. The window of
// the shell opens none: there the link goes to the system browser.
export function webLinks(onError) {
  if (!tauri) return;
  document.addEventListener("click", (e) => {
    const a = e.target.closest?.('a[target="_blank"]');
    if (!a) return;
    e.preventDefault();
    tauri.invoke("open_url", { url: a.href }).catch((err) => onError(String(err)));
  });
}

// A picked or dropped file: its text, and the SHA-256 of its bytes for reports
// (empty where Web Crypto is unavailable).
export async function readFile(f) {
  let bytes = await f.arrayBuffer();
  let name = f.name;
  // A gzip file (the IGS clock products are .gz) is unpacked here: the core reads text.
  if (new Uint8Array(bytes, 0, 2).join() === "31,139") {
    const ds = new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip")));
    bytes = await ds.arrayBuffer();
    name = name.replace(/\.gz$/i, "");
  }
  let sha256 = "";
  if (crypto?.subtle) {
    const d = await crypto.subtle.digest("SHA-256", bytes);
    sha256 = [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
  }
  return { name, text: new TextDecoder().decode(bytes), sha256 };
}

export function openFile(accept = ".txt,.dat,.csv,.tsv,.phd,.frd,.frq,.tim,.log,.clk,.gz,text/*") {
  return new Promise((resolve) => {
    const input = document.createElement("input");
    input.type = "file";
    input.accept = accept;
    input.addEventListener("change", async () => {
      const f = input.files?.[0];
      resolve(f ? await readFile(f) : null);
    });
    input.addEventListener("cancel", () => resolve(null));
    input.click();
  });
}

export async function saveFile(name, blob) {
  if (tauri) {
    const data = Array.from(new Uint8Array(await blob.arrayBuffer()));
    return tauri.invoke("save_file", { name, data });
  }
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
