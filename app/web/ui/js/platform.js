// Host integration in one place, so the pages work the same in a browser
// and in the desktop shell (app/desktop):
//   openFile()  -> { name, text, sha256 } or null   (a file input in both)
//   saveFile(name, blob)                      (browser download, or the shell's save dialog)
// The core itself is reached through core.js.

const tauri = window.__TAURI__?.core;

// A picked or dropped file: its text, and the SHA-256 of its bytes for reports
// (empty where Web Crypto is unavailable).
export async function readFile(f) {
  const bytes = await f.arrayBuffer();
  let sha256 = "";
  if (crypto?.subtle) {
    const d = await crypto.subtle.digest("SHA-256", bytes);
    sha256 = [...new Uint8Array(d)].map((b) => b.toString(16).padStart(2, "0")).join("");
  }
  return { name: f.name, text: new TextDecoder().decode(bytes), sha256 };
}

export function openFile(accept = ".txt,.dat,.csv,.tsv,.phd,.frd,.frq,.tim,.log,text/*") {
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
