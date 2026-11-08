// Host integration in one place, so the pages work the same in a browser
// and in the desktop shell (app/desktop):
//   openFile()  -> { name, text } or null   (a file input in both)
//   saveFile(name, blob)                      (browser download, or the shell's save dialog)
// The core itself is reached through core.js.

const tauri = window.__TAURI__?.core;

export function openFile(accept = ".txt,.dat,.csv,.tsv,.phd,.frq,text/*") {
  return new Promise((resolve) => {
    const input = document.createElement("input");
    input.type = "file";
    input.accept = accept;
    input.addEventListener("change", async () => {
      const f = input.files?.[0];
      resolve(f ? { name: f.name, text: await f.text() } : null);
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
