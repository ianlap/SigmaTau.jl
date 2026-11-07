// Host integration in one place, so a desktop shell (Tauri) can replace it
// without touching the pages:
//   openFile()  -> { name, text } or null   (Tauri: dialog plugin + fs read)
//   saveFile(name, blob)                      (Tauri: dialog save + fs write)
// The core itself is reached through api.js.

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

export function saveFile(name, blob) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
