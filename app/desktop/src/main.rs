// Desktop shell for the SigmaTau explorer.
//
// The window shows the static pages in app/web/ui and answers the same
// requests serve.py does, as Tauri commands instead of HTTP:
//   health()            -> { ok, core }
//   analyze(id, req)    -> the core's JSON (request shape: app/web/README.md)
//   cancel(id)          -> stop a running analysis
//   save_file(name, data) -> save dialog, then write the bytes
//   live_start(args)    -> run `sigmatau live`, returns its id
//   live_poll(id, since) -> { lines, next, running, error }, its JSON lines from `since`
//   live_stop()         -> stop the live measurement
// The compiled sigmatau binary (cli/, a juliac bundle) ships in the app's
// resources as core/.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use serde_json::{json, Map, Value};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::process::Child;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use tauri::{AppHandle, Manager, State};
use tauri_plugin_dialog::DialogExt;

const TIMEOUT: Duration = Duration::from_secs(600);
const CORE_EXE: &str = if cfg!(windows) { "sigmatau.exe" } else { "sigmatau" };

#[derive(Default)]
struct Runs {
    cancelled: Mutex<HashSet<u64>>,
}

static WORK_SEQ: AtomicU64 = AtomicU64::new(0);

/// The live measurement, one at a time: its process and the JSON lines it has
/// written so far, filled by a reader thread.
struct Live {
    id: u64,
    child: Child,
    lines: Arc<Mutex<Vec<String>>>,
    stderr: Arc<Mutex<String>>,
    reader: std::thread::JoinHandle<()>,
}

#[derive(Default)]
struct LiveSlot(Mutex<Option<Live>>);

fn stop_live(slot: &LiveSlot) {
    if let Some(mut live) = slot.0.lock().unwrap().take() {
        let _ = live.child.kill();
        let _ = live.child.wait();
    }
}

/// The core to run, in order: $SIGMATAU_CLI, the bundled core, then (debug
/// builds) the juliac bundle in the source tree.
fn core_path(app: &AppHandle) -> Option<PathBuf> {
    if let Ok(p) = std::env::var("SIGMATAU_CLI") {
        return Some(PathBuf::from(p));
    }
    if let Ok(dir) = app.path().resource_dir() {
        let p = dir.join("core").join("bin").join(CORE_EXE);
        if p.exists() {
            return Some(p);
        }
    }
    if cfg!(debug_assertions) {
        let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../cli/build/bin").join(CORE_EXE);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

const COMMANDS: [&str; 5] = ["check", "analyze", "spectrum", "dynamic", "drift"];

/// Check a request `{args: [command, flags...], file: {name, content} | null}`
/// (the same rules as serve.py) and return the sigmatau arguments, with the
/// record saved under `workdir`.
fn request_to_args(req: &Map<String, Value>, workdir: &Path) -> Result<Vec<String>, String> {
    let list = req.get("args").and_then(Value::as_array).ok_or("args must be a list of strings")?;
    let mut args = Vec::with_capacity(list.len() + 1);
    for v in list {
        args.push(v.as_str().ok_or("args must be a list of strings")?.to_owned());
    }
    let command = args.first().ok_or("args is empty")?;
    if !COMMANDS.contains(&command.as_str()) {
        return Err(format!("unknown command {command:?}"));
    }
    if args.iter().any(|a| a == "--out" || a.starts_with("--out=")) {
        return Err("--out is not allowed here".into());
    }
    if let Some(file) = req.get("file").and_then(Value::as_object) {
        let name = file.get("name").and_then(Value::as_str).unwrap_or("record.txt");
        let name = Path::new(name).file_name().and_then(|s| s.to_str()).unwrap_or("record.txt");
        let path = workdir.join(name);
        let content = file.get("content").and_then(Value::as_str).unwrap_or("");
        std::fs::write(&path, content).map_err(|e| format!("could not write the record: {e}"))?;
        args.insert(1, path.to_string_lossy().into_owned());
    }
    Ok(args)
}

/// Run the core, polling so a cancel or the timeout can stop it.
fn run_core(core: &Path, args: &[String], id: u64, runs: &Runs) -> Result<Value, String> {
    let mut cmd = Command::new(core);
    // The compiled core runs on one thread unless told otherwise; the
    // total deviations and Thêo1 split their work across all cores.
    if std::env::var_os("JULIA_NUM_THREADS").is_none() {
        cmd.env("JULIA_NUM_THREADS", "auto");
    }
    let mut child = cmd
        .args(args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("could not start the analysis core: {e}"))?;
    let mut out_pipe = child.stdout.take().unwrap();
    let mut err_pipe = child.stderr.take().unwrap();
    let out_reader = std::thread::spawn(move || {
        let mut b = Vec::new();
        let _ = out_pipe.read_to_end(&mut b);
        b
    });
    let err_reader = std::thread::spawn(move || {
        let mut b = Vec::new();
        let _ = err_pipe.read_to_end(&mut b);
        b
    });
    let started = Instant::now();
    let status = loop {
        if let Some(st) = child.try_wait().map_err(|e| e.to_string())? {
            break st;
        }
        if runs.cancelled.lock().unwrap().remove(&id) {
            let _ = child.kill();
            let _ = child.wait();
            return Err("cancelled".into());
        }
        if started.elapsed() > TIMEOUT {
            let _ = child.kill();
            let _ = child.wait();
            return Err("analysis timed out".into());
        }
        std::thread::sleep(Duration::from_millis(5));
    };
    let stdout = out_reader.join().unwrap_or_default();
    let stderr = err_reader.join().unwrap_or_default();
    let text = String::from_utf8_lossy(&stdout);
    let text = text.trim();
    if text.starts_with('{') {
        return serde_json::from_str(text).map_err(|e| format!("the core's output is not JSON: {e}"));
    }
    let err = String::from_utf8_lossy(&stderr);
    Err(err
        .lines()
        .rev()
        .find(|l| !l.trim().is_empty())
        .map(str::to_owned)
        .unwrap_or_else(|| format!("core exited with {status}")))
}

#[tauri::command]
fn health(app: AppHandle) -> Value {
    match core_path(&app) {
        Some(p) => json!({ "ok": true, "core": p.to_string_lossy() }),
        None => json!({ "ok": false, "core": Value::Null }),
    }
}

#[tauri::command]
async fn analyze(app: AppHandle, id: u64, req: Value) -> Result<Value, String> {
    let core = core_path(&app).ok_or("no analysis core: the app bundle is missing core/")?;
    let req = req.as_object().cloned().ok_or("the request must be a JSON object")?;
    tauri::async_runtime::spawn_blocking(move || {
        let runs = app.state::<Runs>();
        let seq = WORK_SEQ.fetch_add(1, Ordering::Relaxed);
        let work = std::env::temp_dir().join(format!("sigmatau-desktop-{}-{seq}", std::process::id()));
        std::fs::create_dir_all(&work).map_err(|e| e.to_string())?;
        let started = Instant::now();
        let result = request_to_args(&req, &work).and_then(|args| {
            let out = run_core(&core, &args, id, &runs);
            if std::env::var_os("SIGMATAU_DESKTOP_LOG").is_some() {
                let what = match &out { Ok(_) => "ok".to_owned(), Err(e) => format!("error: {e}") };
                eprintln!("analyze #{id} {:?}: {what} in {} ms", args, started.elapsed().as_millis());
            }
            out
        });
        let _ = std::fs::remove_dir_all(&work);
        result
    })
    .await
    .map_err(|e| e.to_string())?
}

#[tauri::command]
fn cancel(id: u64, runs: State<'_, Runs>) {
    runs.cancelled.lock().unwrap().insert(id);
}

#[tauri::command]
async fn save_file(app: AppHandle, name: String, data: Vec<u8>) -> Result<bool, String> {
    let Some(target) = app.dialog().file().set_file_name(&name).blocking_save_file() else {
        return Ok(false);
    };
    let path = target.into_path().map_err(|e| e.to_string())?;
    std::fs::write(&path, data).map_err(|e| format!("could not save {}: {e}", path.display()))?;
    Ok(true)
}

#[tauri::command]
fn live_start(app: AppHandle, args: Vec<String>, slot: State<'_, LiveSlot>) -> Result<u64, String> {
    let core = core_path(&app).ok_or("no analysis core: the app bundle is missing core/")?;
    if args.first().map(String::as_str) != Some("live") {
        return Err("a live request starts with the live command".into());
    }
    if args.iter().any(|a| a == "--out" || a.starts_with("--out=")) {
        return Err("--out is not allowed here".into());
    }
    stop_live(&slot);
    let mut child = Command::new(core)
        .args(&args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("could not start the analysis core: {e}"))?;
    let lines = Arc::new(Mutex::new(Vec::new()));
    let stderr = Arc::new(Mutex::new(String::new()));
    let (out, mut err) = (child.stdout.take().unwrap(), child.stderr.take().unwrap());
    let (l, e) = (lines.clone(), stderr.clone());
    let reader = std::thread::spawn(move || {
        for line in BufReader::new(out).lines().map_while(Result::ok) {
            l.lock().unwrap().push(line);
        }
        let _ = err.read_to_string(&mut e.lock().unwrap());
    });
    let id = WORK_SEQ.fetch_add(1, Ordering::Relaxed);
    *slot.0.lock().unwrap() = Some(Live { id, child, lines, stderr, reader });
    Ok(id)
}

#[tauri::command]
fn live_poll(id: u64, since: usize, slot: State<'_, LiveSlot>) -> Result<Value, String> {
    let mut guard = slot.0.lock().unwrap();
    let live = guard.as_mut().filter(|l| l.id == id).ok_or("no such live measurement")?;
    let running = !live.reader.is_finished();
    let lines: Vec<String> = live.lines.lock().unwrap().get(since..).unwrap_or_default().to_vec();
    let mut error = Value::Null;
    if !running {
        if let Ok(Some(st)) = live.child.try_wait() {
            if !st.success() {
                let text = live.stderr.lock().unwrap();
                error = json!(text.lines().rev().find(|l| !l.trim().is_empty()).map(str::to_owned)
                    .unwrap_or_else(|| format!("core exited with {st}")));
            }
        }
    }
    Ok(json!({ "next": since + lines.len(), "lines": lines, "running": running, "error": error }))
}

#[tauri::command]
fn live_stop(slot: State<'_, LiveSlot>) {
    stop_live(&slot);
}

fn main() {
    let launched = Instant::now();
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(Runs::default())
        .manage(LiveSlot::default())
        .invoke_handler(tauri::generate_handler![health, analyze, cancel, save_file, live_start, live_poll, live_stop])
        .on_page_load(move |_, payload| {
            if std::env::var_os("SIGMATAU_DESKTOP_LOG").is_some() {
                eprintln!("page {:?} at {} ms", payload.event(), launched.elapsed().as_millis());
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running the SigmaTau explorer");
}
