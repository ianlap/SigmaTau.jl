// Desktop shell for the SigmaTau explorer.
//
// The window shows the static pages in app/web/ui and answers the same
// requests serve.py does, as Tauri commands instead of HTTP:
//   health()            -> { ok, core }
//   analyze(id, req)    -> the core's JSON (request shape: app/web/README.md)
//   cancel(id)          -> stop a running analysis
//   save_file(name, data) -> save dialog, then write the bytes
// The compiled core (juliac bundle) ships in the app's resources as core/.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use serde_json::{json, Map, Value};
use std::collections::HashSet;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};
use tauri::{AppHandle, Manager, State};
use tauri_plugin_dialog::DialogExt;

const TIMEOUT: Duration = Duration::from_secs(600);
const CORE_EXE: &str = if cfg!(windows) { "sigmatau-cli.exe" } else { "sigmatau-cli" };

#[derive(Default)]
struct Runs {
    cancelled: Mutex<HashSet<u64>>,
}

static WORK_SEQ: AtomicU64 = AtomicU64::new(0);

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
        let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("../web/cli/build/bin").join(CORE_EXE);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

fn present<'a>(req: &'a Map<String, Value>, key: &str) -> Option<&'a Value> {
    req.get(key).filter(|v| !v.is_null())
}

fn as_int(v: &Value) -> Result<i64, String> {
    v.as_i64()
        .or_else(|| v.as_f64().map(|f| f as i64))
        .or_else(|| v.as_str().and_then(|s| s.trim().parse().ok()))
        .ok_or_else(|| format!("expected an integer, got {v}"))
}

fn as_float(v: &Value) -> Result<f64, String> {
    v.as_f64()
        .or_else(|| v.as_str().and_then(|s| s.trim().parse().ok()))
        .ok_or_else(|| format!("expected a number, got {v}"))
}

fn as_text(v: &Value) -> String {
    match v {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// Translate a request into sigmatau-cli flags (the same mapping as serve.py).
fn request_to_args(req: &Map<String, Value>, workdir: &Path) -> Result<Vec<String>, String> {
    let mut args: Vec<String> = Vec::new();
    let src = present(req, "source").and_then(Value::as_object).cloned().unwrap_or_default();
    match src.get("type").and_then(Value::as_str) {
        Some("synth") => {
            let n = src.get("n").map(as_int).transpose()?.unwrap_or(8192);
            let seed = src.get("seed").map(as_int).transpose()?.unwrap_or(42);
            args.extend(["--synth".into(), n.to_string(), "--seed".into(), seed.to_string()]);
        }
        Some("file") => {
            let name = src.get("name").and_then(Value::as_str).unwrap_or("record.txt");
            let name = Path::new(name).file_name().and_then(|s| s.to_str()).unwrap_or("record.txt");
            let path = workdir.join(name);
            let content = src.get("content").and_then(Value::as_str).unwrap_or("");
            std::fs::write(&path, content).map_err(|e| format!("could not write the record: {e}"))?;
            args.extend(["--input".into(), path.to_string_lossy().into_owned()]);
        }
        _ => return Err("source.type must be 'file' or 'synth'".into()),
    }

    for (key, flag) in [("kind", "--kind"), ("detrend", "--detrend"), ("grid", "--grid")] {
        if let Some(v) = present(req, key) {
            args.extend([flag.into(), as_text(v)]);
        }
    }
    for (key, flag) in [("tau0", "--tau0"), ("scaling", "--scaling")] {
        if let Some(v) = present(req, key) {
            args.extend([flag.into(), as_float(v)?.to_string()]);
        }
    }
    for (key, flag) in [
        ("time_col", "--time-col"),
        ("value_col", "--value-col"),
        ("reference_col", "--reference-col"),
        ("header", "--header"),
        ("series", "--series"),
    ] {
        if let Some(v) = present(req, key) {
            args.extend([flag.into(), as_int(v)?.to_string()]);
        }
    }
    if let Some(devs) = present(req, "devs").and_then(Value::as_array) {
        let list: Vec<String> = devs.iter().map(as_text).collect();
        let joined = list.join(",");
        args.extend(["--devs".into(), if joined.is_empty() { "none".into() } else { joined }]);
    }
    let truthy = |key: &str| present(req, key).map_or(false, |v| v.as_bool().unwrap_or(true));
    if truthy("fillgaps") {
        args.push("--fillgaps".into());
    }
    if truthy("remove_outliers") {
        args.push("--remove-outliers".into());
    }
    match req.get("confidence") {
        Some(Value::Null) => args.push("--no-ci".into()),
        Some(v) => args.extend(["--confidence".into(), as_float(v)?.to_string()]),
        None => {}
    }
    if req.get("bias_correction") == Some(&Value::Bool(false)) {
        args.push("--no-bias-correction".into());
    }
    if truthy("spectrum") {
        args.push("--spectrum".into());
    }
    if let Some(v) = present(req, "nperseg") {
        args.extend(["--nperseg".into(), as_int(v)?.to_string()]);
    }
    if let Some(v) = present(req, "dynamic").filter(|v| v.as_bool() != Some(false)) {
        args.extend(["--dynamic".into(), as_text(v)]);
    }
    for (key, flag) in [("window", "--window"), ("step", "--step")] {
        if let Some(v) = present(req, key) {
            args.extend([flag.into(), as_int(v)?.to_string()]);
        }
    }
    Ok(args)
}

/// Run the core, polling so a cancel or the timeout can stop it.
fn run_core(core: &Path, args: &[String], id: u64, runs: &Runs) -> Result<Value, String> {
    let mut child = Command::new(core)
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

fn main() {
    let launched = Instant::now();
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(Runs::default())
        .invoke_handler(tauri::generate_handler![health, analyze, cancel, save_file])
        .on_page_load(move |_, payload| {
            if std::env::var_os("SIGMATAU_DESKTOP_LOG").is_some() {
                eprintln!("page {:?} at {} ms", payload.event(), launched.elapsed().as_millis());
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running the SigmaTau explorer");
}
