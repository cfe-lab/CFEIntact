/**
 * Web Worker entry point for WASI execution
 * Runs container2wasm in an isolated context to keep the UI responsive
 */
import {
  WASI,
  File as WasiFile,
  OpenFile,
  PreopenDirectory,
  ConsoleStdout,
  type Fd,
  type Inode,
} from "@bjorn3/browser_wasi_shim";
import type { RunRequest, WorkerMessage, ResultEntry } from "../shared/types";
import { ensurePath, collectFilesRecursive } from "../shared/fs";

const WASM_CACHE_NAME = "cfeintact-wasm-v1";

/**
 * Fetch WASM bytes from a URL, caching in Cache Storage so subsequent
 * page loads reuse the cached copy instead of re-downloading.
 * Posts `loading` progress messages on the way.
 */
async function fetchWasmWithCache(url: string, runId: number): Promise<Uint8Array> {
  const cache = await caches.open(WASM_CACHE_NAME);
  const cached = await cache.match(url);

  if (cached) {
    const msg: WorkerMessage = { type: "loading", runId, progress: -1, message: "Loading VM image from local cache…" };
    self.postMessage(msg);
    const buf = await cached.arrayBuffer();
    return new Uint8Array(buf);
  }

  // Not cached yet – stream-download with progress reporting.
  console.log(`[CFEIntact worker] Fetching WASM from: ${url}`);
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Failed to fetch WASM: HTTP ${response.status}`);
  if (!response.body) throw new Error("Response body is empty");

  const contentLength = response.headers.get("content-length");
  const total = contentLength ? parseInt(contentLength, 10) : 0;
  const reader = response.body.getReader();
  const chunks: Uint8Array[] = [];
  let received = 0;

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    const progress = total > 0 ? received / total : -1;
    const receivedMB = (received / 1_048_576).toFixed(0);
    const totalMB = total > 0 ? (total / 1_048_576).toFixed(0) : "?";
    const msg: WorkerMessage = {
      type: "loading",
      runId,
      progress,
      message: `Downloading VM image… ${receivedMB} / ${totalMB} MB`,
    };
    self.postMessage(msg);
  }

  // Assemble the full buffer.
  const combined = new Uint8Array(received);
  let offset = 0;
  for (const chunk of chunks) {
    combined.set(chunk, offset);
    offset += chunk.length;
  }

  // Persist in Cache Storage for future visits (best-effort).
  try {
    const cacheResponse = new Response(combined, {
      headers: { "content-type": "application/wasm", "content-length": String(received) },
    });
    await cache.put(url, cacheResponse);
  } catch {
    // Non-fatal – just means it won't be cached.
  }

  return combined;
}

// Type predicate for validating WASI exports
function isWasiStartFunction(fn: unknown): fn is () => unknown {
  return typeof fn === "function";
}

self.onmessage = async (ev: MessageEvent<RunRequest>) => {
  const { runId, wasmUrl, wasmBytes: wasmBytesBuffer, fileMappings, mountPaths, cmdArgs, wasiArgs } = ev.data;

  try {
    // Resolve WASM bytes – either fetch from URL (with caching) or use provided buffer.
    let wasmBytes: Uint8Array;
    if (wasmUrl) {
      wasmBytes = await fetchWasmWithCache(wasmUrl, runId);
    } else if (wasmBytesBuffer) {
      wasmBytes = new Uint8Array(wasmBytesBuffer);
    } else {
      throw new Error("No WASM source provided (set wasmUrl or wasmBytes)");
    }

    // Signal: compiling the WASM module (can take minutes on first run)
    {
      const msg: WorkerMessage = {
        type: "loading",
        runId,
        progress: -1,
        message: "Compiling VM (this can take a minute on first run)…",
      };
      self.postMessage(msg);
    }

    // Build virtual filesystem with all mount points and input files
    const fsRoot = new Map();
    const mountedDirs = new Map<string, Map<string, Inode>>();

    // Create all mount points as writable directories
    for (const mountPath of mountPaths) {
      const dirMap = ensurePath(fsRoot, mountPath);
      mountedDirs.set(mountPath, dirMap);
    }

    // Place input files at their specified paths
    for (const { path, bytes } of fileMappings) {
      const parts = path.split("/").filter((p) => p);
      if (parts.length === 0) continue;

      const filename = parts[parts.length - 1];
      const dirPath = "/" + parts.slice(0, -1).join("/");
      const dirMap = ensurePath(fsRoot, dirPath);

      dirMap.set(filename, new WasiFile(new Uint8Array(bytes)));
    }

    const fds: Fd[] = [
      new OpenFile(new WasiFile([])), // stdin
      ConsoleStdout.lineBuffered((line: string) => {
        const msg: WorkerMessage = { type: "stdout", line };
        self.postMessage(msg);
      }), // stdout
      ConsoleStdout.lineBuffered((line: string) => {
        const msg: WorkerMessage = { type: "stderr", line };
        self.postMessage(msg);
      }), // stderr
    ];

    // Add each mount point as a separate PreopenDirectory
    for (const [mountPath, dirMap] of mountedDirs) {
      fds.push(new PreopenDirectory(mountPath, dirMap));
    }

    const wasi = new WASI(
      cmdArgs.length > 0
        ? ["app.wasm", "--no-stdin", ...wasiArgs, "--", ...cmdArgs]
        : ["app.wasm", "--no-stdin", ...wasiArgs],
      [],
      fds
    );

    // Ensure we pass an ArrayBuffer (not Uint8Array<ArrayBufferLike>) so TypeScript
    // resolves the BufferSource overload of WebAssembly.instantiate correctly.
    const wasmBuffer: ArrayBuffer =
      wasmBytes.buffer instanceof ArrayBuffer
        ? wasmBytes.buffer
        : (wasmBytes.buffer.slice(0) as unknown as ArrayBuffer);
    const { instance } = await WebAssembly.instantiate(wasmBuffer, {
      wasi_snapshot_preview1: wasi.wasiImport,
    });

    // Extract and validate exports with proper type guards
    const { memory, _start } = instance.exports;
    
    if (!(memory instanceof WebAssembly.Memory)) {
      throw new Error("Expected WebAssembly.Memory in exports");
    }
    if (!isWasiStartFunction(_start)) {
      throw new Error("Expected _start function in exports");
    }

    // TypeScript now knows _start is () => unknown via type predicate
    const wasiInstance: { exports: { memory: WebAssembly.Memory; _start: () => unknown } } = {
      exports: { memory, _start },
    };

    // Signal: VM is up, CFEIntact is starting
    {
      const msg: WorkerMessage = { type: "loading", runId, progress: -1, message: "Running CFEIntact…" };
      self.postMessage(msg);
    }

    const exitCode = wasi.start(wasiInstance);

    // Collect all files from all mount points recursively
    const results: ResultEntry[] = [];
    const transfers: Transferable[] = [];

    for (const [mountPath, dirMap] of mountedDirs) {
      // Collect all files recursively under this mount
      const files = collectFilesRecursive(dirMap, mountPath);
      
      for (const { path, bytes } of files) {
        results.push({
          path,
          bytes,
        });
        // Transfer the underlying buffer for performance
        // bytes.buffer is ArrayBufferLike which includes ArrayBuffer and SharedArrayBuffer
        if (bytes.buffer instanceof ArrayBuffer) {
          transfers.push(bytes.buffer);
        }
      }
    }

    const msg: WorkerMessage = { type: "done", runId, exitCode, results };
    self.postMessage(msg, transfers);
  } catch (e: unknown) {
    const message = e instanceof Error ? e.message : String(e);
    const msg: WorkerMessage = {
      type: "error",
      runId,
      message,
    };
    self.postMessage(msg);
  }
};
