/**
 * Types shared between UI and Worker
 */

// Message types from Main → Worker
export interface RunRequest {
  runId: number;
  /** URL to fetch the WASM from (worker will cache it in Cache Storage). */
  wasmUrl?: string;
  /** Raw WASM bytes (alternative to wasmUrl). */
  wasmBytes?: ArrayBuffer;
  mountPaths: string[];
  fileMappings: FileMapping[];
  cmdArgs: string[];
  wasiArgs: string[];
}

export interface FileMapping {
  path: string;
  bytes: ArrayBuffer;
}

// Message types from Worker → Main
export type WorkerMessage =
  | { type: "stdout"; line: string }
  | { type: "stderr"; line: string }
  | { type: "loading"; runId: number; progress: number; message: string }
  | { type: "done"; runId: number; exitCode: number; results: ResultEntry[] }
  | { type: "error"; runId: number; message: string };

// Result entries returned from worker
export interface ResultEntry {
  path: string;
  bytes: Uint8Array;
}

// UI-side file mapping (before transfer to worker)
export interface UIFileMapping {
  file: File;
  path: string;
}
