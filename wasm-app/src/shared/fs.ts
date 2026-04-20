/**
 * Filesystem utilities for working with browser_wasi_shim structures
 */
import { Directory, File as WasiFile, type Inode } from "@bjorn3/browser_wasi_shim";

// Type guard for objects with string index signature
function isObjectWithStringKeys(obj: unknown): obj is Record<string, unknown> {
  return typeof obj === "object" && obj !== null;
}

/**
 * Try to find a Map of entries inside a Directory/PreopenDirectory-like object.
 */
export function extractEntriesMap(dirObj: unknown): Map<string, unknown> | null {
  if (!dirObj || typeof dirObj !== "object") return null;
  if (dirObj instanceof Map) return dirObj;

  if (!isObjectWithStringKeys(dirObj)) return null;

  // Common places
  for (const k of ["entries", "contents", "map", "_entries", "_contents", "dir"]) {
    const v = dirObj[k];
    if (v instanceof Map) return v;
  }
  // Otherwise: scan any property for a Map
  for (const v of Object.values(dirObj)) {
    if (v instanceof Map) return v;
  }
  return null;
}

/**
 * Try hard to extract the current byte buffer from a browser_wasi_shim File instance
 * without depending on a single private field name.
 */
export function extractFileBytes(fileObj: unknown): Uint8Array {
  if (!fileObj || typeof fileObj !== "object") return new Uint8Array();
  
  if (!isObjectWithStringKeys(fileObj)) return new Uint8Array();

  const candidates = ["data", "buffer", "contents", "content", "_data", "_buffer", "buf"];
  for (const k of candidates) {
    const v = fileObj[k];
    if (v instanceof Uint8Array) return v;
    if (v instanceof ArrayBuffer) return new Uint8Array(v);
  }
  // Fallback: scan for a Uint8Array-valued property
  for (const v of Object.values(fileObj)) {
    if (v instanceof Uint8Array) return v;
    if (v instanceof ArrayBuffer) return new Uint8Array(v);
  }
  return new Uint8Array();
}

/**
 * Recursively collect all files under a directory-ish node.
 * Returns flattened list of {path, bytes} entries.
 */
export function collectFilesRecursive(node: unknown, prefix = ""): Array<{ path: string; bytes: Uint8Array }> {
  const results: Array<{ path: string; bytes: Uint8Array }> = [];
  const m = extractEntriesMap(node);
  if (!m) return results;

  for (const [name, child] of m.entries()) {
    const path = prefix ? `${prefix}/${name}` : name;

    // Heuristic: if it "looks like" a file (has bytes), treat as file; else treat as directory.
    const bytes = extractFileBytes(child);
    const hasBytes = bytes && bytes.length !== undefined && (bytes.length > 0 || child instanceof WasiFile);

    if (hasBytes) {
      results.push({ path, bytes });
    } else {
      // Might be Directory or something wrapping it
      results.push(...collectFilesRecursive(child, path));
    }
  }
  return results;
}

/**
 * Ensure directory path exists in a virtual filesystem Map.
 * Returns the Map representing the deepest directory.
 */
export function ensurePath(fsRoot: Map<string, Inode>, path: string): Map<string, Inode> {
  const parts = path.split("/").filter((p) => p);
  let current = fsRoot;

  for (const part of parts) {
    if (!current.has(part)) {
      const newDirMap = new Map<string, Inode>();
      current.set(part, new Directory(newDirMap));
      current = newDirMap;
    } else {
      const entry = current.get(part);
      const extractedMap = extractEntriesMap(entry);
      if (extractedMap && isInodeMap(extractedMap)) {
        current = extractedMap;
      } else {
        current = new Map<string, Inode>();
      }
    }
  }
  return current;
}

// Type guard to check if a Map has Inode values
function isInodeMap(map: Map<string, unknown>): map is Map<string, Inode> {
  // If empty, assume it's a valid Inode map
  if (map.size === 0) return true;
  // Check a sample entry
  for (const value of map.values()) {
    return value instanceof WasiFile || value instanceof Directory;
  }
  return true;
}
