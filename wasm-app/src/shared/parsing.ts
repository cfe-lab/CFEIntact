/**
 * Parsing utilities for command-line arguments and paths
 */

/**
 * A small shell-ish splitter: handles quotes (single/double) and backslash escapes in double quotes.
 */
export function shellSplit(s: string): string[] {
  const out: string[] = [];
  let cur = "";
  let q: "'" | '"' | null = null;

  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if (!q && /\s/.test(ch)) {
      if (cur) {
        out.push(cur);
        cur = "";
      }
      continue;
    }
    if (!q && (ch === "'" || ch === '"')) {
      q = ch;
      continue;
    }
    if (q && ch === q) {
      q = null;
      continue;
    }
    if (q === '"' && ch === "\\" && i + 1 < s.length) {
      const nxt = s[++i];
      cur += nxt;
      continue;
    }
    cur += ch;
  }
  if (cur) out.push(cur);
  return out;
}

/**
 * Validate that a path is absolute (starts with /)
 */
export function isAbsolutePath(path: string): boolean {
  return path.startsWith("/");
}

/**
 * Normalize a path by removing duplicate slashes and trailing slashes
 */
export function normalizePath(path: string): string {
  return path.replace(/\/+/g, "/").replace(/\/$/, "") || "/";
}
