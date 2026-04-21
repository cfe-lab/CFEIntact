import { useState, useEffect, useRef, useCallback } from "react";
import {
  Box,
  Container,
  VStack,
  Heading,
  Button,
  Text,
  HStack,
  Badge,
} from "@chakra-ui/react";
import { toaster } from "../ui/toaster";
import type { WorkerMessage, ResultEntry, RunRequest } from "../shared/types";

// ─── Constants ────────────────────────────────────────────────────────────────

// Resolve to an absolute URL on the main thread before sending to the Web Worker.
// A bare relative string like "./cfeintact.wasm" would be resolved against the
// *worker script* URL (e.g. assets/runner-xxx.js), not the page URL, causing a 404.
// document.baseURI is always absolute and already ends with a trailing slash.
const WASM_URL = new URL("cfeintact.wasm", document.baseURI).href;
const INPUT_MOUNT = "/mnt/in";
/**
 * Temporary directory INSIDE the container's native (non-9P) filesystem where
 * CFEIntact writes its outputs. Using the container's internal /tmp avoids the
 * 9P-backed /w path, which triggers EOPNOTSUPP when Python truncates an existing
 * file (e.g. blast.csv on the second subtype iteration).
 */
const CONTAINER_WORK_DIR = "/tmp/work";
/** WASI-preopened output mount where cp deposits files for the host to read. */
const OUTPUT_MOUNT = "/mnt/out";
const INPUT_PATH = `${INPUT_MOUNT}/input.fasta`;

/**
 * Default command shown to the user in the Advanced options field.
 * The app automatically wraps this with /bin/sh and appends the cp step.
 */
const DEFAULT_CMD = `cfeintact check ${INPUT_PATH}`;

// ─── Helpers ──────────────────────────────────────────────────────────────────

function formatBytes(n: number): string {
  if (n < 1024) return `${n} B`;
  if (n < 1_048_576) return `${(n / 1024).toFixed(1)} KB`;
  return `${(n / 1_048_576).toFixed(1)} MB`;
}

function downloadBlob(bytes: Uint8Array, filename: string) {
  // .slice() produces a Uint8Array<ArrayBuffer> which satisfies BlobPart
  const blob = new Blob([bytes.slice()]);
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function basename(path: string): string {
  return path.split("/").filter(Boolean).at(-1) ?? path;
}

// ─── DropZone ─────────────────────────────────────────────────────────────────

interface DropZoneProps {
  file: File | null;
  onFile: (f: File) => void;
  disabled?: boolean;
}

function DropZone({ file, onFile, disabled }: DropZoneProps) {
  const [dragging, setDragging] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setDragging(false);
    if (disabled) return;
    const dropped = e.dataTransfer.files[0];
    if (dropped) onFile(dropped);
  };

  return (
    <Box
      border="2px dashed"
      borderColor={dragging ? "teal.400" : file ? "teal.300" : "gray.300"}
      borderRadius="xl"
      p={8}
      textAlign="center"
      cursor={disabled ? "not-allowed" : "pointer"}
      bg={dragging ? "teal.50" : file ? "teal.50" : "gray.50"}
      transition="all 0.15s"
      _hover={disabled ? {} : { borderColor: "teal.400", bg: "teal.50" }}
      onDragOver={(e) => { e.preventDefault(); if (!disabled) setDragging(true); }}
      onDragLeave={() => setDragging(false)}
      onDrop={handleDrop}
      onClick={() => { if (!disabled) inputRef.current?.click(); }}
    >
      <input
        ref={inputRef}
        type="file"
        accept=".fasta,.fa,.fna,.fas"
        style={{ display: "none" }}
        onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
          const f = e.target.files?.[0];
          if (f) onFile(f);
        }}
        disabled={disabled}
      />
      {file ? (
        <VStack gap={1}>
          <Text fontSize="2xl">🧬</Text>
          <Text fontWeight="semibold" color="teal.700">{file.name}</Text>
          <Text fontSize="sm" color="gray.500">{formatBytes(file.size)}</Text>
          <Text fontSize="xs" color="teal.500" mt={1}>Click or drop to replace</Text>
        </VStack>
      ) : (
        <VStack gap={2}>
          <Text fontSize="3xl">📂</Text>
          <Text fontWeight="semibold" color="gray.600">Drop your FASTA file here</Text>
          <Text fontSize="sm" color="gray.400">or click to browse</Text>
          <Text fontSize="xs" color="gray.400" mt={1}>.fasta / .fa / .fna / .fas</Text>
        </VStack>
      )}
    </Box>
  );
}

// ─── ProgressBar ──────────────────────────────────────────────────────────────

interface ProgressBarProps {
  /** 0–1, or −1 to show an indeterminate animation */
  progress: number;
  label: string;
}

function ProgressBar({ progress, label }: ProgressBarProps) {
  const indeterminate = progress < 0;
  const pct = indeterminate ? 40 : Math.round(progress * 100);

  return (
    <VStack gap={2} align="stretch">
      <HStack justify="space-between">
        <Text fontSize="sm" color="gray.600">{label}</Text>
        {!indeterminate && <Text fontSize="sm" color="gray.500">{pct}%</Text>}
      </HStack>
      <Box h="8px" bg="gray.200" borderRadius="full" overflow="hidden">
        <Box
          h="100%"
          bg="teal.400"
          borderRadius="full"
          style={{ width: `${pct}%` }}
          transition={indeterminate ? "none" : "width 0.3s ease"}
        />
      </Box>
    </VStack>
  );
}

// ─── ResultCard ───────────────────────────────────────────────────────────────

interface ResultCardProps {
  entry: ResultEntry;
}

function ResultCard({ entry }: ResultCardProps) {
  const name = basename(entry.path);
  const isCSV = name.toLowerCase().endsWith(".csv");
  const [preview, setPreview] = useState<string | null>(null);

  return (
    <Box
      bg="white"
      borderRadius="lg"
      borderWidth={1}
      borderColor={isCSV ? "teal.200" : "gray.200"}
      p={4}
      shadow="sm"
    >
      <HStack justify="space-between" align="start">
        <VStack align="start" gap={0}>
          <HStack gap={2}>
            <Text fontWeight="semibold" fontSize="sm">{name}</Text>
            {isCSV && <Badge colorScheme="teal" size="sm">CSV</Badge>}
          </HStack>
          <Text fontSize="xs" color="gray.500">{formatBytes(entry.bytes.length)}</Text>
        </VStack>
        <HStack gap={2}>
          {isCSV && (
            <Button
              size="sm"
              variant="ghost"
              colorScheme="gray"
              onClick={() =>
                setPreview((p) =>
                  p === null
                    ? new TextDecoder().decode(entry.bytes).slice(0, 3000)
                    : null
                )
              }
            >
              {preview !== null ? "Hide" : "Preview"}
            </Button>
          )}
          <Button
            size="sm"
            colorScheme="teal"
            onClick={() => downloadBlob(entry.bytes, name)}
          >
            Download
          </Button>
        </HStack>
      </HStack>
      {preview !== null && (
        <Box
          mt={3}
          bg="gray.50"
          borderRadius="md"
          p={3}
          maxH="220px"
          overflowY="auto"
          borderWidth={1}
          borderColor="gray.200"
        >
          <pre
            style={{
              margin: 0,
              fontSize: "0.75rem",
              whiteSpace: "pre-wrap",
              wordBreak: "break-all",
              fontFamily: "monospace",
            }}
          >
            {preview}
            {entry.bytes.length > 3000 && "\n…(truncated)"}
          </pre>
        </Box>
      )}
    </Box>
  );
}

// ─── App ──────────────────────────────────────────────────────────────────────

const App = () => {
  // Worker
  const workerRef = useRef<Worker | null>(null);
  const pendingRunsRef = useRef(
    new Map<
      number,
      {
        resolve: (v: { exitCode: number; results: ResultEntry[] }) => void;
        reject: (e: Error) => void;
      }
    >()
  );
  const nextRunIdRef = useRef(1);

  // UI state
  const [fastaFile, setFastaFile] = useState<File | null>(null);
  const [cmdOverride, setCmdOverride] = useState("");
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [loadingProgress, setLoadingProgress] = useState<number | null>(null);
  const [loadingMessage, setLoadingMessage] = useState("");
  const [executionTime, setExecutionTime] = useState(0);
  const [stdout, setStdout] = useState("");
  const [stderr, setStderr] = useState("");
  const [results, setResults] = useState<ResultEntry[]>([]);
  const [exitCode, setExitCode] = useState<number | null>(null);
  const [showLog, setShowLog] = useState(true);  // open by default

  // Output buffering
  const stdoutBufRef = useRef<string[]>([]);
  const stderrBufRef = useRef<string[]>([]);
  const flushTimerRef = useRef<number | null>(null);
  const execTimerRef = useRef<number | null>(null);
  const startTimeRef = useRef<number | null>(null);
  const logBottomRef = useRef<HTMLDivElement | null>(null);

  // ── Worker ────────────────────────────────────────────────────────────────

  useEffect(() => {
    const worker = new Worker(new URL("../worker/runner.ts", import.meta.url), {
      type: "module",
    });

    worker.onmessage = (ev: MessageEvent<WorkerMessage>) => {
      const m = ev.data;
      if (!m?.type) return;

      if (m.type === "stdout") { bufferLine("stdout", m.line); return; }
      if (m.type === "stderr") { bufferLine("stderr", m.line); return; }

      if (m.type === "loading") {
        setLoadingProgress(m.progress);
        setLoadingMessage(m.message);
        return;
      }

      const pending = pendingRunsRef.current.get(m.runId);
      if (!pending) return;

      if (m.type === "done") {
        pendingRunsRef.current.delete(m.runId);
        flushBuffers();
        pending.resolve({ exitCode: m.exitCode, results: m.results });
        return;
      }
      if (m.type === "error") {
        pendingRunsRef.current.delete(m.runId);
        flushBuffers();
        pending.reject(new Error(m.message));
      }
    };

    workerRef.current = worker;
    return () => worker.terminate();
  }, []);

  // ── Execution timer ───────────────────────────────────────────────────────

  useEffect(() => {
    if (isRunning) {
      startTimeRef.current = Date.now();
      execTimerRef.current = window.setInterval(() => {
        if (startTimeRef.current)
          setExecutionTime(Math.floor((Date.now() - startTimeRef.current) / 1000));
      }, 1000);
    } else {
      if (execTimerRef.current) { clearInterval(execTimerRef.current); execTimerRef.current = null; }
    }
    return () => { if (execTimerRef.current) clearInterval(execTimerRef.current); };
  }, [isRunning]);

  // ── Buffered output ───────────────────────────────────────────────────────

  const bufferLine = (type: "stdout" | "stderr", line: string) => {
    (type === "stdout" ? stdoutBufRef : stderrBufRef).current.push(line);
    if (!flushTimerRef.current)
      flushTimerRef.current = window.setTimeout(flushBuffers, 16);
  };

  const flushBuffers = () => {
    if (stdoutBufRef.current.length) {
      const chunk = stdoutBufRef.current.join("\n");
      setStdout((p) => p + chunk + "\n");
      stdoutBufRef.current = [];
    }
    if (stderrBufRef.current.length) {
      const chunk = stderrBufRef.current.join("\n");
      setStderr((p) => p + chunk + "\n");
      stderrBufRef.current = [];
    }
    flushTimerRef.current = null;
    // Scroll log to bottom after flushing
    requestAnimationFrame(() => {
      logBottomRef.current?.scrollIntoView({ behavior: "smooth", block: "nearest" });
    });
  };

  // ── Run ───────────────────────────────────────────────────────────────────

  const handleRun = useCallback(async () => {
    if (!fastaFile) return;

    // Reset state
    setStdout("");
    setStderr("");
    setResults([]);
    setExitCode(null);
    setExecutionTime(0);
    setLoadingProgress(null);
    setLoadingMessage("");
    stdoutBufRef.current = [];
    stderrBufRef.current = [];

    setIsRunning(true);
    console.log(`[CFEIntact] Starting run. WASM URL: ${WASM_URL}`);

    try {
      const fastaBytes = await fastaFile.arrayBuffer();
      const runId = nextRunIdRef.current++;

      // Override the container entrypoint to /bin/sh using c2w's -entrypoint flag,
      // then run the user's command followed by a cp to collect outputs.
      // The cp step is hidden from the user; only the cfeintact exit code is surfaced.
      const userCmd = (cmdOverride.trim() || DEFAULT_CMD).trim();
      // Run cfeintact in the container's internal /tmp (NOT the 9P-backed /w).
      // That avoids EOPNOTSUPP when Python truncates blast.csv on the 2nd subtype
      // iteration. Afterwards, cp the results to the WASI-preopened /mnt/out/.
      const shellCmd =
        `mkdir -p ${CONTAINER_WORK_DIR} && cd ${CONTAINER_WORK_DIR} && ${userCmd}; rc=$?; cp -- ${CONTAINER_WORK_DIR}/*.csv ${CONTAINER_WORK_DIR}/*.fasta ${OUTPUT_MOUNT}/ 2>/dev/null; exit $rc`;
      const request: RunRequest = {
        runId,
        wasmUrl: WASM_URL,
        mountPaths: [INPUT_MOUNT, OUTPUT_MOUNT],
        fileMappings: [{ path: INPUT_PATH, bytes: fastaBytes }],
        wasiArgs: ["-entrypoint", "/bin/sh"], // overrides the container entrypoint
        cmdArgs: ["-c", shellCmd],             // CMD passed to /bin/sh
      };

      const p = new Promise<{ exitCode: number; results: ResultEntry[] }>(
        (resolve, reject) => pendingRunsRef.current.set(runId, { resolve, reject })
      );

      workerRef.current?.postMessage(request, [fastaBytes]);

      const { exitCode: code, results: files } = await p;

      setExitCode(code);
      // Collect only files written to the output mount by the cp step.
      setResults(files.filter((f) => f.path.startsWith(OUTPUT_MOUNT + "/")));

      if (code === 0) {
        toaster.success({ title: "Analysis complete" });
      } else {
        toaster.warning({ title: "Analysis finished", description: `Exit code ${code}` });
      }
    } catch (e: unknown) {
      const msg = e instanceof Error ? e.message : String(e);
      console.error(`[CFEIntact] Run failed. WASM URL was: ${WASM_URL}\nError: ${msg}`);
      toaster.error({
        title: "Error",
        description: `${msg} (WASM: ${WASM_URL})`,
      });
    } finally {
      setIsRunning(false);
      setLoadingProgress(null);
      setLoadingMessage("");
    }
  }, [fastaFile, cmdOverride]);

  // ── Derived ───────────────────────────────────────────────────────────────

  const csvResults = results.filter((r) => r.path.toLowerCase().endsWith(".csv"));
  const otherResults = results.filter((r) => !r.path.toLowerCase().endsWith(".csv"));
  const hasLog = isRunning || stdout.trim() || stderr.trim();
  const hasResults = results.length > 0;

  // ── Render ────────────────────────────────────────────────────────────────

  return (
    <Box minH="100vh" bg="gray.50">

      {/* Header */}
      <Box
        style={{ background: "linear-gradient(135deg, #234e52 0%, #2c7a7b 100%)" }}
        color="white"
        py={10}
        px={4}
      >
        <Container maxW="container.md">
          <VStack gap={2} align="start">
            <HStack gap={3} align="center">
              <Text fontSize="3xl" lineHeight={1}>🧬</Text>
              <Heading size="2xl" letterSpacing="tight">CFEIntact</Heading>
            </HStack>
            <Text fontSize="md" style={{ opacity: 0.85 }}>
              HIV-1 Genome Intactness Analysis
            </Text>
            <HStack gap={2} mt={2} flexWrap="wrap">
              <Box
                bg="rgba(255,255,255,0.15)"
                borderRadius="full"
                px={3}
                py={1}
                fontSize="xs"
                backdropFilter="blur(4px)"
              >
                Browser-powered · WebAssembly
              </Box>
              <Box
                bg="rgba(255,255,255,0.15)"
                borderRadius="full"
                px={3}
                py={1}
                fontSize="xs"
                backdropFilter="blur(4px)"
              >
                No data leaves your device
              </Box>
            </HStack>
          </VStack>
        </Container>
      </Box>

      {/* Content */}
      <Container maxW="container.md" py={8}>
        <VStack gap={6} align="stretch">

          {/* Upload */}
          <Box bg="white" borderRadius="2xl" p={6} shadow="sm" borderWidth={1} borderColor="gray.200">
            <Text fontWeight="bold" fontSize="lg" mb={4} color="gray.700">
              Input FASTA File
            </Text>
            <DropZone file={fastaFile} onFile={setFastaFile} disabled={isRunning} />
          </Box>

          {/* Advanced options */}
          <Box bg="white" borderRadius="2xl" shadow="sm" borderWidth={1} borderColor="gray.200" overflow="hidden">
            <Button
              variant="ghost"
              size="sm"
              w="100%"
              justifyContent="space-between"
              onClick={() => setShowAdvanced((v) => !v)}
              px={6}
              py={4}
              h="auto"
              borderRadius={0}
              color="gray.500"
            >
              <Text fontWeight="medium">Advanced options</Text>
              <Text>{showAdvanced ? "▲" : "▼"}</Text>
            </Button>

            {showAdvanced && (
              <Box px={6} pb={5} borderTopWidth={1} borderColor="gray.100">
                <Text fontSize="sm" fontWeight="medium" color="gray.600" mt={4} mb={2}>
                  Command (executed inside the container)
                </Text>
                <textarea
                  value={cmdOverride || DEFAULT_CMD}
                  onChange={(e) => setCmdOverride(e.target.value)}
                  rows={2}
                  style={{
                    width: "100%",
                    fontFamily: "monospace",
                    fontSize: "0.8rem",
                    padding: "8px 12px",
                    border: "1px solid #CBD5E0",
                    borderRadius: "8px",
                    resize: "vertical",
                    background: "#F7FAFC",
                    outline: "none",
                    boxSizing: "border-box",
                  }}
                />
                <Text fontSize="xs" color="gray.400" mt={1}>
                  Command run inside the container. Input at <code>{INPUT_PATH}</code>.
                  CFEIntact outputs are automatically collected from the container and made available for download when the run completes.
                </Text>

                {/* Diagnostic info */}
                <Box mt={4} p={3} bg="gray.100" borderRadius="md" fontSize="xs" fontFamily="monospace">
                  <Text fontWeight="semibold" color="gray.600" mb={1}>Diagnostic info</Text>
                  <Text color="gray.600">Page: {window.location.href}</Text>
                  <Text color="gray.600">
                    WASM:{" "}
                    <a href={WASM_URL} target="_blank" rel="noreferrer" style={{ color: "#2b6cb0", wordBreak: "break-all" }}>
                      {WASM_URL}
                    </a>
                  </Text>
                </Box>
              </Box>
            )}
          </Box>

          {/* Run button */}
          <Button
            colorScheme="teal"
            size="lg"
            borderRadius="xl"
            h={14}
            fontSize="md"
            fontWeight="bold"
            onClick={handleRun}
            disabled={!fastaFile || isRunning}
            loading={isRunning}
            loadingText={
              loadingProgress !== null
                ? loadingMessage || "Loading…"
                : `Running… ${executionTime}s`
            }
            shadow="md"
          >
            Run Analysis
          </Button>

          {/* Progress */}
          {isRunning && (
            <Box bg="white" borderRadius="2xl" p={5} shadow="sm" borderWidth={1} borderColor="gray.200">
              {loadingProgress !== null ? (
                <ProgressBar progress={loadingProgress} label={loadingMessage || "Preparing…"} />
              ) : (
                <ProgressBar progress={-1} label={`Executing… ${executionTime}s`} />
              )}
            </Box>
          )}

          {/* Results */}
          {hasResults && (
            <Box bg="white" borderRadius="2xl" p={6} shadow="sm" borderWidth={1} borderColor="gray.200">
              <HStack justify="space-between" mb={4}>
                <Text fontWeight="bold" fontSize="lg" color="gray.700">Results</Text>
                {exitCode !== null && (
                  <Badge
                    colorScheme={exitCode === 0 ? "green" : "orange"}
                    borderRadius="full"
                    px={3}
                    py={1}
                  >
                    Exit {exitCode}
                  </Badge>
                )}
              </HStack>

              {csvResults.length > 0 ? (
                <VStack gap={3} align="stretch" mb={otherResults.length > 0 ? 4 : 0}>
                  {csvResults.map((r) => <ResultCard key={r.path} entry={r} />)}
                </VStack>
              ) : (
                <Box
                  bg="orange.50"
                  borderRadius="lg"
                  p={4}
                  borderWidth={1}
                  borderColor="orange.200"
                  mb={otherResults.length > 0 ? 4 : 0}
                >
                  <Text fontSize="sm" color="orange.700">
                    No CSV files were produced. Check the execution log below for details.
                  </Text>
                </Box>
              )}

              {otherResults.length > 0 && (
                <VStack gap={3} align="stretch">
                  <Text fontSize="sm" fontWeight="medium" color="gray.500">Other output files</Text>
                  {otherResults.map((r) => <ResultCard key={r.path} entry={r} />)}
                </VStack>
              )}
            </Box>
          )}

          {/* Log */}
          {hasLog && (
            <Box bg="white" borderRadius="2xl" shadow="sm" borderWidth={1} borderColor="gray.200" overflow="hidden">
              <Button
                variant="ghost"
                size="sm"
                w="100%"
                justifyContent="space-between"
                onClick={() => setShowLog((v) => !v)}
                px={6}
                py={4}
                h="auto"
                borderRadius={0}
                color="gray.500"
              >
                <HStack gap={2}>
                  <Text fontWeight="medium">Execution log</Text>
                  {isRunning && (
                    <Box
                      w={2} h={2} borderRadius="full" bg="teal.400"
                      style={{ animation: "pulse-dot 1s ease-in-out infinite" }}
                    />
                  )}
                </HStack>
                <Text>{showLog ? "▲" : "▼"}</Text>
              </Button>

              {showLog && (
                <Box px={6} pb={5} borderTopWidth={1} borderColor="gray.100">
                  {isRunning && !stdout.trim() && !stderr.trim() && (
                    <Box mt={4}>
                      <Text fontSize="xs" color="gray.400" fontStyle="italic">Waiting for output…</Text>
                    </Box>
                  )}
                  {stdout.trim() && (
                    <Box mt={4} mb={stderr.trim() ? 3 : 0}>
                      <Text fontSize="xs" fontWeight="bold" color="green.600" mb={1}>stdout</Text>
                      <Box bg="gray.50" borderRadius="md" p={3} maxH="320px" overflowY="auto" borderWidth={1} borderColor="gray.200">
                        <pre style={{ margin: 0, fontSize: "0.75rem", whiteSpace: "pre-wrap", fontFamily: "monospace" }}>
                          {stdout}
                        </pre>
                      </Box>
                    </Box>
                  )}
                  {stderr.trim() && (
                    <Box mt={stdout.trim() ? 0 : 4}>
                      <Text fontSize="xs" fontWeight="bold" color="red.500" mb={1}>stderr</Text>
                      <Box bg="red.50" borderRadius="md" p={3} maxH="320px" overflowY="auto" borderWidth={1} borderColor="red.200">
                        <pre style={{ margin: 0, fontSize: "0.75rem", whiteSpace: "pre-wrap", fontFamily: "monospace" }}>
                          {stderr}
                        </pre>
                      </Box>
                    </Box>
                  )}
                  {/* Invisible sentinel – scrolled into view after each flush */}
                  <div ref={logBottomRef} style={{ height: 1 }} />
                </Box>
              )}
            </Box>
          )}

          {/* Footer */}
          <Text fontSize="xs" color="gray.400" textAlign="center" pb={4}>
            Powered by container2wasm · Analysis runs entirely in your browser
          </Text>

        </VStack>
      </Container>
    </Box>
  );
};

export default App;
