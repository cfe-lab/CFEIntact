import { useState } from "react";
import {
  Box,
  Text,
  VStack,
  HStack,
  Button,
  Link,
  Collapsible,
} from "@chakra-ui/react";
import JSZip from "jszip";
import type { ResultEntry } from "../shared/types";

interface FileResultsProps {
  results: ResultEntry[];
}

const FileResults = ({ results }: FileResultsProps) => {
  if (results.length === 0) {
    return (
      <Box>
        <Text fontWeight="bold" mb={2}>Output Files</Text>
        <Text color="gray.500" fontStyle="italic">
          (none)
        </Text>
      </Box>
    );
  }

  // Group files by mount point
  const mountGroups = new Map<string, ResultEntry[]>();

  for (const entry of results) {
    const parts = entry.path.split("/").filter((p) => p);
    const mountPoint = parts.length > 0 ? `/${parts[0]}` : "/";

    if (!mountGroups.has(mountPoint)) {
      mountGroups.set(mountPoint, []);
    }
    mountGroups.get(mountPoint)!.push(entry);
  }

  return (
    <Box>
      <Text fontWeight="bold" mb={2}>Output Files</Text>
      <VStack gap={4} align="stretch">
        {Array.from(mountGroups.entries()).map(([mountPoint, entries]) => (
          <MountPointGroup key={mountPoint} mountPoint={mountPoint} entries={entries} />
        ))}
      </VStack>
    </Box>
  );
};

interface MountPointGroupProps {
  mountPoint: string;
  entries: ResultEntry[];
}

const MountPointGroup = ({ mountPoint, entries }: MountPointGroupProps) => {
  const [isOpen, setIsOpen] = useState(true);
  const [isCreatingZip, setIsCreatingZip] = useState(false);

  const handleDownloadZip = async () => {
    setIsCreatingZip(true);
    try {
      const zip = new JSZip();

      for (const entry of entries) {
        let zipPath = entry.path;
        if (zipPath.startsWith(mountPoint + "/")) {
          zipPath = zipPath.slice(mountPoint.length + 1);
        } else if (zipPath === mountPoint) {
          zipPath = "root";
        }
        zip.file(zipPath, entry.bytes);
      }

      const blob = await zip.generateAsync({ type: "blob" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `${mountPoint.replace(/\//g, "_")}.zip`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (e: unknown) {
      const message = e instanceof Error ? e.message : String(e);
      alert(`Error creating ZIP: ${message}`);
    } finally {
      setIsCreatingZip(false);
    }
  };

  return (
    <Box borderWidth={1} borderRadius="md" p={4} bg="gray.50">
      <HStack justify="space-between" mb={2}>
        <HStack cursor="pointer" onClick={() => setIsOpen(!isOpen)}>
          <Text fontSize="lg">{isOpen ? "▼" : "▶"}</Text>
          <Text fontWeight="bold" fontSize="lg">
            {mountPoint}/
          </Text>
          <Text color="gray.600" fontSize="sm">
            ({entries.length} {entries.length === 1 ? "file" : "files"})
          </Text>
        </HStack>
        <Button
          size="sm"
          colorScheme="blue"
          onClick={handleDownloadZip}
          loading={isCreatingZip}
          loadingText="Creating..."
        >
          ⬇ Download as ZIP
        </Button>
      </HStack>

      <Collapsible.Root open={isOpen}>
        <Collapsible.Content>
          <VStack align="stretch" gap={1} ml={6} mt={2}>
            {entries.map((entry: ResultEntry, idx: number) => {
              const displayPath =
                entry.path.startsWith(mountPoint)
                  ? entry.path.slice(mountPoint.length + 1) || entry.path
                  : entry.path;
              const fileName = displayPath.split("/").pop() || "file";
              const blobBytes = new Uint8Array(entry.bytes);
              const url = URL.createObjectURL(new Blob([blobBytes]));

              return (
                <HStack key={idx} gap={2}>
                  <Link href={url} download={fileName} color="blue.600" fontFamily="monospace" fontSize="sm">
                    {displayPath}
                  </Link>
                  <Text color="gray.500" fontSize="xs">
                    ({entry.bytes.length} bytes)
                  </Text>
                </HStack>
              );
            })}
          </VStack>
        </Collapsible.Content>
      </Collapsible.Root>
    </Box>
  );
};

export default FileResults;
