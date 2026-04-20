/**
 * Output rendering and download functionality
 */
import JSZip from "jszip";
import type { ResultEntry } from "../shared/types";

/**
 * Render output files with download links and ZIP capabilities
 */
export async function renderOutFiles(fileList: ResultEntry[], container: HTMLElement): Promise<void> {
  if (!fileList.length) {
    container.textContent = "(none)";
    return;
  }

  // Group files by mount point (top-level directory)
  const mountGroups = new Map<string, ResultEntry[]>();

  for (const entry of fileList) {
    const parts = entry.path.split("/").filter((p) => p);
    const mountPoint = parts.length > 0 ? `/${parts[0]}` : "/";

    if (!mountGroups.has(mountPoint)) {
      mountGroups.set(mountPoint, []);
    }
    mountGroups.get(mountPoint)!.push(entry);
  }

  const ul = document.createElement("ul");

  for (const [mountPoint, entries] of mountGroups) {
    const mountLi = document.createElement("li");
    mountLi.style.marginBottom = "16px";

    const mountHeader = document.createElement("strong");
    mountHeader.textContent = `${mountPoint}/ `;
    mountLi.appendChild(mountHeader);

    // Add "Download all as ZIP" button for the mount
    const zipBtn = document.createElement("button");
    zipBtn.textContent = "Download all as ZIP";
    zipBtn.style.padding = "4px 8px";
    zipBtn.style.fontSize = "12px";
    zipBtn.style.marginLeft = "8px";
    zipBtn.onclick = async () => {
      zipBtn.disabled = true;
      zipBtn.textContent = "Creating ZIP...";
      try {
        const zipBlob = await createZipFromEntries(entries, mountPoint);
        const a = document.createElement("a");
        a.href = URL.createObjectURL(zipBlob);
        a.download = `${mountPoint.replace(/\//g, "_")}.zip`;
        a.click();
        zipBtn.textContent = "Download all as ZIP";
        zipBtn.disabled = false;
      } catch (e: unknown) {
        const message = e instanceof Error ? e.message : String(e);
        alert(`Error creating zip: ${message}`);
        zipBtn.textContent = "Download all as ZIP";
        zipBtn.disabled = false;
      }
    };
    mountLi.appendChild(zipBtn);

    // List individual files
    const fileUl = document.createElement("ul");
    fileUl.style.marginTop = "8px";

    for (const entry of entries) {
      const fileLi = document.createElement("li");
      const a = document.createElement("a");
      const displayPath = entry.path.startsWith(mountPoint) 
        ? entry.path.slice(mountPoint.length + 1) || entry.path 
        : entry.path;
      a.textContent = `${displayPath} (${entry.bytes.length} bytes)`;
      // Create a new Uint8Array to ensure ArrayBuffer (not SharedArrayBuffer) backing
      const blobBytes = new Uint8Array(entry.bytes);
      a.href = URL.createObjectURL(new Blob([blobBytes]));
      a.download = displayPath.split("/").pop() || "file";
      fileLi.appendChild(a);
      fileUl.appendChild(fileLi);
    }

    mountLi.appendChild(fileUl);
    ul.appendChild(mountLi);
  }

  container.innerHTML = "";
  container.appendChild(ul);
}

/**
 * Create a ZIP file from a list of result entries
 */
async function createZipFromEntries(entries: ResultEntry[], basePrefix: string): Promise<Blob> {
  const zip = new JSZip();

  for (const entry of entries) {
    // Remove the base prefix from the path for the ZIP structure
    let zipPath = entry.path;
    if (zipPath.startsWith(basePrefix + "/")) {
      zipPath = zipPath.slice(basePrefix.length + 1);
    } else if (zipPath === basePrefix) {
      zipPath = "root";
    }

    zip.file(zipPath, entry.bytes);
  }

  return await zip.generateAsync({ type: "blob" });
}
