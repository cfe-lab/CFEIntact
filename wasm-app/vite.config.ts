import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  base: "./",
  plugins: [react()],
  build: {
    outDir: "dist",
    emptyOutDir: true,
    rollupOptions: {
      output: {
        // Ensure stable chunk naming
        manualChunks: undefined,
      },
    },
  },
  worker: {
    format: "es", // Use ES modules for workers
  },
  server: {
    port: 3000,
    open: true,
  },
});
