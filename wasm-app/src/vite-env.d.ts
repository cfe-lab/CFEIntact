/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_WASM_VERSION?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
