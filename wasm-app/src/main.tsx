/**
 * React app entry point
 */
import React from "react";
import ReactDOM from "react-dom/client";
import { ChakraProvider, defaultSystem } from "@chakra-ui/react";
import { Toaster } from "./ui/toaster";
import App from "./components/App";

ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <ChakraProvider value={defaultSystem}>
      <App />
      <Toaster />
    </ChakraProvider>
  </React.StrictMode>
);
