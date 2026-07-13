import React from "react";
import { createRoot } from "react-dom/client";
import "ketcher-react/dist/index.css";
import { KetcherApp } from "./KetcherApp.jsx";

// Raphael is provided by public/raphael.min.js (script tag before this module).
if (typeof window !== "undefined" && !window.Raphael) {
  console.error("Raphael global missing; Ketcher shell will fail to init");
}

createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <KetcherApp />
  </React.StrictMode>
);
