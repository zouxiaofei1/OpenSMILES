import React from "react";
import { createRoot } from "react-dom/client";
import "ketcher-react/dist/index.css";
import { KetcherApp } from "./KetcherApp.jsx";

createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <KetcherApp />
  </React.StrictMode>
);
