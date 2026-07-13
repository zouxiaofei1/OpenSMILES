import path from "node:path";
import { fileURLToPath } from "node:url";
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

/**
 * ketcher-core does require('raphael') at module top-level.
 * Map that to global Raphael (loaded via public/raphael.min.js script tag).
 */
function raphaelRequirePlugin() {
  return {
    name: "raphael-require-to-window",
    transform(code, id) {
      if (!id.includes("node_modules")) return null;
      if (!code.includes('require("raphael")') && !code.includes("require('raphael')")) {
        return null;
      }
      const next = code
        .replace(
          /require\(["']raphael["']\)/g,
          "(globalThis.Raphael || window.Raphael)"
        )
        .replace(
          /typeof\s+(\w+)\s*===\s*["']function["']\s*\?\s*\1\s*:\s*\1\[["']default["']\]/g,
          "(typeof $1==='function'?$1:($1&&$1.default))"
        )
        .replace(
          /typeof\s+(\w+)\s*==\s*["']function["']\s*\?\s*\1\s*:\s*\1\.default/g,
          "(typeof $1==='function'?$1:($1&&$1.default))"
        );
      return { code: next, map: null };
    },
  };
}

export default defineConfig({
  plugins: [react(), raphaelRequirePlugin()],
  base: "/vendor/ketcher/",
  publicDir: "public",
  resolve: {
    alias: {
      // Avoid bundling a second Raphael via import paths
      raphael: path.resolve(__dirname, "node_modules/raphael/raphael.min.js"),
    },
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
    assetsDir: "assets",
    commonjsOptions: {
      transformMixedEsModules: true,
      include: [/node_modules/],
    },
    chunkSizeWarningLimit: 30000,
  },
  define: {
    "process.env.NODE_ENV": JSON.stringify("production"),
  },
});
